#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <stdexcept>

#include "scene.h"
#include "json.hpp"
#include "utils.h"
#include "overloads.h"
#include "happly.h"
#include "bvh.h"

using json = nlohmann::json;

#define VERBOSE 1

void verbose(const std::string& message) {
    if (VERBOSE) {
        std::cout << "[SCENE] " << message << std::endl;
    }
}

template <typename T> 
T parseSingleValue(const std::string& value) {
    T parsedValue;
    std::istringstream stream(value);
    stream >> parsedValue;
    stream.clear();
    return parsedValue;
}

template <typename T> 
T parsePair(const std::string& value) {
    T parsedValue;
    std::istringstream stream(value);
    stream >> parsedValue.x >> parsedValue.y;
    stream.clear();
    return parsedValue;
}

template <typename T> 
T parseTriplet(const std::string& value) {
    T triplet;
    std::istringstream stream(value);
    stream >> triplet.x >> triplet.y >> triplet.z;
    stream.clear();
    return triplet;
}

template <typename T> 
T parseQuad(const std::string& value) {
    T quad;
    std::istringstream stream(value);
    stream >> quad.x >> quad.y >> quad.z >> quad.w;
    stream.clear();
    return quad;
}


void scene::Scene::loadSceneFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error: The json file cannot be loaded.");
    }

    // Extract base directory from filename
    size_t lastSlash = filename.find_last_of("/\\");
    if (lastSlash != std::string::npos) {
        this->baseDirectory = filename.substr(0, lastSlash + 1);
    } else {
        this->baseDirectory = "";
    }

    json j;
    file >> j;
    file.close();

    if (!j.contains("Scene")) {
        throw std::runtime_error("Error: Scene root is not found.");
    }

    auto scene = j["Scene"];
    if (VERBOSE) {
        verbose("================================================");
        verbose("Parsing Scene File: " + filename);
        verbose("Base Directory: " + this->baseDirectory);
        verbose("================================================");
    }

    if (scene.contains("BackgroundColor") && !scene["BackgroundColor"].is_null()) {
        std::string bgColor = scene["BackgroundColor"].get<std::string>();
        this->backgroundColor = parseTriplet<VectorFloatTriplet>(bgColor);

        verbose("[+] BackgroundColor Parsed: " + std::to_string(this->backgroundColor.x) + " " + std::to_string(this->backgroundColor.y) + " " + std::to_string(this->backgroundColor.z));
    } else {
        this->backgroundColor.x = this->backgroundColor.y = this->backgroundColor.z = 0;
        verbose("[!] Skipping BackgroundColor Parsing. Reason: Not found in the scene file. Assigning default value: " + std::to_string(this->backgroundColor.x) + " " + std::to_string(this->backgroundColor.y) + " " + std::to_string(this->backgroundColor.z));
    }

    if (scene.contains("MaxRecursionDepth") && !scene["MaxRecursionDepth"].is_null()) {
        std::string maxRecursionDepth = scene["MaxRecursionDepth"].get<std::string>();
        this->maxRecursionDepth = parseSingleValue<int>(maxRecursionDepth);
        verbose("[+] MaxRecursionDepth parsed: " + std::to_string(this->maxRecursionDepth)); 
    } else {
        this->maxRecursionDepth = 0;
        verbose("[!] Skipping MaxRecursionDepth parsing. Reason: Not found in the scene file " + std::to_string(this->maxRecursionDepth)); 
    }

    if (scene.contains("ShadowRayEpsilon") && !scene["ShadowRayEpsilon"].is_null()) {
        std::string shadowRayEpsilon = scene["ShadowRayEpsilon"].get<std::string>();
        this->shadowRayEpsilon = parseSingleValue<float>(shadowRayEpsilon);
        verbose("[+] ShadowRayEpsilon Parsed: " + std::to_string(this->shadowRayEpsilon));
    } else {
        this->shadowRayEpsilon = 0.001f;
        verbose("[!] Skipping ShadowRayEpsilon Parsing. Reason: Not found in the scene file. Assigning default value: " + std::to_string(shadowRayEpsilon));
    }

    if (scene.contains("IntersectionTestEpsilon") && !scene["IntersectionTestEpsilon"].is_null()) {
        std::string intersectionTestEpsilonStr = scene["IntersectionTestEpsilon"].get<std::string>();
        this->intersectionTestEpsilon = parseSingleValue<float>(intersectionTestEpsilonStr);
        verbose("[+] IntersectionTestEpsilon Parsed: " + std::to_string(this->intersectionTestEpsilon));
    } else {
        this->intersectionTestEpsilon = 1e-6f;
        verbose("[!] Skipping IntersectionTestEpsilon Parsing. Reason: Not found in the scene file. Assigning default value: " + std::to_string(this->intersectionTestEpsilon));
    }

    if (scene.contains("Cameras") && !scene["Cameras"].is_null()) {
        auto cameras = scene["Cameras"];
        auto cameraDataArray = cameras["Camera"];
        if (cameraDataArray.is_array()) {
            verbose("parsing camera array");
            for (auto cameraData : cameraDataArray) {
                scene::Camera newCamera = parseCamera(cameraData);
                if (newCamera._id != 0) {
                    this->cameras.push_back(newCamera);
                    verbose("[+] Camera Parsed: " + std::to_string(newCamera._id));
                }
            }
        } else {
            scene::Camera newCamera = parseCamera(cameraDataArray);
            if (newCamera._id != 0) {
                this->cameras.push_back(newCamera);
                verbose("[+] Camera Parsed: " + std::to_string(newCamera._id));
            }
        }
    } else {
        verbose("[!] Skipping Cameras Parsing. Reason: Not found in the scene file. Assigning default value: 0");
    }

    if (scene.contains("Lights") && !scene["Lights"].is_null()) {
        auto lights = scene["Lights"];
        if (lights.contains("AmbientLight") && !lights["AmbientLight"].is_null()) {
            std::string ambientLight = lights["AmbientLight"].get<std::string>();
            this->ambientLight.intensity = parseTriplet<VectorFloatTriplet>(ambientLight);
            verbose("[+] AmbientLight Parsed: " + std::to_string(this->ambientLight.intensity.x) + " " + std::to_string(this->ambientLight.intensity.y) + " " + std::to_string(this->ambientLight.intensity.z));
        } else {
            this->ambientLight.intensity.x = this->ambientLight.intensity.y = this->ambientLight.intensity.z = 0;
            verbose("[!] Skipping AmbientLight Parsing. Reason: Not found in the scene file. Assigning default value: " + std::to_string(this->ambientLight.intensity.x) + " " + std::to_string(this->ambientLight.intensity.y) + " " + std::to_string(this->ambientLight.intensity.z));
        }

        if (lights.contains("PointLight") && !lights["PointLight"].is_null()) {
            auto pointLightArray = lights["PointLight"];
            if (pointLightArray.is_array()) {
                for (auto pointLightData : pointLightArray) {
                    scene::PointLight newPointLight = parsePointLight(pointLightData);
                    if (newPointLight._id != 0) {
                        this->pointLights.push_back(newPointLight);
                        verbose("[+] PointLight Parsed: " + std::to_string(newPointLight._id));
                    }
                }
            } else {
                scene::PointLight newPointLight = parsePointLight(pointLightArray);
                if (newPointLight._id != 0) {
                    this->pointLights.push_back(newPointLight);
                    verbose("[+] PointLight Parsed: " + std::to_string(newPointLight._id));
                }
            }
        } else {
            verbose("[!] Skipping PointLight Parsing. Reason: Not found in the scene file. Assigning default value: 0");
        }
    } else {
        verbose("[!] Skipping Lights Parsing. Reason: Not found in the scene file. Assigning default value: 0");
    }

    if (scene.contains("Materials") && !scene["Materials"].is_null()) {
        auto materials = scene["Materials"];
        auto materialDataArray = materials["Material"];
        if (materialDataArray.is_array()) {
            for (auto materialData : materialDataArray) {
                scene::Material newMaterial = parseMaterial(materialData);
                if (newMaterial._id != 0) {
                    this->materials.push_back(newMaterial);
                    this->materialIdToIndex[newMaterial._id] = this->materials.size() - 1;
                    verbose("[+] Material Parsed: " + std::to_string(newMaterial._id));
                } else {
                    verbose("[!] Skipping Material Parsing. Reason: Material ID is 0");
                }
            }
        } else {
            scene::Material newMaterial = parseMaterial(materialDataArray);
            if (newMaterial._id != 0) {
                this->materials.push_back(newMaterial);
                this->materialIdToIndex[newMaterial._id] = this->materials.size() - 1;
                verbose("[+] Material Parsed: " + std::to_string(newMaterial._id));
            }
        }
    } else {
        verbose("[!] Skipping Materials Parsing. Reason: Not found in the scene file. Assigning default value: 0");
    }

    if (scene.contains("VertexData") && !scene["VertexData"].is_null()) {
        auto vertexData = scene["VertexData"];
        auto vertexDataArray = vertexData["_data"];
        if (vertexDataArray.is_array()) {
            for (auto vertexData : vertexDataArray) {
                try {
                    std::vector<VectorFloatTriplet> vertices = parseVertex(vertexData);
                    this->vertices.insert(this->vertices.end(), vertices.begin(), vertices.end());
                    verbose("[+] Vertex Parsed: " + std::to_string(vertices.size()));
                } catch (const std::exception& e) {
                    verbose("[!] Skipping Vertex Parsing. Reason: " + std::string(e.what()));
                }   
            }
        } else {
            std::vector<VectorFloatTriplet> vertices = parseVertex(vertexDataArray);
            this->vertices.insert(this->vertices.end(), vertices.begin(), vertices.end());
            verbose("[+] Vertex Parsed: " + std::to_string(vertices.size()));
        }
    }

    if (scene.contains("Objects") && !scene["Objects"].is_null()) {
        auto objects = scene["Objects"];
        auto meshes = objects["Mesh"];
        auto triangles = objects["Triangle"];
        auto spheres = objects["Sphere"];
        auto planes = objects["Plane"];

        this->meshes = parseObjects<scene::Mesh>(meshes);
        this->triangles = parseObjects<scene::Triangle>(triangles);
        this->spheres = parseObjects<scene::Sphere>(spheres);
        this->planes = parseObjects<scene::Plane>(planes);

        verbose("[+] Meshes Parsed: " + std::to_string(this->meshes.size()));
        verbose("[+] Triangles Parsed: " + std::to_string(this->triangles.size()));
        verbose("[+] Spheres Parsed: " + std::to_string(this->spheres.size()));
        verbose("[+] Planes Parsed: " + std::to_string(this->planes.size()));
    }

    verbose("================================================");
    verbose("Scene File Parsed Successfully");
    verbose("================================================");
}

scene::Camera scene::parseCamera(const json& cameraData) {
    char cameraType = 0;
    if (cameraData.contains("_type") && !cameraData["_type"].is_null()) {
        if (cameraData["_type"].get<std::string>() == "lookAt") {
            cameraType = 1;
        } 
    }
    scene::Camera newCamera = scene::Camera();
    switch (cameraType) {
        case 1: {
            newCamera._id = parseSingleValue<unsigned int>(cameraData["_id"]);
            newCamera.position = parseTriplet<VectorFloatTriplet>(cameraData["Position"]);
            newCamera.up = parseTriplet<VectorFloatTriplet>(cameraData["Up"]);
            newCamera.nearDistance = parseSingleValue<float>(cameraData["NearDistance"]);
            newCamera.imageResolution = parsePair<VectorIntPair>(cameraData["ImageResolution"]);
            newCamera.imageName = cameraData["ImageName"].get<std::string>();
            // Calculate the gaze vector, nearplane
            // Gaze: GazePoint - Position
            VectorFloatTriplet gazePoint = parseTriplet<VectorFloatTriplet>(cameraData["GazePoint"]);
            VectorFloatTriplet gaze = normalize(gazePoint - newCamera.position);
            newCamera.gaze = gaze;

            // NearPlane: 
            float fovY = parseSingleValue<float>(cameraData["FovY"]);
            float fovX = fovY * newCamera.imageResolution.x / newCamera.imageResolution.y;
            float nearPlaneWidth = 2 * tan(fovX / 2) * newCamera.nearDistance;
            float nearPlaneHeight = 2 * tan(fovY / 2) * newCamera.nearDistance;
            newCamera.nearPlane = VectorFloatQuad{nearPlaneWidth, nearPlaneHeight, -nearPlaneWidth / 2, nearPlaneHeight / 2};
            verbose("[+!] Camera Type: lookAt Parsed Successfully");
            verbose("[+] Gaze vector calculated: " + std::to_string(newCamera.gaze.x) + " " + std::to_string(newCamera.gaze.y) + " " + std::to_string(newCamera.gaze.z));
            verbose("[+] NearPlane calculated: " + std::to_string(newCamera.nearPlane.x) + " " + std::to_string(newCamera.nearPlane.y) + " " + std::to_string(newCamera.nearPlane.z) + " " + std::to_string(newCamera.nearPlane.w));
            verbose("[+] Camera Parsed Successfully");
            break;
        }
        case 0:
        default: {
            newCamera._id = parseSingleValue<unsigned int>(cameraData["_id"]);
            newCamera.position = parseTriplet<VectorFloatTriplet>(cameraData["Position"]);
            newCamera.gaze = parseTriplet<VectorFloatTriplet>(cameraData["Gaze"]);
            newCamera.up = parseTriplet<VectorFloatTriplet>(cameraData["Up"]);
            newCamera.nearPlane = parseQuad<VectorFloatQuad>(cameraData["NearPlane"]);
            newCamera.nearDistance = parseSingleValue<float>(cameraData["NearDistance"]);
            newCamera.imageResolution = parsePair<VectorIntPair>(cameraData["ImageResolution"]);
            newCamera.imageName = cameraData["ImageName"].get<std::string>();
            break;
        }
    }
    if (dotProduct(newCamera.gaze, newCamera.up) != 0) {
        VectorFloatTriplet w = normalize(-newCamera.gaze);
        VectorFloatTriplet vPrime = normalize(newCamera.up);
        VectorFloatTriplet u = crossProduct(vPrime, w);
        VectorFloatTriplet v = crossProduct(w, u);
        newCamera.up = v;
        verbose("[+] Gaze and Up vectors are not perpendicular. Correcting the up vector.");
    }
    return newCamera;
}

scene::PointLight scene::parsePointLight(const json& pointLightData) {
    scene::PointLight newPointLight;
    newPointLight._id = parseSingleValue<unsigned int>(pointLightData["_id"]);
    newPointLight.position = parseTriplet<VectorFloatTriplet>(pointLightData["Position"]);
    newPointLight.intensity = parseTriplet<VectorFloatTriplet>(pointLightData["Intensity"]);
    return newPointLight;
}

scene::Material scene::parseMaterial(const json& materialData) {
    scene::Material newMaterial;
    newMaterial._id = parseSingleValue<unsigned int>(materialData["_id"]);
    newMaterial.ambientReflectance = parseTriplet<VectorFloatTriplet>(materialData["AmbientReflectance"]);
    newMaterial.diffuseReflectance = parseTriplet<VectorFloatTriplet>(materialData["DiffuseReflectance"]);
    newMaterial.specularReflectance = parseTriplet<VectorFloatTriplet>(materialData["SpecularReflectance"]);
    if (materialData.contains("PhongExponent")) {
        newMaterial.phongExponent = parseSingleValue<float>(materialData["PhongExponent"]);
    } else {
        newMaterial.phongExponent = 0;
    }
    // Parse material type and Fresnel properties
    if (materialData.contains("_type") && !materialData["_type"].is_null()) {
        newMaterial.type = materialData["_type"].get<std::string>();
        newMaterial.isMirror = (newMaterial.type == "mirror");
        if (materialData.contains("MirrorReflectance")) {
            newMaterial.mirrorReflectance = parseTriplet<VectorFloatTriplet>(materialData["MirrorReflectance"]);
        }
        if (materialData.contains("RefractionIndex")) {
            newMaterial.refractionIndex = parseSingleValue<float>(materialData["RefractionIndex"]);
        }
        if (materialData.contains("AbsorptionIndex")) {
            newMaterial.absorptionIndex = parseSingleValue<float>(materialData["AbsorptionIndex"]);
        }
        if (materialData.contains("AbsorptionCoefficient")) {
            newMaterial.absorptionCoefficient = parseTriplet<VectorFloatTriplet>(materialData["AbsorptionCoefficient"]);
        }
    }
    return newMaterial;
}

std::vector<scene::VectorFloatTriplet> scene::parseVertex(const json& vertexData) {
    std::stringstream stream(vertexData.get<std::string>());
    std::vector<scene::VectorFloatTriplet> vertices;
    scene::VectorFloatTriplet vertex;
    while (stream >> vertex.x >> vertex.y >> vertex.z) {
        vertices.push_back(vertex);
    }
    stream.clear();
    return vertices;   
}

std::vector<scene::VectorIntTriplet> scene::Scene::parseFaces(const json& facesData) {
    if (facesData.contains("_data")) {
        auto facesDataArray = facesData["_data"];
        std::stringstream stream(facesDataArray.get<std::string>());
        std::vector<scene::VectorIntTriplet> faces;
        scene::VectorIntTriplet face;
        while (stream >> face.x >> face.y >> face.z) {
            // Convert from 1-based to 0-based indexing
            face.x -= 1;
            face.y -= 1;
            face.z -= 1;
            faces.push_back(face);
        }
        stream.clear();
        return faces;
    } else if (facesData.contains("_plyFile")) {
        std::string plyFile = facesData["_plyFile"].get<std::string>();
        std::string fullPath = this->baseDirectory + plyFile;
        // parsePLYFile will add vertices to this->vertices and return adjusted faces
        std::vector<scene::VectorIntTriplet> faces = parsePLYFile(fullPath, this->vertices);
        return faces;
    } else {
        verbose("[!] Skipping Faces Parsing. Reason: Not found in the scene file.");
        return std::vector<scene::VectorIntTriplet>();
    }
}

template<typename T> 
std::vector<T> scene::Scene::parseObjects(const json& objectsData) {
    std::vector<T> objects;
    if (objectsData.is_array()) {
        for (auto objectData : objectsData) {
            T newObject;
            newObject._id = parseSingleValue<unsigned int>(objectData["_id"]);
            unsigned int materialId = parseSingleValue<unsigned int>(objectData["Material"]);
            newObject.material = getMaterialById(materialId);
            parseSpecificAttributes<T>(newObject, objectData);
            objects.push_back(newObject);
        }
    } else if (!objectsData.is_null()) {
        T newObject;
        newObject._id = parseSingleValue<unsigned int>(objectsData["_id"]);
        unsigned int materialId = parseSingleValue<unsigned int>(objectsData["Material"]);
        newObject.material = getMaterialById(materialId);
        parseSpecificAttributes<T>(newObject, objectsData);
        objects.push_back(newObject);
    }
    return objects;
}

// Specialized parsing for each object type
template<>
void scene::Scene::parseSpecificAttributes<scene::Mesh>(scene::Mesh& object, const json& objectData) {
    char shadingMode = 'f';
    if (objectData.contains("_shadingMode") && !objectData["_shadingMode"].is_null()) {
        object.shadingMode = tolower(objectData["_shadingMode"].get<std::string>()[0]);
    }
    if (objectData.contains("Faces") && !objectData["Faces"].is_null()) {
        object.faces = parseFaces(objectData["Faces"]);
    }
}

template<>
void scene::Scene::parseSpecificAttributes<scene::Triangle>(scene::Triangle& object, const json& objectData) {
    if (objectData.contains("Indices")) {
        object.indices = parseTriplet<VectorIntTriplet>(objectData["Indices"]);
        // Convert from 1-based to 0-based indexing
        object.indices.x -= 1;
        object.indices.y -= 1;
        object.indices.z -= 1;
    }
}

template<>
void scene::Scene::parseSpecificAttributes<scene::Sphere>(scene::Sphere& object, const json& objectData) {
    if (objectData.contains("Center")) {
        object.center = parseSingleValue<unsigned int>(objectData["Center"]);
        // Convert from 1-based to 0-based indexing
        object.center -= 1;
    }
    if (objectData.contains("Radius")) {
        object.radius = parseSingleValue<float>(objectData["Radius"]);
    }
}

template<>
void scene::Scene::parseSpecificAttributes<scene::Plane>(scene::Plane& object, const json& objectData) {
    if (objectData.contains("Point")) {
        object.point = parseSingleValue<unsigned int>(objectData["Point"]);
        // Convert from 1-based to 0-based indexing
        object.point -= 1;
    }
    if (objectData.contains("Normal")) {
        object.normal = parseTriplet<VectorFloatTriplet>(objectData["Normal"]);
    }
}

scene::Material* scene::Scene::getMaterialById(unsigned int id) {
    auto it = materialIdToIndex.find(id);
    if (it != materialIdToIndex.end()) {
        return &materials[it->second];
    }
    return nullptr;
}


void scene::Scene::getSummary() {
    verbose("Scene:");
    verbose("BackgroundColor: " + std::to_string(this->backgroundColor.x) + " " + std::to_string(this->backgroundColor.y) + " " + std::to_string(this->backgroundColor.z));
    verbose("MaxRecursionDepth " + std::to_string(this->maxRecursionDepth));
    verbose("ShadowRayEpsilon: " + std::to_string(this->shadowRayEpsilon));
    verbose("IntersectionTestEpsilon: " + std::to_string(this->intersectionTestEpsilon));
    verbose("Cameras: " + std::to_string(this->cameras.size()));
    for (auto camera : this->cameras) {
        verbose("\t Camera: " + std::to_string(camera._id) + "| Position: " + std::to_string(camera.position.x) + " " + std::to_string(camera.position.y) + " " + std::to_string(camera.position.z) + "| Gaze: " + std::to_string(camera.gaze.x) + " " + std::to_string(camera.gaze.y) + " " + std::to_string(camera.gaze.z) + "| Up: " + std::to_string(camera.up.x) + " " + std::to_string(camera.up.y) + " " + std::to_string(camera.up.z) + "| NearPlane: " + std::to_string(camera.nearPlane.x) + " " + std::to_string(camera.nearPlane.y) + " " + std::to_string(camera.nearPlane.z) + " " + std::to_string(camera.nearPlane.w) + "| NearDistance: " + std::to_string(camera.nearDistance) + "| ImageResolution: " + std::to_string(camera.imageResolution.x) + " " + std::to_string(camera.imageResolution.y) + "| ImageName: " + camera.imageName);
    }
    verbose("Lights: " + std::to_string(this->pointLights.size()) + "| lights: ");
    for (auto light : this->pointLights) {
        verbose("\t Light: " + std::to_string(light._id) + "| Position: " + std::to_string(light.position.x) + " " + std::to_string(light.position.y) + " " + std::to_string(light.position.z) + "| Intensity: " + std::to_string(light.intensity.x) + " " + std::to_string(light.intensity.y) + " " + std::to_string(light.intensity.z));
    }
    verbose("Materials: " + std::to_string(this->materials.size()) + "| materials: ");
    for (auto material : this->materials) {
        verbose("\t Material: " + std::to_string(material._id) + "| AmbientReflectance: " + std::to_string(material.ambientReflectance.x) + " " + std::to_string(material.ambientReflectance.y) + " " + std::to_string(material.ambientReflectance.z) + "| DiffuseReflectance: " + std::to_string(material.diffuseReflectance.x) + " " + std::to_string(material.diffuseReflectance.y) + " " + std::to_string(material.diffuseReflectance.z) + "| SpecularReflectance: " + std::to_string(material.specularReflectance.x) + " " + std::to_string(material.specularReflectance.y) + " " + std::to_string(material.specularReflectance.z) + "| PhongExponent: " + std::to_string(material.phongExponent) + "| isMirror: " + std::to_string(material.isMirror));
    }
    verbose("VertexData: " + std::to_string(this->vertices.size()));
    verbose("Meshes: " + std::to_string(this->meshes.size()));
    for (auto mesh : this->meshes) {
        verbose("\t Mesh: " + std::to_string(mesh._id) + "| Faces: " + std::to_string(mesh.faces.size()));
    }
    verbose("Triangles: " + std::to_string(this->triangles.size()));
    for (auto triangle : this->triangles) {
        verbose("\t Triangle: " + std::to_string(triangle._id) + "| Indices: " + std::to_string(triangle.indices.x) + " " + std::to_string(triangle.indices.y) + " " + std::to_string(triangle.indices.z));
    }
    verbose("Spheres: " + std::to_string(this->spheres.size()));
    for (auto sphere : this->spheres) {
        verbose("\t Sphere: " + std::to_string(sphere._id) + "| Center: " + std::to_string(sphere.center) + " " + std::to_string(sphere.radius));
    }
    verbose("Planes: " + std::to_string(this->planes.size()));
    for (auto plane : this->planes) {
        verbose("\t Plane: " + std::to_string(plane._id) + "| Point Index: " + std::to_string(plane.point) + "| Normal: " + std::to_string(plane.normal.x) + " " + std::to_string(plane.normal.y) + " " + std::to_string(plane.normal.z));
    }
}

void scene::Scene::writePPM(const std::string& filename, unsigned char* image, int width, int height) {
    // @TODO: Remove following lines to write png file, I did this for debugging
    size_t dotPos = filename.find_last_of('.');
    std::string ppmFilename;
    if (dotPos != std::string::npos) {
        ppmFilename = filename.substr(0, dotPos) + ".ppm";
    }
    FILE *outfile;
    if ((outfile = fopen(ppmFilename.c_str(), "w")) == NULL) {
        throw std::runtime_error("Error: The ppm file cannot be opened for writing: " + filename);
    }
    fprintf(outfile, "P3\n%d %d\n255\n", width, height);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < 3; c++) {
                fprintf(outfile, "%u ", image[(y * width + x) * 3 + c]);
            }
        }
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
    fclose(outfile);
}

std::vector<scene::VectorIntTriplet> scene::parsePLYFile(const std::string& plyFile, std::vector<scene::VectorFloatTriplet>& vertexList) {
    if (!std::ifstream(plyFile).good()) {
        throw std::runtime_error("Error: PLY file does not exist: " + plyFile);
    }
    happly::PLYData plyData(plyFile);
    std::vector<std::vector<unsigned long>> faces = plyData.getFaceIndices();
    std::vector<scene::VectorIntTriplet> facesVector = std::vector<scene::VectorIntTriplet>(faces.size());
    for (size_t i = 0; i < faces.size(); i++) {
        auto face = faces[i];
        scene::VectorIntTriplet faceVector;
        faceVector.x = face[0];
        faceVector.y = face[1];
        faceVector.z = face[2];
        facesVector[i] = faceVector;
    }
    std::vector<std::array<double, 3>> vertices = plyData.getVertexPositions();
    std::vector<scene::VectorFloatTriplet> verticesVector = std::vector<scene::VectorFloatTriplet>(vertices.size());
    for (size_t i = 0; i < vertices.size(); i++) {
        auto vertex = vertices[i];
        scene::VectorFloatTriplet vertexVector;
        vertexVector.x = vertex[0];
        vertexVector.y = vertex[1];
        vertexVector.z = vertex[2];
        verticesVector[i] = vertexVector;
    }
    vertexList.insert(vertexList.end(), verticesVector.begin(), verticesVector.end());
    return facesVector;
}


void scene::Scene::buildBVH() {
    verbose("================================================");
    verbose("Building BVH acceleration structures...");
    verbose("================================================");
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    for (auto bvh : meshBVHs) {
        delete bvh;
    }
    meshBVHs.clear();
    
    meshBVHs.resize(meshes.size());
    for (size_t i = 0; i < meshes.size(); i++) {
        if (meshes[i].faces.size() > 0) {
            verbose("[BVH] Building BVH for mesh " + std::to_string(i) + " (" + 
                   std::to_string(meshes[i].faces.size()) + " faces)...");
            
            meshBVHs[i] = new scene::MeshBVH();
            meshBVHs[i]->build(meshes[i], vertices);
        } else {
            meshBVHs[i] = nullptr;
        }
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    verbose("================================================");
    verbose("BVH construction complete in " + std::to_string(duration.count()) + " milliseconds");
    verbose("================================================");
}