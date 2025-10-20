#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "scene.h"
#include "json.hpp"

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
            std::cout << "parsing camera array" << std::endl;
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
    switch (cameraType) {
        case 1:
            // TODO: I'll implement this later i.e. other_dragon.json
            verbose("[-] Skipping camera parsing type is: " + cameraData["_type"].get<std::string>());
            return scene::Camera();
        case 0:
        default:
            scene::Camera newCamera;
            newCamera._id = parseSingleValue<unsigned int>(cameraData["_id"]);
            newCamera.position = parseTriplet<VectorFloatTriplet>(cameraData["Position"]);
            newCamera.gaze = parseTriplet<VectorFloatTriplet>(cameraData["Gaze"]);
            newCamera.up = parseTriplet<VectorFloatTriplet>(cameraData["Up"]);
            newCamera.nearPlane = parseQuad<VectorFloatQuad>(cameraData["NearPlane"]);
            newCamera.nearDistance = parseSingleValue<float>(cameraData["NearDistance"]);
            newCamera.imageResolution = parsePair<VectorIntPair>(cameraData["ImageResolution"]);
            newCamera.imageName = cameraData["ImageName"].get<std::string>();
            return newCamera;
    }
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
    newMaterial.phongExponent = parseSingleValue<float>(materialData["PhongExponent"]);
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

std::vector<int> scene::parseFaces(const json& facesData) {
    auto facesDataArray = facesData["_data"]; // "122 163 1640 623 ..."
    std::stringstream stream(facesDataArray.get<std::string>());
    std::vector<int> faces;
    int face;
    while (stream >> face) {
        faces.push_back(face);
    }
    stream.clear();
    return faces;
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
    if (objectData.contains("Faces")) {
        object.faces = parseFaces(objectData["Faces"]);
    }
}

template<>
void scene::Scene::parseSpecificAttributes<scene::Triangle>(scene::Triangle& object, const json& objectData) {
    if (objectData.contains("Indices")) {
        object.indices = parseTriplet<VectorFloatTriplet>(objectData["Indices"]);
    }
}

template<>
void scene::Scene::parseSpecificAttributes<scene::Sphere>(scene::Sphere& object, const json& objectData) {
    if (objectData.contains("Center")) {
        object.center = parseSingleValue<unsigned int>(objectData["Center"]);
    }
    if (objectData.contains("Radius")) {
        object.radius = parseSingleValue<float>(objectData["Radius"]);
    }
}

template<>
void scene::Scene::parseSpecificAttributes<scene::Plane>(scene::Plane& object, const json& objectData) {
    if (objectData.contains("Point")) {
        object.point = parseSingleValue<unsigned int>(objectData["Point"]);
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
    std::cout << "Scene:" << std::endl;
    std::cout << "BackgroundColor: " << this->backgroundColor.x << " " << this->backgroundColor.y << " " << this->backgroundColor.z << std::endl;
    std::cout << "ShadowRayEpsilon: " << this->shadowRayEpsilon << std::endl;
    std::cout << "IntersectionTestEpsilon: " << this->intersectionTestEpsilon << std::endl;
    std::cout << "Cameras: " << this->cameras.size() << std::endl;
    for (auto camera : this->cameras) {
        std::cout << "\t Camera: " << camera._id << "| Position: " << camera.position.x << " " << camera.position.y << " " << camera.position.z << "| Gaze: " << camera.gaze.x << " " << camera.gaze.y << " " << camera.gaze.z << "| Up: " << camera.up.x << " " << camera.up.y << " " << camera.up.z << "| NearPlane: " << camera.nearPlane.x << " " << camera.nearPlane.y << " " << camera.nearPlane.z << " " << camera.nearPlane.w << "| NearDistance: " << camera.nearDistance << "| ImageResolution: " << camera.imageResolution.x << " " << camera.imageResolution.y << "| ImageName: " << camera.imageName << std::endl;
    }
    std::cout << "Lights: " << this->pointLights.size() << "| lights: " << std::endl;
    for (auto light : this->pointLights) {
        std::cout << "\t Light: " << light._id << "| Position: " << light.position.x << " " << light.position.y << " " << light.position.z << "| Intensity: " << light.intensity.x << " " << light.intensity.y << " " << light.intensity.z << std::endl;
    }
    std::cout << "Materials: " << this->materials.size() << "| materials: " << std::endl;
    for (auto material : this->materials) {
        std::cout << "\t Material: " << material._id << "| AmbientReflectance: " << material.ambientReflectance.x << " " << material.ambientReflectance.y << " " << material.ambientReflectance.z << "| DiffuseReflectance: " << material.diffuseReflectance.x << " " << material.diffuseReflectance.y << " " << material.diffuseReflectance.z << "| SpecularReflectance: " << material.specularReflectance.x << " " << material.specularReflectance.y << " " << material.specularReflectance.z << "| PhongExponent: " << material.phongExponent << std::endl;
    }
    std::cout << "VertexData: " << this->vertices.size() << std::endl;
    std::cout << "Meshes: " << this->meshes.size() << std::endl;
    for (auto mesh : this->meshes) {
        std::cout << "\t Mesh: " << mesh._id << "| Faces: " << mesh.faces.size() << std::endl;
    }
    std::cout << "Triangles: " << this->triangles.size() << std::endl;
    for (auto triangle : this->triangles) {
        std::cout << "\t Triangle: " << triangle._id << "| Indices: " << triangle.indices.x << " " << triangle.indices.y << " " << triangle.indices.z << std::endl;
    }
    std::cout << "Spheres: " << this->spheres.size() << std::endl;
    for (auto sphere : this->spheres) {
        std::cout << "\t Sphere: " << sphere._id << "| Center: " << sphere.center << " " << sphere.radius << std::endl;
    }
    std::cout << "Planes: " << this->planes.size() << std::endl;
    for (auto plane : this->planes) {
        std::cout << "\t Plane: " << plane._id << "| Point: " << plane.point << " " << plane.normal.x << " " << plane.normal.y << " " << plane.normal.z << std::endl;
    }
}

