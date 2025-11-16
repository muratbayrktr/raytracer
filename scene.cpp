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
        this->shadowRayEpsilon = parseSingleValue<double>(shadowRayEpsilon);
        verbose("[+] ShadowRayEpsilon Parsed: " + std::to_string(this->shadowRayEpsilon));
    } else {
        this->shadowRayEpsilon = 0.001f;
        verbose("[!] Skipping ShadowRayEpsilon Parsing. Reason: Not found in the scene file. Assigning default value: " + std::to_string(shadowRayEpsilon));
    }

    if (scene.contains("IntersectionTestEpsilon") && !scene["IntersectionTestEpsilon"].is_null()) {
        std::string intersectionTestEpsilonStr = scene["IntersectionTestEpsilon"].get<std::string>();
        this->intersectionTestEpsilon = parseSingleValue<double>(intersectionTestEpsilonStr);
        verbose("[+] IntersectionTestEpsilon Parsed: " + std::to_string(this->intersectionTestEpsilon));
    } else {
        this->intersectionTestEpsilon = 1e-10f;
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

    // Parse Transformations
    if (scene.contains("Transformations") && !scene["Transformations"].is_null()) {
        auto transformations = scene["Transformations"];
        
        // Parse Translations
        if (transformations.contains("Translation") && !transformations["Translation"].is_null()) {
            auto translationData = transformations["Translation"];
            if (translationData.is_array()) {
                for (auto trans : translationData) {
                    scene::Translation newTranslation = parseTranslation(trans);
                    this->translations.push_back(newTranslation);
                    this->translationIdToIndex[newTranslation._id] = this->translations.size() - 1;
                    verbose("[+] Translation Parsed: " + std::to_string(newTranslation._id));
                }
            } else {
                scene::Translation newTranslation = parseTranslation(translationData);
                this->translations.push_back(newTranslation);
                this->translationIdToIndex[newTranslation._id] = this->translations.size() - 1;
                verbose("[+] Translation Parsed: " + std::to_string(newTranslation._id));
            }
        }
        
        // Parse Scalings
        if (transformations.contains("Scaling") && !transformations["Scaling"].is_null()) {
            auto scalingData = transformations["Scaling"];
            if (scalingData.is_array()) {
                for (auto scale : scalingData) {
                    scene::Scaling newScaling = parseScaling(scale);
                    this->scalings.push_back(newScaling);
                    this->scalingIdToIndex[newScaling._id] = this->scalings.size() - 1;
                    verbose("[+] Scaling Parsed: " + std::to_string(newScaling._id));
                }
            } else {
                scene::Scaling newScaling = parseScaling(scalingData);
                this->scalings.push_back(newScaling);
                this->scalingIdToIndex[newScaling._id] = this->scalings.size() - 1;
                verbose("[+] Scaling Parsed: " + std::to_string(newScaling._id));
            }
        }
        
        // Parse Rotations
        if (transformations.contains("Rotation") && !transformations["Rotation"].is_null()) {
            auto rotationData = transformations["Rotation"];
            if (rotationData.is_array()) {
                for (auto rot : rotationData) {
                    scene::Rotation newRotation = parseRotation(rot);
                    this->rotations.push_back(newRotation);
                    this->rotationIdToIndex[newRotation._id] = this->rotations.size() - 1;
                    verbose("[+] Rotation Parsed: " + std::to_string(newRotation._id));
                }
            } else {
                scene::Rotation newRotation = parseRotation(rotationData);
                this->rotations.push_back(newRotation);
                this->rotationIdToIndex[newRotation._id] = this->rotations.size() - 1;
                verbose("[+] Rotation Parsed: " + std::to_string(newRotation._id));
            }
        }
        
        // Parse Composites
        if (transformations.contains("Composite") && !transformations["Composite"].is_null()) {
            auto compositeData = transformations["Composite"];
            if (compositeData.is_array()) {
                for (auto comp : compositeData) {
                    scene::Composite newComposite = parseComposite(comp);
                    this->composites.push_back(newComposite);
                    this->compositeIdToIndex[newComposite._id] = this->composites.size() - 1;
                    verbose("[+] Composite Parsed: " + std::to_string(newComposite._id));
                }
            } else {
                scene::Composite newComposite = parseComposite(compositeData);
                this->composites.push_back(newComposite);
                this->compositeIdToIndex[newComposite._id] = this->composites.size() - 1;
                verbose("[+] Composite Parsed: " + std::to_string(newComposite._id));
            }
        }
        
        verbose("[+] Transformations Parsed: T=" + std::to_string(this->translations.size()) + 
                " S=" + std::to_string(this->scalings.size()) + 
                " R=" + std::to_string(this->rotations.size()) + 
                " C=" + std::to_string(this->composites.size()));
    } else {
        verbose("[!] Skipping Transformations Parsing. Reason: Not found in the scene file.");
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
        auto meshInstances = objects["MeshInstance"];

        this->meshes = parseObjects<scene::Mesh>(meshes);
        this->triangles = parseObjects<scene::Triangle>(triangles);
        this->spheres = parseObjects<scene::Sphere>(spheres);
        this->planes = parseObjects<scene::Plane>(planes);
        this->meshInstances = parseObjects<scene::MeshInstance>(meshInstances);

        verbose("[+] Meshes Parsed: " + std::to_string(this->meshes.size()));
        verbose("[+] Triangles Parsed: " + std::to_string(this->triangles.size()));
        verbose("[+] Spheres Parsed: " + std::to_string(this->spheres.size()));
        verbose("[+] Planes Parsed: " + std::to_string(this->planes.size()));
        verbose("[+] MeshInstances Parsed: " + std::to_string(this->meshInstances.size()));
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
            newCamera.nearDistance = parseSingleValue<double>(cameraData["NearDistance"]);
            newCamera.imageResolution = parsePair<VectorIntPair>(cameraData["ImageResolution"]);
            newCamera.imageName = cameraData["ImageName"].get<std::string>();
            // Calculate the gaze vector, nearplane
            // Gaze: GazePoint - Position
            VectorFloatTriplet gazePoint = parseTriplet<VectorFloatTriplet>(cameraData["GazePoint"]);
            VectorFloatTriplet gaze = normalize(gazePoint - newCamera.position);
            newCamera.gaze = gaze;

            // NearPlane:
            // FovY is the vertical field of view in degrees. We compute the near-plane
            // height from it, then use the aspect ratio to get the width.
            double fovY = parseSingleValue<double>(cameraData["FovY"]) * M_PI / 180.0; // convert degrees to radians
            double nearPlaneHeight = 2.0 * tan(fovY / 2.0) * newCamera.nearDistance;
            double aspect = static_cast<double>(newCamera.imageResolution.x) / static_cast<double>(newCamera.imageResolution.y);
            double nearPlaneWidth = nearPlaneHeight * aspect;

            // The renderer expects nearPlane = (l, r, b, t)
            double l = -nearPlaneWidth / 2.0;
            double r =  nearPlaneWidth / 2.0;
            double b = -nearPlaneHeight / 2.0;
            double t =  nearPlaneHeight / 2.0;
            newCamera.nearPlane = VectorFloatQuad{l, r, b, t};
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
            newCamera.nearDistance = parseSingleValue<double>(cameraData["NearDistance"]);
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
    
    // Parse transformations if present
    if (cameraData.contains("Transformations") && !cameraData["Transformations"].is_null()) {
        std::string transformStr = cameraData["Transformations"].get<std::string>();
        newCamera.transformations = parseTransformationString(transformStr);
        verbose("[+] Camera Transformations parsed: " + transformStr + " (" + std::to_string(newCamera.transformations.size()) + " transforms)");
    }
    
    return newCamera;
}

scene::PointLight scene::parsePointLight(const json& pointLightData) {
    scene::PointLight newPointLight;
    newPointLight._id = parseSingleValue<unsigned int>(pointLightData["_id"]);
    newPointLight.position = parseTriplet<VectorFloatTriplet>(pointLightData["Position"]);
    newPointLight.intensity = parseTriplet<VectorFloatTriplet>(pointLightData["Intensity"]);
    
    // Parse transformations if present
    if (pointLightData.contains("Transformations") && !pointLightData["Transformations"].is_null()) {
        std::string transformStr = pointLightData["Transformations"].get<std::string>();
        newPointLight.transformations = parseTransformationString(transformStr);
        verbose("[+] PointLight Transformations parsed: " + transformStr + " (" + std::to_string(newPointLight.transformations.size()) + " transforms)");
    }
    
    return newPointLight;
}

scene::Material scene::parseMaterial(const json& materialData) {
    scene::Material newMaterial;
    newMaterial._id = parseSingleValue<unsigned int>(materialData["_id"]);
    newMaterial.ambientReflectance = parseTriplet<VectorFloatTriplet>(materialData["AmbientReflectance"]);
    newMaterial.diffuseReflectance = parseTriplet<VectorFloatTriplet>(materialData["DiffuseReflectance"]);
    newMaterial.specularReflectance = parseTriplet<VectorFloatTriplet>(materialData["SpecularReflectance"]);
    if (materialData.contains("PhongExponent")) {
        newMaterial.phongExponent = parseSingleValue<double>(materialData["PhongExponent"]);
    } else {
        newMaterial.phongExponent = 0;
    }
    
    // Parse MirrorReflectance if present (independent of _type)
    if (materialData.contains("MirrorReflectance") && !materialData["MirrorReflectance"].is_null()) {
        newMaterial.mirrorReflectance = parseTriplet<VectorFloatTriplet>(materialData["MirrorReflectance"]);
        // Set isMirror if mirror reflectance is non-zero
        if (newMaterial.mirrorReflectance.x > 0.0 || newMaterial.mirrorReflectance.y > 0.0 || newMaterial.mirrorReflectance.z > 0.0)  {
            newMaterial.isMirror = true;
        }
    }
    
    // Parse material type and other Fresnel properties
    if (materialData.contains("_type") && !materialData["_type"].is_null()) {
        newMaterial.type = materialData["_type"].get<std::string>();
        // Override isMirror if type is explicitly "mirror"
        if (newMaterial.type == "mirror") {
            newMaterial.isMirror = true;
        }
        if (materialData.contains("RefractionIndex")) {
            newMaterial.refractionIndex = parseSingleValue<double>(materialData["RefractionIndex"]);
        }
        if (materialData.contains("AbsorptionIndex")) {
            newMaterial.absorptionIndex = parseSingleValue<double>(materialData["AbsorptionIndex"]);
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
            
            // Parse Material if present (optional for MeshInstance)
            if (objectData.contains("Material") && !objectData["Material"].is_null()) {
                unsigned int materialId = parseSingleValue<unsigned int>(objectData["Material"]);
                newObject.material = getMaterialById(materialId);
            } else {
                newObject.material = nullptr;
            }
            
            // Parse transformations if present
            if (objectData.contains("Transformations") && !objectData["Transformations"].is_null()) {
                std::string transformStr = objectData["Transformations"].get<std::string>();
                newObject.transformations = parseTransformationString(transformStr);
                verbose("[+] Object Transformations parsed: " + transformStr + " (" + std::to_string(newObject.transformations.size()) + " transforms)");
            }
            
            parseSpecificAttributes<T>(newObject, objectData);
            objects.push_back(newObject);
        }
    } else if (!objectsData.is_null()) {
        T newObject;
        newObject._id = parseSingleValue<unsigned int>(objectsData["_id"]);
        
        // Parse Material if present (optional for MeshInstance)
        if (objectsData.contains("Material") && !objectsData["Material"].is_null()) {
            unsigned int materialId = parseSingleValue<unsigned int>(objectsData["Material"]);
            newObject.material = getMaterialById(materialId);
        } else {
            const Mesh* baseMesh = findMeshOrInstanceById(newObject._id);
            if (baseMesh) {
                newObject.material = baseMesh->material;
            } else {
                verbose("[!] Skipping Object Material Parsing. Reason: Base mesh not found");
                newObject.material = nullptr;
            }
        }
        
        // Parse transformations if present
        if (objectsData.contains("Transformations") && !objectsData["Transformations"].is_null()) {
            std::string transformStr = objectsData["Transformations"].get<std::string>();
            newObject.transformations = parseTransformationString(transformStr);
            verbose("[+] Object Transformations parsed: " + transformStr + " (" + std::to_string(newObject.transformations.size()) + " transforms)");
        }
        
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
        object.radius = parseSingleValue<double>(objectData["Radius"]);
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

template<>
void scene::Scene::parseSpecificAttributes<scene::MeshInstance>(scene::MeshInstance& object, const json& objectData) {
    if (objectData.contains("_baseMeshId")) {
        object.baseMeshId = parseSingleValue<unsigned int>(objectData["_baseMeshId"]);
    }
    if (objectData.contains("_resetTransform") && !objectData["_resetTransform"].is_null()) {
        std::string resetTransformStr = objectData["_resetTransform"].get<std::string>();
        object.resetTransform = (resetTransformStr == "true" || resetTransformStr == "True" || resetTransformStr == "TRUE" || resetTransformStr == "1");
    } else {
        object.resetTransform = false;
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
    
    // Remember how many vertices we already have so we can offset indices
    // from the PLY file to point into the combined vertex list.
    size_t baseIndex = vertexList.size();

    std::vector<std::vector<unsigned long>> faces = plyData.getFaceIndices();
    std::vector<scene::VectorIntTriplet> facesVector = std::vector<scene::VectorIntTriplet>(faces.size());
    for (size_t i = 0; i < faces.size(); i++) {
        auto face = faces[i];
        scene::VectorIntTriplet faceVector;
        // PLY indices are 0-based; shift them by baseIndex so they refer to
        // the vertices we append below.
        faceVector.x = static_cast<int>(baseIndex + face[0]);
        faceVector.y = static_cast<int>(baseIndex + face[1]);
        faceVector.z = static_cast<int>(baseIndex + face[2]);
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

scene::Translation scene::parseTranslation(const json& translationData) {
    scene::Translation newTranslation;
    newTranslation._id = parseSingleValue<unsigned int>(translationData["_id"]);
    newTranslation.data = parseTriplet<VectorFloatTriplet>(translationData["_data"]);
    return newTranslation;
}

scene::Scaling scene::parseScaling(const json& scalingData) {
    scene::Scaling newScaling;
    newScaling._id = parseSingleValue<unsigned int>(scalingData["_id"]);
    newScaling.data = parseTriplet<VectorFloatTriplet>(scalingData["_data"]);
    return newScaling;
}

scene::Rotation scene::parseRotation(const json& rotationData) {
    scene::Rotation newRotation;
    newRotation._id = parseSingleValue<unsigned int>(rotationData["_id"]);
    std::string rotationStr = rotationData["_data"].get<std::string>();
    std::istringstream stream(rotationStr);
    stream >> newRotation.angle >> newRotation.axis.x >> newRotation.axis.y >> newRotation.axis.z;
    stream.clear();
    return newRotation;
}

scene::Composite scene::parseComposite(const json& compositeData) {
    scene::Composite newComposite;
    newComposite._id = parseSingleValue<unsigned int>(compositeData["_id"]);
    std::string compositeStr = compositeData["_data"].get<std::string>();
    std::istringstream stream(compositeStr);
    for (int i = 0; i < 16; i++) {
        stream >> newComposite.data[i];
    }
    stream.clear();
    return newComposite;
}


std::vector<scene::TransformationRef> scene::parseTransformationString(const std::string& transformStr) {
    std::vector<scene::TransformationRef> transformations;
    if (transformStr.empty()) {
        return transformations;
    }
    
    std::istringstream stream(transformStr);
    std::string token;
    
    while (stream >> token) {
        if (token.empty()) continue;
        
        scene::TransformationRef ref;
        ref.type = tolower(token[0]); // First character is the type: t, s, r, or c
        
        // Parse the id from the rest of the string
        std::string idStr = token.substr(1);
        ref.id = parseSingleValue<unsigned int>(idStr);
        
        transformations.push_back(ref);
        
        verbose("[+] Parsed transformation: type=" + std::string(1, ref.type) + " id=" + std::to_string(ref.id));
    }
    
    return transformations;
}

const Mesh* scene::Scene::findMeshOrInstanceById(unsigned int id) const {
    for (const auto& mesh : meshes) {
        if (mesh._id == id) {
            return &mesh;
        }
    }
    for (const auto& instance : meshInstances) {
        if (instance._id == id) {
            return instance.baseMesh;
        }
    }
    return nullptr;
}

int scene::Scene::findBaseMeshIndex(unsigned int meshOrInstanceId) const {
    for (size_t i = 0; i < meshes.size(); i++) {
        if (meshes[i]._id == meshOrInstanceId) {
            return i;
        }
    }
    for (const auto& instance : meshInstances) {
        if (instance._id == meshOrInstanceId && instance.baseMesh) {
            return findBaseMeshIndex(instance.baseMeshId);
        }
    }
    return -1;
}

template<typename T>
void processObjectTransformations(std::vector<T>& objects, Scene& scene) {
    for (auto& obj : objects) {
        if (!obj.transformations.empty()) {
            obj.transformMatrix = new Matrix4x4(buildObjectTransformMatrix(scene, obj.transformations));
            obj.inverseTransformMatrix = new Matrix4x4(invertMatrix(*obj.transformMatrix));
            Matrix4x4 invTrans = transposeMatrix(*obj.inverseTransformMatrix);
            obj.normalMatrix = new Matrix4x4(invTrans);
            obj.hasTransformation = true;
            obj.hasNegativeScale = hasNegativeScale(*obj.transformMatrix);
            if (obj.hasNegativeScale) {
                verbose("[+] Object " + std::to_string(obj._id) + " has negative scale (reflection)");
            }
            verbose("[+] Computed transformation for object " + std::to_string(obj._id));
        } else {
            obj.transformMatrix = nullptr;
            obj.inverseTransformMatrix = nullptr;
            obj.normalMatrix = nullptr;
            obj.hasTransformation = false;
            obj.hasNegativeScale = false;
        }
    }
}

void scene::Scene::precomputeTransformations() {
    verbose("================================================");
    verbose("Precomputing object transformations...");
    verbose("================================================");
    
    // BVH must be built first so we can use it to get bounding boxes
    if (meshBVHs.empty()) {
        verbose("[WARNING] BVH not built yet, world-space bounds optimization will be limited");
    }
    
    processObjectTransformations(meshes, *this);
    processObjectTransformations(triangles, *this);
    processObjectTransformations(spheres, *this);
    processObjectTransformations(planes, *this);
    
    // Compute world-space bounding boxes for transformed meshes
    for (size_t i = 0; i < meshes.size(); i++) {
        auto& mesh = meshes[i];
        if (mesh.hasTransformation && mesh.transformMatrix) {
            AABB localBounds;
            if (i < meshBVHs.size() && meshBVHs[i] != nullptr) {
                localBounds = meshBVHs[i]->getRootBounds();
            } else {
                localBounds = computeMeshAABB(mesh, vertices);
            }
            mesh.worldSpaceBounds = new AABB(localBounds.transform(*mesh.transformMatrix));
            verbose("[+] Computed world-space bounds for mesh " + std::to_string(mesh._id));
        }
    }
    
    for (auto& instance : meshInstances) {
        instance.baseMesh = findMeshOrInstanceById(instance.baseMeshId);
        instance.baseMeshIndex = findBaseMeshIndex(instance.baseMeshId);
        
        if (instance.baseMesh) {
            // Inherit material from base mesh if not specified
            if (instance.material == nullptr) {
                instance.material = instance.baseMesh->material;
                verbose("[+] Mesh instance " + std::to_string(instance._id) + " inherited material from base mesh");
            }
            
            Matrix4x4 finalMatrix;
            bool hasActualTransform = false;
            
            if (instance.resetTransform) {
                if (!instance.transformations.empty()) {
                    finalMatrix = buildObjectTransformMatrix(*this, instance.transformations);
                    hasActualTransform = true;
                }
            } else {
                // Check if baseMeshId refers to another instance (chained instance)
                MeshInstance* parentInstance = nullptr;
                for (auto& otherInst : meshInstances) {
                    if (otherInst._id == instance.baseMeshId) {
                        parentInstance = &otherInst;
                        break;
                    }
                }
                
                bool baseHasTransform = false;
                Matrix4x4 baseMatrix = identityMatrix();
                
                if (parentInstance) {
                    // This instance references another instance
                    baseHasTransform = parentInstance->hasTransformation;
                    if (baseHasTransform && parentInstance->transformMatrix) {
                        baseMatrix = *parentInstance->transformMatrix;
                    }
                } else {
                    // This instance references a mesh directly
                    baseHasTransform = instance.baseMesh->hasTransformation;
                    if (baseHasTransform && instance.baseMesh->transformMatrix) {
                        baseMatrix = *instance.baseMesh->transformMatrix;
                    }
                }
                
                bool instanceHasTransform = !instance.transformations.empty();
                
                if (baseHasTransform || instanceHasTransform) {
                    Matrix4x4 instanceMatrix = instanceHasTransform ? buildObjectTransformMatrix(*this, instance.transformations) : identityMatrix();
                    finalMatrix = multiplyMatrices(instanceMatrix, baseMatrix);
                    hasActualTransform = true;
                }
            }
            
            if (hasActualTransform) {
                instance.transformMatrix = new Matrix4x4(finalMatrix);
                instance.inverseTransformMatrix = new Matrix4x4(invertMatrix(*instance.transformMatrix));
                Matrix4x4 invTrans = transposeMatrix(*instance.inverseTransformMatrix);
                instance.normalMatrix = new Matrix4x4(invTrans);
                instance.hasTransformation = true;
                instance.hasNegativeScale = hasNegativeScale(*instance.transformMatrix);
                if (instance.hasNegativeScale) {
                    verbose("[+] Mesh instance " + std::to_string(instance._id) + " has negative scale (reflection)");
                }
                
                // Compute world-space bounding box for the instance
                int baseMeshIdx = instance.baseMeshIndex;
                if (baseMeshIdx >= 0) {
                    AABB localBounds;
                    if (baseMeshIdx < (int)meshBVHs.size() && meshBVHs[baseMeshIdx] != nullptr) {
                        localBounds = meshBVHs[baseMeshIdx]->getRootBounds();
                    } else {
                        localBounds = computeMeshAABB(*instance.baseMesh, vertices);
                    }
                    instance.worldSpaceBounds = new AABB(localBounds.transform(*instance.transformMatrix));
                    verbose("[+] Computed world-space bounds for mesh instance " + std::to_string(instance._id));
                }
                
                verbose("[+] Computed transformation for mesh instance " + std::to_string(instance._id));
            } else {
                instance.transformMatrix = nullptr;
                instance.inverseTransformMatrix = nullptr;
                instance.normalMatrix = nullptr;
                instance.hasTransformation = false;
                instance.hasNegativeScale = false;
            }
        }
    }
    
    for (auto& camera : cameras) {
        if (!camera.transformations.empty()) {
            Matrix4x4 cameraTransform = buildObjectTransformMatrix(*this, camera.transformations);
            camera.position = transformPoint(cameraTransform, camera.position);
            camera.gaze = normalize(transformDirection(cameraTransform, camera.gaze));
            camera.up = normalize(transformDirection(cameraTransform, camera.up));
            verbose("[+] Applied transformation to camera " + std::to_string(camera._id));
        }
    }
    
    for (auto& light : pointLights) {
        if (!light.transformations.empty()) {
            Matrix4x4 lightTransform = buildObjectTransformMatrix(*this, light.transformations);
            light.position = transformPoint(lightTransform, light.position);
            verbose("[+] Applied transformation to light " + std::to_string(light._id));
        }
    }
    
    verbose("================================================");
}

void scene::Scene::buildBVH() {
    verbose("================================================");
    verbose("Building BVH acceleration structures...");
    verbose("================================================");
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    std::vector<MeshBVH*> newMeshBVHs(meshes.size());
    for (size_t i = 0; i < meshes.size(); i++) {
        if (meshes[i].faces.size() > 0) {
            verbose("[BVH] Building BVH for mesh " + std::to_string(i) + " (" + 
                   std::to_string(meshes[i].faces.size()) + " faces)...");
            
            newMeshBVHs[i] = new scene::MeshBVH();
            newMeshBVHs[i]->build(meshes[i], vertices);
        } else {
            newMeshBVHs[i] = nullptr;
        }
    }
    this->meshBVHs = std::move(newMeshBVHs);
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    verbose("================================================");
    verbose("BVH construction complete in " + std::to_string(duration.count()) + " milliseconds");
    verbose("================================================");
}