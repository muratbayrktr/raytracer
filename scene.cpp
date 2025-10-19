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

template <typename T> T parseSingleValue(const std::string& value) {
    T parsedValue;
    std::istringstream stream(value);
    stream >> parsedValue;
    stream.clear();
    return parsedValue;
}

template <typename T> T parsePair(const std::string& value) {
    T parsedValue;
    std::istringstream stream(value);
    stream >> parsedValue.x >> parsedValue.y;
    stream.clear();
    return parsedValue;
}

template <typename T> T parseTriplet(const std::string& value) {
    T triplet;
    std::istringstream stream(value);
    stream >> triplet.x >> triplet.y >> triplet.z;
    stream.clear();
    return triplet;
}

template <typename T> T parseQuad(const std::string& value) {
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
                    verbose("[+] Material Parsed: " + std::to_string(newMaterial._id));
                } else {
                    verbose("[!] Skipping Material Parsing. Reason: Material ID is 0");
                }
            }
        } else {
            scene::Material newMaterial = parseMaterial(materialDataArray);
            if (newMaterial._id != 0) {
                this->materials.push_back(newMaterial);
                verbose("[+] Material Parsed: " + std::to_string(newMaterial._id));
            }
        }
    } else {
        verbose("[!] Skipping Materials Parsing. Reason: Not found in the scene file. Assigning default value: 0");
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