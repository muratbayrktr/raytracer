#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <limits>

#include "scene.h"
#include "json.hpp"
#include "utils.h"
#include "overloads.h"
#include "happly.h"
#include "bvh.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define TINYEXR_USE_STB_ZLIB 1
#define TINYEXR_USE_MINIZ 0
#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"

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
        this->maxRecursionDepth = 1;
        verbose("[!] Skipping MaxRecursionDepth parsing. Reason: Not found in the scene file " + std::to_string(this->maxRecursionDepth)); 
    }

    if (scene.contains("ShadowRayEpsilon") && !scene["ShadowRayEpsilon"].is_null()) {
        std::string shadowRayEpsilon = scene["ShadowRayEpsilon"].get<std::string>();
        this->shadowRayEpsilon = parseSingleValue<double>(shadowRayEpsilon);
        verbose("[+] ShadowRayEpsilon Parsed: " + std::to_string(this->shadowRayEpsilon));
    } else {
        this->shadowRayEpsilon = 0.01;
        verbose("[!] Skipping ShadowRayEpsilon Parsing. Reason: Not found in the scene file. Assigning default value: " + std::to_string(shadowRayEpsilon));
    }

    if (scene.contains("IntersectionTestEpsilon") && !scene["IntersectionTestEpsilon"].is_null()) {
        std::string intersectionTestEpsilonStr = scene["IntersectionTestEpsilon"].get<std::string>();
        this->intersectionTestEpsilon = parseSingleValue<double>(intersectionTestEpsilonStr);
        verbose("[+] IntersectionTestEpsilon Parsed: " + std::to_string(this->intersectionTestEpsilon));
    } else {
        this->intersectionTestEpsilon = 0.0;
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

        if (lights.contains("AreaLight") && !lights["AreaLight"].is_null()) {
            auto areaLightArray = lights["AreaLight"];
            if (areaLightArray.is_array()) {
                for (auto areaLightData : areaLightArray) {
                    scene::AreaLight newAreaLight = parseAreaLight(areaLightData);
                    if (newAreaLight._id != 0) {
                        this->areaLights.push_back(newAreaLight);
                        verbose("[+] AreaLight Parsed: " + std::to_string(newAreaLight._id));
                    }
                }
            } else {
                scene::AreaLight newAreaLight = parseAreaLight(areaLightArray);
                if (newAreaLight._id != 0) {
                    this->areaLights.push_back(newAreaLight);
                    verbose("[+] AreaLight Parsed: " + std::to_string(newAreaLight._id));
                }
            }
        } else {
            verbose("[!] Skipping AreaLight Parsing. Reason: Not found in the scene file.");
        }

        if (lights.contains("DirectionalLight") && !lights["DirectionalLight"].is_null()) {
            auto directionalLightArray = lights["DirectionalLight"];
            if (directionalLightArray.is_array()) {
                for (auto directionalLightData : directionalLightArray) {
                    scene::DirectionalLight newDirectionalLight = parseDirectionalLight(directionalLightData);
                    if (newDirectionalLight._id != 0) {
                        this->directionalLights.push_back(newDirectionalLight);
                        verbose("[+] DirectionalLight Parsed: " + std::to_string(newDirectionalLight._id));
                    }
                }
            } else {
                scene::DirectionalLight newDirectionalLight = parseDirectionalLight(directionalLightArray);
                if (newDirectionalLight._id != 0) {
                    this->directionalLights.push_back(newDirectionalLight);
                    verbose("[+] DirectionalLight Parsed: " + std::to_string(newDirectionalLight._id));
                }
            }
        } else {
            verbose("[!] Skipping DirectionalLight Parsing. Reason: Not found in the scene file.");
        }

        if (lights.contains("SpotLight") && !lights["SpotLight"].is_null()) {
            auto spotLightArray = lights["SpotLight"];
            if (spotLightArray.is_array()) {
                for (auto spotLightData : spotLightArray) {
                    scene::SpotLight newSpotLight = parseSpotLight(spotLightData);
                    if (newSpotLight._id != 0) {
                        this->spotLights.push_back(newSpotLight);
                        verbose("[+] SpotLight Parsed: " + std::to_string(newSpotLight._id));
                    }
                }
            } else {
                scene::SpotLight newSpotLight = parseSpotLight(spotLightArray);
                if (newSpotLight._id != 0) {
                    this->spotLights.push_back(newSpotLight);
                    verbose("[+] SpotLight Parsed: " + std::to_string(newSpotLight._id));
                }
            }
        } else {
            verbose("[!] Skipping SpotLight Parsing. Reason: Not found in the scene file.");
        }

        if (lights.contains("SphericalDirectionalLight") && !lights["SphericalDirectionalLight"].is_null()) {
            auto sphericalDirectionalLightArray = lights["SphericalDirectionalLight"];
            if (sphericalDirectionalLightArray.is_array()) {
                for (auto sphericalDirectionalLightData : sphericalDirectionalLightArray) {
                    scene::SphericalDirectionalLight newSphericalDirectionalLight = parseSphericalDirectionalLight(sphericalDirectionalLightData);
                    if (newSphericalDirectionalLight._id != 0) {
                        this->sphericalDirectionalLights.push_back(newSphericalDirectionalLight);
                        verbose("[+] SphericalDirectionalLight Parsed: " + std::to_string(newSphericalDirectionalLight._id));
                    }
                }
            } else {
                scene::SphericalDirectionalLight newSphericalDirectionalLight = parseSphericalDirectionalLight(sphericalDirectionalLightArray);
                if (newSphericalDirectionalLight._id != 0) {
                    this->sphericalDirectionalLights.push_back(newSphericalDirectionalLight);
                    verbose("[+] SphericalDirectionalLight Parsed: " + std::to_string(newSphericalDirectionalLight._id));
                }
            }
        } else {
            verbose("[!] Skipping SphericalDirectionalLight Parsing. Reason: Not found in the scene file.");
        }
    } else {
        verbose("[!] Skipping Lights Parsing. Reason: Not found in the scene file. Assigning default value: 0");
    }

    // Parse BRDFs
    if (scene.contains("BRDFs") && !scene["BRDFs"].is_null()) {
        auto brdfs = scene["BRDFs"];
        
        // Parse OriginalBlinnPhong
        if (brdfs.contains("OriginalBlinnPhong") && !brdfs["OriginalBlinnPhong"].is_null()) {
            auto brdfArray = brdfs["OriginalBlinnPhong"];
            if (brdfArray.is_array()) {
                for (auto brdfData : brdfArray) {
                    scene::BRDF newBRDF = parseBRDF(brdfData, scene::BRDFType::OriginalBlinnPhong);
                    if (newBRDF._id != 0) {
                        this->brdfs.push_back(newBRDF);
                        this->brdfIdToIndex[newBRDF._id] = this->brdfs.size() - 1;
                        verbose("[+] BRDF OriginalBlinnPhong Parsed: " + std::to_string(newBRDF._id));
                    }
                }
            } else {
                scene::BRDF newBRDF = parseBRDF(brdfArray, scene::BRDFType::OriginalBlinnPhong);
                if (newBRDF._id != 0) {
                    this->brdfs.push_back(newBRDF);
                    this->brdfIdToIndex[newBRDF._id] = this->brdfs.size() - 1;
                    verbose("[+] BRDF OriginalBlinnPhong Parsed: " + std::to_string(newBRDF._id));
                }
            }
        }
        
        // Parse OriginalPhong
        if (brdfs.contains("OriginalPhong") && !brdfs["OriginalPhong"].is_null()) {
            auto brdfArray = brdfs["OriginalPhong"];
            if (brdfArray.is_array()) {
                for (auto brdfData : brdfArray) {
                    scene::BRDF newBRDF = parseBRDF(brdfData, scene::BRDFType::OriginalPhong);
                    if (newBRDF._id != 0) {
                        this->brdfs.push_back(newBRDF);
                        this->brdfIdToIndex[newBRDF._id] = this->brdfs.size() - 1;
                        verbose("[+] BRDF OriginalPhong Parsed: " + std::to_string(newBRDF._id));
                    }
                }
            } else {
                scene::BRDF newBRDF = parseBRDF(brdfArray, scene::BRDFType::OriginalPhong);
                if (newBRDF._id != 0) {
                    this->brdfs.push_back(newBRDF);
                    this->brdfIdToIndex[newBRDF._id] = this->brdfs.size() - 1;
                    verbose("[+] BRDF OriginalPhong Parsed: " + std::to_string(newBRDF._id));
                }
            }
        }
        
        // Parse ModifiedBlinnPhong
        if (brdfs.contains("ModifiedBlinnPhong") && !brdfs["ModifiedBlinnPhong"].is_null()) {
            auto brdfArray = brdfs["ModifiedBlinnPhong"];
            if (brdfArray.is_array()) {
                for (auto brdfData : brdfArray) {
                    scene::BRDF newBRDF = parseBRDF(brdfData, scene::BRDFType::ModifiedBlinnPhong);
                    if (newBRDF._id != 0) {
                        this->brdfs.push_back(newBRDF);
                        this->brdfIdToIndex[newBRDF._id] = this->brdfs.size() - 1;
                        verbose("[+] BRDF ModifiedBlinnPhong Parsed: " + std::to_string(newBRDF._id));
                    }
                }
            } else {
                scene::BRDF newBRDF = parseBRDF(brdfArray, scene::BRDFType::ModifiedBlinnPhong);
                if (newBRDF._id != 0) {
                    this->brdfs.push_back(newBRDF);
                    this->brdfIdToIndex[newBRDF._id] = this->brdfs.size() - 1;
                    verbose("[+] BRDF ModifiedBlinnPhong Parsed: " + std::to_string(newBRDF._id));
                }
            }
        }
        
        // Parse ModifiedPhong
        if (brdfs.contains("ModifiedPhong") && !brdfs["ModifiedPhong"].is_null()) {
            auto brdfArray = brdfs["ModifiedPhong"];
            if (brdfArray.is_array()) {
                for (auto brdfData : brdfArray) {
                    scene::BRDF newBRDF = parseBRDF(brdfData, scene::BRDFType::ModifiedPhong);
                    if (newBRDF._id != 0) {
                        this->brdfs.push_back(newBRDF);
                        this->brdfIdToIndex[newBRDF._id] = this->brdfs.size() - 1;
                        verbose("[+] BRDF ModifiedPhong Parsed: " + std::to_string(newBRDF._id));
                    }
                }
            } else {
                scene::BRDF newBRDF = parseBRDF(brdfArray, scene::BRDFType::ModifiedPhong);
                if (newBRDF._id != 0) {
                    this->brdfs.push_back(newBRDF);
                    this->brdfIdToIndex[newBRDF._id] = this->brdfs.size() - 1;
                    verbose("[+] BRDF ModifiedPhong Parsed: " + std::to_string(newBRDF._id));
                }
            }
        }
        
        // Parse TorranceSparrow
        if (brdfs.contains("TorranceSparrow") && !brdfs["TorranceSparrow"].is_null()) {
            auto brdfArray = brdfs["TorranceSparrow"];
            if (brdfArray.is_array()) {
                for (auto brdfData : brdfArray) {
                    scene::BRDF newBRDF = parseBRDF(brdfData, scene::BRDFType::TorranceSparrow);
                    if (newBRDF._id != 0) {
                        this->brdfs.push_back(newBRDF);
                        this->brdfIdToIndex[newBRDF._id] = this->brdfs.size() - 1;
                        verbose("[+] BRDF TorranceSparrow Parsed: " + std::to_string(newBRDF._id));
                    }
                }
            } else {
                scene::BRDF newBRDF = parseBRDF(brdfArray, scene::BRDFType::TorranceSparrow);
                if (newBRDF._id != 0) {
                    this->brdfs.push_back(newBRDF);
                    this->brdfIdToIndex[newBRDF._id] = this->brdfs.size() - 1;
                    verbose("[+] BRDF TorranceSparrow Parsed: " + std::to_string(newBRDF._id));
                }
            }
        }
    } else {
        verbose("[!] Skipping BRDFs Parsing. Reason: Not found in the scene file.");
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
        if (vertexData.contains("_data") && !vertexData["_data"].is_null()) {
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
        } else if (vertexData.contains("_binaryFile") && !vertexData["_binaryFile"].is_null()) {
            // Parse binary vertex data: 4-byte count + count*12 bytes (xyz floats)
            std::string binaryFile = vertexData["_binaryFile"].get<std::string>();
            std::string fullPath = this->baseDirectory + binaryFile;
            std::ifstream file(fullPath, std::ios::binary);
            if (file.is_open()) {
                uint32_t count;
                file.read(reinterpret_cast<char*>(&count), sizeof(count));
                this->vertices.reserve(this->vertices.size() + count);
                for (uint32_t i = 0; i < count; i++) {
                    float x, y, z;
                    file.read(reinterpret_cast<char*>(&x), sizeof(x));
                    file.read(reinterpret_cast<char*>(&y), sizeof(y));
                    file.read(reinterpret_cast<char*>(&z), sizeof(z));
                    this->vertices.push_back(VectorFloatTriplet{(double)x, (double)y, (double)z});
                }
                file.close();
                verbose("[+] Binary Vertex Parsed: " + std::to_string(count) + " vertices from " + binaryFile);
            } else {
                verbose("[!] Failed to open binary vertex file: " + fullPath);
            }
        } else {
            verbose("[!] Skipping VertexData Parsing. Reason: VertexData missing _data or _binaryFile");
        }
    }

    if (scene.contains("TexCoordData") && !scene["TexCoordData"].is_null()) {
        auto texCoordData = scene["TexCoordData"];
        if (texCoordData.contains("_data") && !texCoordData["_data"].is_null()) {
            try {
                this->texCoords = parseTexCoordData(texCoordData);
                verbose("[+] TexCoordData Parsed: " + std::to_string(this->texCoords.size()) + " UV pairs");
            } catch (const std::exception& e) {
                verbose("[!] Skipping TexCoordData Parsing. Reason: " + std::string(e.what()));
            }
        } else if (texCoordData.contains("_binaryFile") && !texCoordData["_binaryFile"].is_null()) {
            // Parse binary texcoord data: 4-byte count + count*8 bytes (uv floats)
            std::string binaryFile = texCoordData["_binaryFile"].get<std::string>();
            std::string fullPath = this->baseDirectory + binaryFile;
            std::ifstream file(fullPath, std::ios::binary);
            if (file.is_open()) {
                uint32_t count;
                file.read(reinterpret_cast<char*>(&count), sizeof(count));
                this->texCoords.reserve(this->texCoords.size() + count);
                for (uint32_t i = 0; i < count; i++) {
                    float u, v;
                    file.read(reinterpret_cast<char*>(&u), sizeof(u));
                    file.read(reinterpret_cast<char*>(&v), sizeof(v));
                    this->texCoords.push_back(VectorFloatPair{(double)u, (double)v});
                }
                file.close();
                verbose("[+] Binary TexCoord Parsed: " + std::to_string(count) + " UV pairs from " + binaryFile);
            } else {
                verbose("[!] Failed to open binary texcoord file: " + fullPath);
            }
        } else {
            verbose("[!] Skipping TexCoordData Parsing. Reason: TexCoordData missing _data or _binaryFile");
        }
    }

    if (scene.contains("Textures") && !scene["Textures"].is_null()) {
        auto textures = scene["Textures"];
        
        if (textures.contains("Images") && !textures["Images"].is_null()) {
            auto imagesData = textures["Images"];
            if (imagesData.contains("Image") && !imagesData["Image"].is_null()) {
                auto imageArray = imagesData["Image"];
                if (imageArray.is_array()) {
                    for (auto imageData : imageArray) {
                        try {
                            Image newImage = parseImage(imageData);
                            if (newImage._id != 0) {
                                this->images.push_back(newImage);
                                this->imageIdToIndex[newImage._id] = this->images.size() - 1;
                                verbose("[+] Image Parsed: " + std::to_string(newImage._id));
                            }
                        } catch (const std::exception& e) {
                            verbose("[!] Skipping Image Parsing. Reason: " + std::string(e.what()));
                        }
                    }
                } else {
                    try {
                        Image newImage = parseImage(imageArray);
                        if (newImage._id != 0) {
                            this->images.push_back(newImage);
                            this->imageIdToIndex[newImage._id] = this->images.size() - 1;
                            verbose("[+] Image Parsed: " + std::to_string(newImage._id));
                        }
                    } catch (const std::exception& e) {
                        verbose("[!] Skipping Image Parsing. Reason: " + std::string(e.what()));
                    }
                }
            }
        }
        
        // Load image data after parsing all image definitions
        for (auto& image : this->images) {
            std::string fullPath = this->baseDirectory + image.filename;
            int width, height, channels;
            
            // Check file extension to determine if it's HDR
            std::string lowerFilename = image.filename;
            std::transform(lowerFilename.begin(), lowerFilename.end(), lowerFilename.begin(), ::tolower);
            bool isEXR = lowerFilename.length() >= 4 && lowerFilename.substr(lowerFilename.length() - 4) == ".exr";
            bool isHDR = lowerFilename.length() >= 4 && lowerFilename.substr(lowerFilename.length() - 4) == ".hdr";
            
            if (isEXR) {
                // Load EXR using tinyexr
                float* rgba = nullptr;
                const char* err = nullptr;
                int ret = LoadEXR(&rgba, &width, &height, fullPath.c_str(), &err);
                if (ret == TINYEXR_SUCCESS) {
                    image.hdrData = rgba;
                    image.width = width;
                    image.height = height;
                    image.channels = 4;  // EXR loads as RGBA
                    image.isHDR = true;
                    verbose("[+] EXR image loaded: " + image.filename + " (" + std::to_string(width) + "x" + std::to_string(height) + ")");
                } else {
                    if (err) {
                        verbose("[!] Failed to load EXR image: " + fullPath + " - " + std::string(err));
                        FreeEXRErrorMessage(err);
                    } else {
                        verbose("[!] Failed to load EXR image: " + fullPath);
                    }
                    throw std::runtime_error("Failed to load EXR image: " + fullPath);
                }
            } else if (isHDR) {
                // Load HDR using stbi_loadf
                float* data = stbi_loadf(fullPath.c_str(), &width, &height, &channels, 0);
                if (data) {
                    image.hdrData = data;
                    image.width = width;
                    image.height = height;
                    image.channels = channels;
                    image.isHDR = true;
                    verbose("[+] HDR image loaded: " + image.filename + " (" + std::to_string(width) + "x" + std::to_string(height) + ", " + std::to_string(channels) + " channels)");
                } else {
                    verbose("[!] Failed to load HDR image: " + fullPath);
                    throw std::runtime_error("Failed to load HDR image: " + fullPath);
                }
            } else {
                // Load LDR image using stbi_load
                unsigned char* data = stbi_load(fullPath.c_str(), &width, &height, &channels, 0);
                if (data) {
                    image.data = data;
                    image.width = width;
                    image.height = height;
                    image.channels = channels;
                    image.isHDR = false;
                    verbose("[+] Image loaded: " + image.filename + " (" + std::to_string(width) + "x" + std::to_string(height) + ", " + std::to_string(channels) + " channels)");
                } else {
                    verbose("[!] Failed to load image: " + fullPath);
                    throw std::runtime_error("Failed to load image: " + fullPath);
                }
            }
        }
        
        if (textures.contains("TextureMap") && !textures["TextureMap"].is_null()) {
            auto textureMapArray = textures["TextureMap"];
            if (textureMapArray.is_array()) {
                for (auto textureMapData : textureMapArray) {
                    try {
                        TextureMap newTextureMap = parseTextureMap(textureMapData);
                        if (newTextureMap._id != 0) {
                            this->textureMaps.push_back(newTextureMap);
                            this->textureMapIdToIndex[newTextureMap._id] = this->textureMaps.size() - 1;
                            if (newTextureMap.decalMode == DecalMode::ReplaceBackground) {
                                // Only accept if it points to a valid loaded image (for image textures)
                                bool ok = true;
                                if (newTextureMap.type == "image") {
                                    ok = (this->getImageById(newTextureMap.imageId) != nullptr);
                                }
                                if (ok) {
                                    this->backgroundTextureId = newTextureMap._id;
                                    verbose("[+] BackgroundTexture set from TextureMap " + std::to_string(newTextureMap._id) +
                                            " (imageId=" + std::to_string(newTextureMap.imageId) + ")");
                                } else {
                                    verbose("[!] Ignoring replace_background TextureMap " + std::to_string(newTextureMap._id) +
                                            " because ImageId " + std::to_string(newTextureMap.imageId) + " was not loaded");
                                }
                            }
                            verbose("[+] TextureMap Parsed: " + std::to_string(newTextureMap._id));
                        }
                    } catch (const std::exception& e) {
                        verbose("[!] Skipping TextureMap Parsing. Reason: " + std::string(e.what()));
                    }
                }
            } else {
                try {
                    TextureMap newTextureMap = parseTextureMap(textureMapArray);
                    if (newTextureMap._id != 0) {
                        this->textureMaps.push_back(newTextureMap);
                        this->textureMapIdToIndex[newTextureMap._id] = this->textureMaps.size() - 1;
                        if (newTextureMap.decalMode == DecalMode::ReplaceBackground) {
                            bool ok = true;
                            if (newTextureMap.type == "image") {
                                ok = (this->getImageById(newTextureMap.imageId) != nullptr);
                            }
                            if (ok) {
                                this->backgroundTextureId = newTextureMap._id;
                                verbose("[+] BackgroundTexture set from TextureMap " + std::to_string(newTextureMap._id) +
                                        " (imageId=" + std::to_string(newTextureMap.imageId) + ")");
                            } else {
                                verbose("[!] Ignoring replace_background TextureMap " + std::to_string(newTextureMap._id) +
                                        " because ImageId " + std::to_string(newTextureMap.imageId) + " was not loaded");
                            }
                        }
                        verbose("[+] TextureMap Parsed: " + std::to_string(newTextureMap._id));
                    }
                } catch (const std::exception& e) {
                    verbose("[!] Skipping TextureMap Parsing. Reason: " + std::string(e.what()));
                }
            }
        }
    }

    // Fallback: if no explicit BackgroundTexture field and no ReplaceBackground texture was seen while parsing,
    // scan all texture maps and pick the first ReplaceBackground.
    if (this->backgroundTextureId == 0) {
        for (const auto& tm : this->textureMaps) {
            if (tm.decalMode == DecalMode::ReplaceBackground) {
                this->backgroundTextureId = tm._id;
                verbose("[+] BackgroundTexture fallback set to TextureMap " + std::to_string(tm._id));
                break;
            }
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
        
        // Parse LightSphere (emissive spheres)
        if (objects.contains("LightSphere") && !objects["LightSphere"].is_null()) {
            auto lightSphereArray = objects["LightSphere"];
            if (lightSphereArray.is_array()) {
                for (auto lightSphereData : lightSphereArray) {
                    scene::LightSphere newLightSphere;
                    newLightSphere._id = parseSingleValue<unsigned int>(lightSphereData["_id"]);
                    
                    // Parse Material if present
                    if (lightSphereData.contains("Material") && !lightSphereData["Material"].is_null()) {
                        unsigned int materialId = parseSingleValue<unsigned int>(lightSphereData["Material"]);
                        newLightSphere.material = getMaterialById(materialId);
                    } else {
                        newLightSphere.material = nullptr;
                    }
                    
                    // Parse Radiance
                    if (lightSphereData.contains("Radiance") && !lightSphereData["Radiance"].is_null()) {
                        newLightSphere.radiance = parseTriplet<VectorFloatTriplet>(lightSphereData["Radiance"]);
                    } else {
                        newLightSphere.radiance = VectorFloatTriplet{0, 0, 0};
                    }
                    
                    // Parse Center and Radius (same as regular Sphere)
                    if (lightSphereData.contains("Center")) {
                        newLightSphere.center = parseSingleValue<unsigned int>(lightSphereData["Center"]);
                        newLightSphere.center -= 1;  // Convert to 0-based
                    }
                    if (lightSphereData.contains("Radius")) {
                        newLightSphere.radius = parseSingleValue<double>(lightSphereData["Radius"]);
                    }
                    
                    // Parse transformations if present
                    if (lightSphereData.contains("Transformations") && !lightSphereData["Transformations"].is_null()) {
                        std::string transformStr = lightSphereData["Transformations"].get<std::string>();
                        newLightSphere.transformations = parseTransformationString(transformStr);
                        verbose("[+] LightSphere Transformations parsed: " + transformStr);
                    }
                    
                    // Parse Textures if present
                    if (lightSphereData.contains("Textures") && !lightSphereData["Textures"].is_null()) {
                        std::string texturesStr = lightSphereData["Textures"].get<std::string>();
                        std::istringstream stream(texturesStr);
                        unsigned int textureId;
                        while (stream >> textureId) {
                            newLightSphere.textureIds.push_back(textureId);
                        }
                    }
                    
                    this->lightSpheres.push_back(newLightSphere);
                    verbose("[+] LightSphere Parsed: " + std::to_string(newLightSphere._id));
                }
            } else {
                scene::LightSphere newLightSphere;
                newLightSphere._id = parseSingleValue<unsigned int>(lightSphereArray["_id"]);
                
                if (lightSphereArray.contains("Material") && !lightSphereArray["Material"].is_null()) {
                    unsigned int materialId = parseSingleValue<unsigned int>(lightSphereArray["Material"]);
                    newLightSphere.material = getMaterialById(materialId);
                } else {
                    newLightSphere.material = nullptr;
                }
                
                if (lightSphereArray.contains("Radiance") && !lightSphereArray["Radiance"].is_null()) {
                    newLightSphere.radiance = parseTriplet<VectorFloatTriplet>(lightSphereArray["Radiance"]);
                } else {
                    newLightSphere.radiance = VectorFloatTriplet{0, 0, 0};
                }
                
                if (lightSphereArray.contains("Center")) {
                    newLightSphere.center = parseSingleValue<unsigned int>(lightSphereArray["Center"]);
                    newLightSphere.center -= 1;
                }
                if (lightSphereArray.contains("Radius")) {
                    newLightSphere.radius = parseSingleValue<double>(lightSphereArray["Radius"]);
                }
                
                if (lightSphereArray.contains("Transformations") && !lightSphereArray["Transformations"].is_null()) {
                    std::string transformStr = lightSphereArray["Transformations"].get<std::string>();
                    newLightSphere.transformations = parseTransformationString(transformStr);
                }
                
                if (lightSphereArray.contains("Textures") && !lightSphereArray["Textures"].is_null()) {
                    std::string texturesStr = lightSphereArray["Textures"].get<std::string>();
                    std::istringstream stream(texturesStr);
                    unsigned int textureId;
                    while (stream >> textureId) {
                        newLightSphere.textureIds.push_back(textureId);
                    }
                }
                
                this->lightSpheres.push_back(newLightSphere);
                verbose("[+] LightSphere Parsed: " + std::to_string(newLightSphere._id));
            }
        }
        
        // Parse LightMesh (emissive meshes)
        if (objects.contains("LightMesh") && !objects["LightMesh"].is_null()) {
            auto lightMeshArray = objects["LightMesh"];
            if (lightMeshArray.is_array()) {
                for (auto lightMeshData : lightMeshArray) {
                    scene::LightMesh newLightMesh;
                    newLightMesh._id = parseSingleValue<unsigned int>(lightMeshData["_id"]);
                    
                    // Parse Material if present
                    if (lightMeshData.contains("Material") && !lightMeshData["Material"].is_null()) {
                        unsigned int materialId = parseSingleValue<unsigned int>(lightMeshData["Material"]);
                        newLightMesh.material = getMaterialById(materialId);
                    } else {
                        newLightMesh.material = nullptr;
                    }
                    
                    // Parse Radiance
                    if (lightMeshData.contains("Radiance") && !lightMeshData["Radiance"].is_null()) {
                        newLightMesh.radiance = parseTriplet<VectorFloatTriplet>(lightMeshData["Radiance"]);
                    } else {
                        newLightMesh.radiance = VectorFloatTriplet{0, 0, 0};
                    }
                    
                    // Parse shading mode
                    if (lightMeshData.contains("_shadingMode") && !lightMeshData["_shadingMode"].is_null()) {
                        newLightMesh.shadingMode = tolower(lightMeshData["_shadingMode"].get<std::string>()[0]);
                    }
                    
                    // Parse Faces (same as regular Mesh)
                    if (lightMeshData.contains("Faces") && !lightMeshData["Faces"].is_null()) {
                        const json& facesData = lightMeshData["Faces"];
                        if (facesData.contains("_data") && !facesData["_data"].is_null()) {
                            FaceParseResult faceResult = parseFacesWithOffsets(facesData);
                            newLightMesh.faces = faceResult.vertexFaces;
                            newLightMesh.texCoordIndices = faceResult.texCoordFaces;
                        } else if (facesData.contains("_plyFile") && !facesData["_plyFile"].is_null()) {
                            std::string plyFile = facesData["_plyFile"].get<std::string>();
                            std::string fullPath = this->baseDirectory + plyFile;
                            FaceParseResult faceResult = parsePLYFile(fullPath, this->vertices, this->texCoords);
                            newLightMesh.faces = faceResult.vertexFaces;
                            newLightMesh.texCoordIndices = faceResult.texCoordFaces;
                        }
                    }
                    
                    // Parse transformations if present
                    if (lightMeshData.contains("Transformations") && !lightMeshData["Transformations"].is_null()) {
                        std::string transformStr = lightMeshData["Transformations"].get<std::string>();
                        newLightMesh.transformations = parseTransformationString(transformStr);
                        verbose("[+] LightMesh Transformations parsed: " + transformStr);
                    }
                    
                    // Parse Textures if present
                    if (lightMeshData.contains("Textures") && !lightMeshData["Textures"].is_null()) {
                        std::string texturesStr = lightMeshData["Textures"].get<std::string>();
                        std::istringstream stream(texturesStr);
                        unsigned int textureId;
                        while (stream >> textureId) {
                            newLightMesh.textureIds.push_back(textureId);
                        }
                    }
                    
                    this->lightMeshes.push_back(newLightMesh);
                    verbose("[+] LightMesh Parsed: " + std::to_string(newLightMesh._id));
                }
            } else {
                scene::LightMesh newLightMesh;
                newLightMesh._id = parseSingleValue<unsigned int>(lightMeshArray["_id"]);
                
                if (lightMeshArray.contains("Material") && !lightMeshArray["Material"].is_null()) {
                    unsigned int materialId = parseSingleValue<unsigned int>(lightMeshArray["Material"]);
                    newLightMesh.material = getMaterialById(materialId);
                } else {
                    newLightMesh.material = nullptr;
                }
                
                if (lightMeshArray.contains("Radiance") && !lightMeshArray["Radiance"].is_null()) {
                    newLightMesh.radiance = parseTriplet<VectorFloatTriplet>(lightMeshArray["Radiance"]);
                } else {
                    newLightMesh.radiance = VectorFloatTriplet{0, 0, 0};
                }
                
                if (lightMeshArray.contains("_shadingMode") && !lightMeshArray["_shadingMode"].is_null()) {
                    newLightMesh.shadingMode = tolower(lightMeshArray["_shadingMode"].get<std::string>()[0]);
                }
                
                if (lightMeshArray.contains("Faces") && !lightMeshArray["Faces"].is_null()) {
                    const json& facesData = lightMeshArray["Faces"];
                    if (facesData.contains("_data") && !facesData["_data"].is_null()) {
                        FaceParseResult faceResult = parseFacesWithOffsets(facesData);
                        newLightMesh.faces = faceResult.vertexFaces;
                        newLightMesh.texCoordIndices = faceResult.texCoordFaces;
                    } else if (facesData.contains("_plyFile") && !facesData["_plyFile"].is_null()) {
                        std::string plyFile = facesData["_plyFile"].get<std::string>();
                        std::string fullPath = this->baseDirectory + plyFile;
                        FaceParseResult faceResult = parsePLYFile(fullPath, this->vertices, this->texCoords);
                        newLightMesh.faces = faceResult.vertexFaces;
                        newLightMesh.texCoordIndices = faceResult.texCoordFaces;
                    }
                }
                
                if (lightMeshArray.contains("Transformations") && !lightMeshArray["Transformations"].is_null()) {
                    std::string transformStr = lightMeshArray["Transformations"].get<std::string>();
                    newLightMesh.transformations = parseTransformationString(transformStr);
                }
                
                if (lightMeshArray.contains("Textures") && !lightMeshArray["Textures"].is_null()) {
                    std::string texturesStr = lightMeshArray["Textures"].get<std::string>();
                    std::istringstream stream(texturesStr);
                    unsigned int textureId;
                    while (stream >> textureId) {
                        newLightMesh.textureIds.push_back(textureId);
                    }
                }
                
                this->lightMeshes.push_back(newLightMesh);
                verbose("[+] LightMesh Parsed: " + std::to_string(newLightMesh._id));
            }
        }
        
        verbose("[+] LightSpheres Parsed: " + std::to_string(this->lightSpheres.size()));
        verbose("[+] LightMeshes Parsed: " + std::to_string(this->lightMeshes.size()));
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
            if (cameraData.contains("_id") && !cameraData["_id"].is_null()) {
                newCamera._id = parseSingleValue<unsigned int>(cameraData["_id"]);
            }
            if (cameraData.contains("Position") && !cameraData["Position"].is_null()) {
                newCamera.position = parseTriplet<VectorFloatTriplet>(cameraData["Position"]);
            }
            if (cameraData.contains("Up") && !cameraData["Up"].is_null()) {
                newCamera.up = parseTriplet<VectorFloatTriplet>(cameraData["Up"]);
            }
            if (cameraData.contains("NearDistance") && !cameraData["NearDistance"].is_null()) {
                newCamera.nearDistance = parseSingleValue<double>(cameraData["NearDistance"]);
            }
            if (cameraData.contains("ImageResolution") && !cameraData["ImageResolution"].is_null()) {
                newCamera.imageResolution = parsePair<VectorIntPair>(cameraData["ImageResolution"]);
            }
            if (cameraData.contains("ImageName") && !cameraData["ImageName"].is_null()) {
                newCamera.imageName = cameraData["ImageName"].get<std::string>();
            }
            // Calculate the gaze vector, nearplane
            // Gaze: GazePoint - Position OR directly from Gaze field
            if (cameraData.contains("GazePoint") && !cameraData["GazePoint"].is_null()) {
                VectorFloatTriplet gazePoint = parseTriplet<VectorFloatTriplet>(cameraData["GazePoint"]);
                VectorFloatTriplet gazeVec = gazePoint - newCamera.position;
                double gazeLen = std::sqrt(dotProduct(gazeVec, gazeVec));
                if (gazeLen > 1e-10) {
                    newCamera.gaze = gazeVec * (1.0 / gazeLen);
                } else {
                    // Fallback: use default gaze direction
                    newCamera.gaze = VectorFloatTriplet{0, 0, -1};
                    verbose("[!] GazePoint equals Position, using default gaze direction");
                }
            } else if (cameraData.contains("Gaze") && !cameraData["Gaze"].is_null()) {
                // Direct gaze vector specified
                newCamera.gaze = parseTriplet<VectorFloatTriplet>(cameraData["Gaze"]);
                // Normalize the gaze vector
                double gazeLen = std::sqrt(dotProduct(newCamera.gaze, newCamera.gaze));
                if (gazeLen > 1e-10) {
                    newCamera.gaze = newCamera.gaze * (1.0 / gazeLen);
                } else {
                    newCamera.gaze = VectorFloatTriplet{0, 0, -1};
                    verbose("[!] Gaze vector is zero, using default gaze direction");
                }
            }

            // NearPlane:
            // FovY is the vertical field of view in degrees. We compute the near-plane
            // height from it, then use the aspect ratio to get the width.
            double fovY = 45.0;  // Default
            if (cameraData.contains("FovY") && !cameraData["FovY"].is_null()) {
                fovY = parseSingleValue<double>(cameraData["FovY"]);
            }
            fovY = fovY * M_PI / 180.0; // convert degrees to radians
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
            if (cameraData.contains("Position") && !cameraData["Position"].is_null()) {
                newCamera.position = parseTriplet<VectorFloatTriplet>(cameraData["Position"]);
            }
            if (cameraData.contains("Gaze") && !cameraData["Gaze"].is_null()) {
                newCamera.gaze = parseTriplet<VectorFloatTriplet>(cameraData["Gaze"]);
            }
            if (cameraData.contains("Up") && !cameraData["Up"].is_null()) {
                newCamera.up = parseTriplet<VectorFloatTriplet>(cameraData["Up"]);
            }
            if (cameraData.contains("NearPlane") && !cameraData["NearPlane"].is_null()) {
                newCamera.nearPlane = parseQuad<VectorFloatQuad>(cameraData["NearPlane"]);
            }
            if (cameraData.contains("NearDistance") && !cameraData["NearDistance"].is_null()) {
                newCamera.nearDistance = parseSingleValue<double>(cameraData["NearDistance"]);
            }
            if (cameraData.contains("ImageResolution") && !cameraData["ImageResolution"].is_null()) {
                newCamera.imageResolution = parsePair<VectorIntPair>(cameraData["ImageResolution"]);
            }
            if (cameraData.contains("ImageName") && !cameraData["ImageName"].is_null()) {
                newCamera.imageName = cameraData["ImageName"].get<std::string>();
            }
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
    
    // Parse NumSamples if present
    if (cameraData.contains("NumSamples") && !cameraData["NumSamples"].is_null()) {
        newCamera.numSamples = parseSingleValue<int>(cameraData["NumSamples"]);
        verbose("[+] Camera NumSamples parsed: " + std::to_string(newCamera.numSamples));
    } else {
        newCamera.numSamples = 1;
        verbose("[!] Camera NumSamples not found, using default: 1");
    }
    
    // Parse ApertureSize if present (optional, enables depth-of-field)
    if (cameraData.contains("ApertureSize") && !cameraData["ApertureSize"].is_null()) {
        newCamera.apertureSize = parseSingleValue<double>(cameraData["ApertureSize"]);
        verbose("[+] Camera ApertureSize parsed: " + std::to_string(newCamera.apertureSize));
    } else {
        newCamera.apertureSize = 0.0;
        verbose("[!] Camera ApertureSize not found, depth-of-field disabled");
    }
    
    // Parse FocusDistance if present (used with ApertureSize)
    if (cameraData.contains("FocusDistance") && !cameraData["FocusDistance"].is_null()) {
        newCamera.focusDistance = parseSingleValue<double>(cameraData["FocusDistance"]);
        verbose("[+] Camera FocusDistance parsed: " + std::to_string(newCamera.focusDistance));
    } else {
        newCamera.focusDistance = 0.0;
        verbose("[!] Camera FocusDistance not found, using default: 0.0");
    }
    
    // Parse Tonemap if present (can be object or array)
    if (cameraData.contains("Tonemap") && !cameraData["Tonemap"].is_null()) {
        auto tonemapData = cameraData["Tonemap"];
        if (tonemapData.is_array()) {
            for (auto tonemapEntry : tonemapData) {
                TonemapSettings settings;
                settings.tmo = tonemapEntry["TMO"].get<std::string>();
                settings.tmoOptions = tonemapEntry["TMOOptions"].get<std::string>();
                settings.saturation = parseSingleValue<double>(tonemapEntry["Saturation"]);
                settings.gamma = parseSingleValue<double>(tonemapEntry["Gamma"]);
                settings.extension = tonemapEntry["Extension"].get<std::string>();
                newCamera.tonemapSettings.push_back(settings);
                verbose("[+] Camera Tonemap parsed: " + settings.tmo + " with extension " + settings.extension);
            }
        } else {
            TonemapSettings settings;
            settings.tmo = tonemapData["TMO"].get<std::string>();
            settings.tmoOptions = tonemapData["TMOOptions"].get<std::string>();
            settings.saturation = parseSingleValue<double>(tonemapData["Saturation"]);
            settings.gamma = parseSingleValue<double>(tonemapData["Gamma"]);
            settings.extension = tonemapData["Extension"].get<std::string>();
            newCamera.tonemapSettings.push_back(settings);
            verbose("[+] Camera Tonemap parsed: " + settings.tmo + " with extension " + settings.extension);
        }
    }
    
    // Parse Renderer if present (PathTracing or empty for default)
    if (cameraData.contains("Renderer") && !cameraData["Renderer"].is_null()) {
        newCamera.renderer = cameraData["Renderer"].get<std::string>();
        verbose("[+] Camera Renderer parsed: " + newCamera.renderer);
    } else {
        newCamera.renderer = "";
        verbose("[!] Camera Renderer not found, using default shading");
    }
    
    // Parse RendererParams if present (space-separated options)
    if (cameraData.contains("RendererParams") && !cameraData["RendererParams"].is_null()) {
        std::string paramsStr = cameraData["RendererParams"].get<std::string>();
        std::istringstream paramsStream(paramsStr);
        std::string param;
        while (paramsStream >> param) {
            if (param == "ImportanceSampling") {
                newCamera.importanceSampling = true;
                verbose("[+] Camera RendererParams: ImportanceSampling enabled");
            } else if (param == "NextEventEstimation") {
                newCamera.nextEventEstimation = true;
                verbose("[+] Camera RendererParams: NextEventEstimation enabled");
            } else if (param == "MIS_BALANCE") {
                newCamera.misHeuristic = "balance";
                verbose("[+] Camera RendererParams: MIS_BALANCE (balance heuristic)");
            } else if (param == "MIS_POWER") {
                newCamera.misHeuristic = "power";
                verbose("[+] Camera RendererParams: MIS_POWER (power heuristic)");
            } else if (param == "MIS_01") {
                newCamera.misHeuristic = "01";
                verbose("[+] Camera RendererParams: MIS_01 (01 heuristic)");
            } else if (param == "RussianRoulette") {
                newCamera.russianRoulette = true;
                verbose("[+] Camera RendererParams: RussianRoulette enabled");
            }
        }
    }
    
    // Parse MaxRecursionDepth if present (per-camera override)
    if (cameraData.contains("MaxRecursionDepth") && !cameraData["MaxRecursionDepth"].is_null()) {
        newCamera.maxRecursionDepth = parseSingleValue<int>(cameraData["MaxRecursionDepth"]);
        verbose("[+] Camera MaxRecursionDepth parsed: " + std::to_string(newCamera.maxRecursionDepth));
    } else {
        // Use scene default if not specified
        verbose("[!] Camera MaxRecursionDepth not found, will use scene default");
    }
    
    // Parse MinRecursionDepth if present (for Russian Roulette)
    if (cameraData.contains("MinRecursionDepth") && !cameraData["MinRecursionDepth"].is_null()) {
        newCamera.minRecursionDepth = parseSingleValue<int>(cameraData["MinRecursionDepth"]);
        verbose("[+] Camera MinRecursionDepth parsed: " + std::to_string(newCamera.minRecursionDepth));
    } else {
        newCamera.minRecursionDepth = 0;
        verbose("[!] Camera MinRecursionDepth not found, using default: 0");
    }
    
    // Parse SplittingFactor if present
    if (cameraData.contains("SplittingFactor") && !cameraData["SplittingFactor"].is_null()) {
        newCamera.splittingFactor = parseSingleValue<int>(cameraData["SplittingFactor"]);
        verbose("[+] Camera SplittingFactor parsed: " + std::to_string(newCamera.splittingFactor));
    } else {
        newCamera.splittingFactor = 1;
        verbose("[!] Camera SplittingFactor not found, using default: 1");
    }
    
    // Parse SampleMaxVal if present (clamping threshold)
    if (cameraData.contains("SampleMaxVal") && !cameraData["SampleMaxVal"].is_null()) {
        newCamera.sampleMaxVal = parseSingleValue<double>(cameraData["SampleMaxVal"]);
        verbose("[+] Camera SampleMaxVal parsed: " + std::to_string(newCamera.sampleMaxVal));
    } else {
        newCamera.sampleMaxVal = 0.0;
        verbose("[!] Camera SampleMaxVal not found, no clamping (0.0)");
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

scene::AreaLight scene::parseAreaLight(const json& areaLightData) {
    scene::AreaLight newAreaLight;
    newAreaLight._id = parseSingleValue<unsigned int>(areaLightData["_id"]);
    newAreaLight.position = parseTriplet<VectorFloatTriplet>(areaLightData["Position"]);
    newAreaLight.normal = parseTriplet<VectorFloatTriplet>(areaLightData["Normal"]);
    newAreaLight.size = parseSingleValue<double>(areaLightData["Size"]);
    newAreaLight.radiance = parseTriplet<VectorFloatTriplet>(areaLightData["Radiance"]);
    
    // Parse transformations if present
    if (areaLightData.contains("Transformations") && !areaLightData["Transformations"].is_null()) {
        std::string transformStr = areaLightData["Transformations"].get<std::string>();
        newAreaLight.transformations = parseTransformationString(transformStr);
        verbose("[+] AreaLight Transformations parsed: " + transformStr + " (" + std::to_string(newAreaLight.transformations.size()) + " transforms)");
    }
    
    return newAreaLight;
}

scene::DirectionalLight scene::parseDirectionalLight(const json& directionalLightData) {
    scene::DirectionalLight newDirectionalLight;
    newDirectionalLight._id = parseSingleValue<unsigned int>(directionalLightData["_id"]);
    newDirectionalLight.direction = parseTriplet<VectorFloatTriplet>(directionalLightData["Direction"]);
    newDirectionalLight.radiance = parseTriplet<VectorFloatTriplet>(directionalLightData["Radiance"]);
    
    // Parse transformations if present
    if (directionalLightData.contains("Transformations") && !directionalLightData["Transformations"].is_null()) {
        std::string transformStr = directionalLightData["Transformations"].get<std::string>();
        newDirectionalLight.transformations = parseTransformationString(transformStr);
        verbose("[+] DirectionalLight Transformations parsed: " + transformStr + " (" + std::to_string(newDirectionalLight.transformations.size()) + " transforms)");
    }
    
    return newDirectionalLight;
}

scene::SpotLight scene::parseSpotLight(const json& spotLightData) {
    scene::SpotLight newSpotLight;
    newSpotLight._id = parseSingleValue<unsigned int>(spotLightData["_id"]);
    newSpotLight.position = parseTriplet<VectorFloatTriplet>(spotLightData["Position"]);
    newSpotLight.direction = parseTriplet<VectorFloatTriplet>(spotLightData["Direction"]);
    newSpotLight.intensity = parseTriplet<VectorFloatTriplet>(spotLightData["Intensity"]);
    newSpotLight.coverageAngle = parseSingleValue<double>(spotLightData["CoverageAngle"]);
    newSpotLight.falloffAngle = parseSingleValue<double>(spotLightData["FalloffAngle"]);
    
    // Parse transformations if present
    if (spotLightData.contains("Transformations") && !spotLightData["Transformations"].is_null()) {
        std::string transformStr = spotLightData["Transformations"].get<std::string>();
        newSpotLight.transformations = parseTransformationString(transformStr);
        verbose("[+] SpotLight Transformations parsed: " + transformStr + " (" + std::to_string(newSpotLight.transformations.size()) + " transforms)");
    }
    
    return newSpotLight;
}

scene::SphericalDirectionalLight scene::parseSphericalDirectionalLight(const json& sphericalDirectionalLightData) {
    scene::SphericalDirectionalLight newSphericalDirectionalLight;
    newSphericalDirectionalLight._id = parseSingleValue<unsigned int>(sphericalDirectionalLightData["_id"]);
    newSphericalDirectionalLight.imageId = parseSingleValue<unsigned int>(sphericalDirectionalLightData["ImageId"]);
    
    if (sphericalDirectionalLightData.contains("_type") && !sphericalDirectionalLightData["_type"].is_null()) {
        newSphericalDirectionalLight.type = sphericalDirectionalLightData["_type"].get<std::string>();
    } else {
        newSphericalDirectionalLight.type = "latlong";  // Default
    }
    
    if (sphericalDirectionalLightData.contains("Sampler") && !sphericalDirectionalLightData["Sampler"].is_null()) {
        newSphericalDirectionalLight.sampler = sphericalDirectionalLightData["Sampler"].get<std::string>();
    } else {
        newSphericalDirectionalLight.sampler = "cosine";  // Default
    }
    
    // Parse transformations if present
    if (sphericalDirectionalLightData.contains("Transformations") && !sphericalDirectionalLightData["Transformations"].is_null()) {
        std::string transformStr = sphericalDirectionalLightData["Transformations"].get<std::string>();
        newSphericalDirectionalLight.transformations = parseTransformationString(transformStr);
        verbose("[+] SphericalDirectionalLight Transformations parsed: " + transformStr + " (" + std::to_string(newSphericalDirectionalLight.transformations.size()) + " transforms)");
    }
    
    return newSphericalDirectionalLight;
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
        // disabled this because we are only setting mirrors for type "mirror"
        if (newMaterial.mirrorReflectance.x > 0.0 || newMaterial.mirrorReflectance.y > 0.0 || newMaterial.mirrorReflectance.z > 0.0)  {
            // newMaterial.isMirror = true;
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
    
    // Parse Roughness if present (for mirrors, conductors, and dielectrics)
    if (materialData.contains("Roughness") && !materialData["Roughness"].is_null()) {
        newMaterial.roughness = parseSingleValue<double>(materialData["Roughness"]);
        verbose("[+] Material Roughness parsed: " + std::to_string(newMaterial.roughness));
    } else {
        newMaterial.roughness = 0.0;
    }
    
    // Parse _BRDF if present (reference to BRDF definition)
    if (materialData.contains("_BRDF") && !materialData["_BRDF"].is_null()) {
        newMaterial.brdfId = parseSingleValue<unsigned int>(materialData["_BRDF"]);
        verbose("[+] Material _BRDF parsed: " + std::to_string(newMaterial.brdfId));
    } else {
        newMaterial.brdfId = 0;  // Default: use OriginalBlinnPhong
    }
    
    return newMaterial;
}

scene::BRDF scene::parseBRDF(const json& brdfData, scene::BRDFType type) {
    scene::BRDF newBRDF;
    newBRDF._id = parseSingleValue<unsigned int>(brdfData["_id"]);
    newBRDF.type = type;
    
    // Parse Exponent
    if (brdfData.contains("Exponent") && !brdfData["Exponent"].is_null()) {
        newBRDF.exponent = parseSingleValue<double>(brdfData["Exponent"]);
    } else {
        newBRDF.exponent = 1.0;
    }
    
    // Parse _normalized flag
    if (brdfData.contains("_normalized") && !brdfData["_normalized"].is_null()) {
        std::string normalizedStr = brdfData["_normalized"].get<std::string>();
        newBRDF.normalized = (normalizedStr == "true" || normalizedStr == "True" || normalizedStr == "TRUE" || normalizedStr == "1");
    } else {
        newBRDF.normalized = false;
    }
    
    // Parse kdfresnel for TorranceSparrow
    if (type == scene::BRDFType::TorranceSparrow && brdfData.contains("kdfresnel") && !brdfData["kdfresnel"].is_null()) {
        std::string kdfresnelStr = brdfData["kdfresnel"].get<std::string>();
        newBRDF.kdFresnel = (kdfresnelStr == "true" || kdfresnelStr == "True" || kdfresnelStr == "TRUE" || kdfresnelStr == "1");
    } else {
        newBRDF.kdFresnel = false;
    }
    
    return newBRDF;
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

scene::Scene::FaceParseResult scene::Scene::parseFacesWithOffsets(const json& facesData) {
    FaceParseResult result;
    
    if (facesData.contains("_data")) {
        auto facesDataArray = facesData["_data"];
        std::stringstream stream(facesDataArray.get<std::string>());
        
        // Parse vertex and texture offsets
        // Offsets are interpreted as a shift applied after converting to 0-based indexing:
        //
        //   finalIndex = rawIndex + offset - 1
        //
        // This convention supports:
        // - **1-based faces with no shift**: offset=0 -> final = raw - 1
        // - **0-based faces with no shift**: offset=1 -> final = raw
        // - **large positive vertex offsets** (e.g. appended vertex blocks): offset=46352
        // - **negative texture offsets** (e.g. mapping into a small TexCoordData array): offset=-4, -8, ...
        int vertexOffset = 0;
        int textureOffset = 0;
        
        if (facesData.contains("_vertexOffset") && !facesData["_vertexOffset"].is_null()) {
            vertexOffset = parseSingleValue<int>(facesData["_vertexOffset"]);
            verbose("[+] Faces _vertexOffset: " + std::to_string(vertexOffset));
        }
        
        if (facesData.contains("_textureOffset") && !facesData["_textureOffset"].is_null()) {
            textureOffset = parseSingleValue<int>(facesData["_textureOffset"]);
            verbose("[+] Faces _textureOffset: " + std::to_string(textureOffset));
        }
        
        // Read and convert faces
        int rawX, rawY, rawZ;  // Store original face values from file
        scene::VectorIntTriplet vertexFace, texCoordFace;
        while (stream >> rawX >> rawY >> rawZ) {
            vertexFace.x = rawX + vertexOffset - 1;
            vertexFace.y = rawY + vertexOffset - 1;
            vertexFace.z = rawZ + vertexOffset - 1;
            result.vertexFaces.push_back(vertexFace);
            
            texCoordFace.x = rawX + textureOffset - 1;
            texCoordFace.y = rawY + textureOffset - 1;
            texCoordFace.z = rawZ + textureOffset - 1;
            result.texCoordFaces.push_back(texCoordFace);
        }
        stream.clear();
        return result;
    }
    
    return result;
}

std::vector<scene::VectorIntTriplet> scene::Scene::parseFaces(const json& facesData) {
    if (facesData.contains("_data")) {
        FaceParseResult faceResult = parseFacesWithOffsets(facesData);
        return faceResult.vertexFaces;
    } else if (facesData.contains("_plyFile")) {
        std::string plyFile = facesData["_plyFile"].get<std::string>();
        std::string fullPath = this->baseDirectory + plyFile;
        // parsePLYFile will add vertices and texture coordinates, and return adjusted faces
        FaceParseResult faceResult = parsePLYFile(fullPath, this->vertices, this->texCoords);
        return faceResult.vertexFaces;
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
            
            // Parse Textures if present
            if (objectData.contains("Textures") && !objectData["Textures"].is_null()) {
                std::string texturesStr = objectData["Textures"].get<std::string>();
                std::istringstream stream(texturesStr);
                unsigned int textureId;
                while (stream >> textureId) {
                    newObject.textureIds.push_back(textureId);
                }
                verbose("[+] Object Textures parsed: " + texturesStr + " (" + std::to_string(newObject.textureIds.size()) + " textures)");
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
        
        // Parse Textures if present
        if (objectsData.contains("Textures") && !objectsData["Textures"].is_null()) {
            std::string texturesStr = objectsData["Textures"].get<std::string>();
            std::istringstream stream(texturesStr);
            unsigned int textureId;
            while (stream >> textureId) {
                newObject.textureIds.push_back(textureId);
            }
            verbose("[+] Object Textures parsed: " + texturesStr + " (" + std::to_string(newObject.textureIds.size()) + " textures)");
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
        const json& facesData = objectData["Faces"];
        if (facesData.contains("_data") && !facesData["_data"].is_null()) {
            FaceParseResult faceResult = parseFacesWithOffsets(facesData);
            object.faces = faceResult.vertexFaces;
            object.texCoordIndices = faceResult.texCoordFaces;
        } else if (facesData.contains("_plyFile") && !facesData["_plyFile"].is_null()) {
            std::string plyFile = facesData["_plyFile"].get<std::string>();
            std::string fullPath = this->baseDirectory + plyFile;
            FaceParseResult faceResult = parsePLYFile(fullPath, this->vertices, this->texCoords);
            object.faces = faceResult.vertexFaces;
            object.texCoordIndices = faceResult.texCoordFaces;
        } else if (facesData.contains("_binaryFile") && !facesData["_binaryFile"].is_null()) {
            // Parse binary face data: 4-byte count + count*12 bytes (3 uint32 per triangle)
            std::string binaryFile = facesData["_binaryFile"].get<std::string>();
            std::string fullPath = this->baseDirectory + binaryFile;
            std::ifstream file(fullPath, std::ios::binary);
            if (file.is_open()) {
                uint32_t count;
                file.read(reinterpret_cast<char*>(&count), sizeof(count));
                object.faces.reserve(count);
                for (uint32_t i = 0; i < count; i++) {
                    uint32_t v0, v1, v2;
                    file.read(reinterpret_cast<char*>(&v0), sizeof(v0));
                    file.read(reinterpret_cast<char*>(&v1), sizeof(v1));
                    file.read(reinterpret_cast<char*>(&v2), sizeof(v2));
                    // Binary files are typically 0-indexed, no offset needed
                    object.faces.push_back(VectorIntTriplet{(int)v0, (int)v1, (int)v2});
                    // Also set texture coordinates if available (same indices for now)
                    object.texCoordIndices.push_back(VectorIntTriplet{(int)v0, (int)v1, (int)v2});
                }
                file.close();
            }
        }
    }
    
    // Parse MotionBlur if present
    if (objectData.contains("MotionBlur") && !objectData["MotionBlur"].is_null()) {
        object.motionBlur = parseTriplet<VectorFloatTriplet>(objectData["MotionBlur"]);
        object.hasMotionBlur = true;
        verbose("[+] Mesh MotionBlur parsed: " + std::to_string(object.motionBlur.x) + " " + std::to_string(object.motionBlur.y) + " " + std::to_string(object.motionBlur.z));
    } else {
        object.motionBlur = {0, 0, 0};
        object.hasMotionBlur = false;
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
    
    // Parse MotionBlur if present
    if (objectData.contains("MotionBlur") && !objectData["MotionBlur"].is_null()) {
        object.motionBlur = parseTriplet<VectorFloatTriplet>(objectData["MotionBlur"]);
        object.hasMotionBlur = true;
        verbose("[+] Triangle MotionBlur parsed: " + std::to_string(object.motionBlur.x) + " " + std::to_string(object.motionBlur.y) + " " + std::to_string(object.motionBlur.z));
    } else {
        object.motionBlur = {0, 0, 0};
        object.hasMotionBlur = false;
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
    
    // Parse MotionBlur if present
    if (objectData.contains("MotionBlur") && !objectData["MotionBlur"].is_null()) {
        object.motionBlur = parseTriplet<VectorFloatTriplet>(objectData["MotionBlur"]);
        object.hasMotionBlur = true;
        verbose("[+] Sphere MotionBlur parsed: " + std::to_string(object.motionBlur.x) + " " + std::to_string(object.motionBlur.y) + " " + std::to_string(object.motionBlur.z));
    } else {
        object.motionBlur = {0, 0, 0};
        object.hasMotionBlur = false;
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
    
    // Parse MotionBlur if present
    if (objectData.contains("MotionBlur") && !objectData["MotionBlur"].is_null()) {
        object.motionBlur = parseTriplet<VectorFloatTriplet>(objectData["MotionBlur"]);
        object.hasMotionBlur = true;
        verbose("[+] Plane MotionBlur parsed: " + std::to_string(object.motionBlur.x) + " " + std::to_string(object.motionBlur.y) + " " + std::to_string(object.motionBlur.z));
    } else {
        object.motionBlur = {0, 0, 0};
        object.hasMotionBlur = false;
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
    
    // Parse MotionBlur if present
    if (objectData.contains("MotionBlur") && !objectData["MotionBlur"].is_null()) {
        object.motionBlur = parseTriplet<VectorFloatTriplet>(objectData["MotionBlur"]);
        object.hasMotionBlur = true;
        verbose("[+] MeshInstance MotionBlur parsed: " + std::to_string(object.motionBlur.x) + " " + std::to_string(object.motionBlur.y) + " " + std::to_string(object.motionBlur.z));
    } else {
        object.motionBlur = {0, 0, 0};
        object.hasMotionBlur = false;
    }
}

scene::Material* scene::Scene::getMaterialById(unsigned int id) {
    auto it = materialIdToIndex.find(id);
    if (it != materialIdToIndex.end()) {
        return &materials[it->second];
    }
    return nullptr;
}

const scene::BRDF* scene::Scene::getBRDFById(unsigned int id) const {
    auto it = brdfIdToIndex.find(id);
    if (it != brdfIdToIndex.end()) {
        return &brdfs[it->second];
    }
    return nullptr;
}

const scene::Image* scene::Scene::getImageById(unsigned int id) const {
    auto it = imageIdToIndex.find(id);
    if (it != imageIdToIndex.end()) {
        return &images[it->second];
    }
    return nullptr;
}

const scene::TextureMap* scene::Scene::getTextureMapById(unsigned int id) const {
    auto it = textureMapIdToIndex.find(id);
    if (it != textureMapIdToIndex.end()) {
        return &textureMaps[it->second];
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
    verbose("Point Lights: " + std::to_string(this->pointLights.size()) + "| lights: ");
    for (auto light : this->pointLights) {
        verbose("\t Light: " + std::to_string(light._id) + "| Position: " + std::to_string(light.position.x) + " " + std::to_string(light.position.y) + " " + std::to_string(light.position.z) + "| Intensity: " + std::to_string(light.intensity.x) + " " + std::to_string(light.intensity.y) + " " + std::to_string(light.intensity.z));
    }
    verbose("Area Lights: " + std::to_string(this->areaLights.size()) + "| lights: ");
    for (auto light : this->areaLights) {
        verbose("\t Light: " + std::to_string(light._id) + "| Position: " + std::to_string(light.position.x) + " " + std::to_string(light.position.y) + " " + std::to_string(light.position.z) + "| Normal: " + std::to_string(light.normal.x) + " " + std::to_string(light.normal.y) + " " + std::to_string(light.normal.z) + "| Size: " + std::to_string(light.size) + "| Radiance: " + std::to_string(light.radiance.x) + " " + std::to_string(light.radiance.y) + " " + std::to_string(light.radiance.z));
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
    std::string outFilename = filename;
    size_t dotPos = outFilename.find_last_of('.');
    if (dotPos != std::string::npos) {
        outFilename = outFilename.substr(0, dotPos) + ".png";
    } else {
        outFilename += ".png";
    }
    if (!stbi_write_png(outFilename.c_str(), width, height, 3, image, width * 3)) {
        throw std::runtime_error("Error: Failed to write PNG file: " + outFilename);
    }
}

scene::Scene::FaceParseResult scene::Scene::parsePLYFile(const std::string& plyFile, 
                                                          std::vector<scene::VectorFloatTriplet>& vertexList,
                                                          std::vector<scene::VectorFloatPair>& texCoordList) {
    FaceParseResult result;
    
    if (!std::ifstream(plyFile).good()) {
        throw std::runtime_error("Error: PLY file does not exist: " + plyFile);
    }
    happly::PLYData plyData(plyFile);

    // Remember how many vertices and texture coordinates we already have so we can offset indices
    // from the PLY file to point into the combined lists.
    size_t baseVertexIndex = vertexList.size();
    size_t baseTexCoordIndex = texCoordList.size();

    // Try to extract texture coordinates from PLY file
    // PLY files can have texture coordinates as vertex properties with various names:
    // "s"/"t", "u"/"v", "texture_u"/"texture_v", "texcoord_u"/"texcoord_v"
    bool hasTexCoords = false;
    std::vector<double> texU, texV;
    
    try {
        happly::Element& vertexElement = plyData.getElement("vertex");
        std::vector<std::string> propertyNames = vertexElement.getPropertyNames();
        
        // Try common texture coordinate property names
        std::vector<std::pair<std::string, std::string>> texCoordNames = {
            {"s", "t"},
            {"u", "v"},
            {"texture_u", "texture_v"},
            {"texcoord_u", "texcoord_v"},
            {"tx", "ty"}
        };
        
        for (const auto& namePair : texCoordNames) {
            bool hasU = false, hasV = false;
            for (const std::string& propName : propertyNames) {
                if (propName == namePair.first) hasU = true;
                if (propName == namePair.second) hasV = true;
            }
            
            if (hasU && hasV) {
                try {
                    texU = vertexElement.getProperty<double>(namePair.first);
                    texV = vertexElement.getProperty<double>(namePair.second);
                    hasTexCoords = true;
                    verbose("[+] Found texture coordinates in PLY file: " + namePair.first + "/" + namePair.second);
                    break;
                } catch (const std::exception&) {
                    // Try as float if double fails
                    try {
                        std::vector<float> uFloat = vertexElement.getProperty<float>(namePair.first);
                        std::vector<float> vFloat = vertexElement.getProperty<float>(namePair.second);
                        texU.resize(uFloat.size());
                        texV.resize(vFloat.size());
                        for (size_t i = 0; i < uFloat.size(); i++) {
                            texU[i] = static_cast<double>(uFloat[i]);
                            texV[i] = static_cast<double>(vFloat[i]);
                        }
                        hasTexCoords = true;
                        verbose("[+] Found texture coordinates in PLY file (as float): " + namePair.first + "/" + namePair.second);
                        break;
                    } catch (const std::exception&) {
                        // Continue to next name pair
                    }
                }
            }
        }
    } catch (const std::exception& e) {
        // No texture coordinates found, that's okay
        verbose("[!] No texture coordinates found in PLY file: " + std::string(e.what()));
    }

    // Extract texture coordinates if found
    if (hasTexCoords && texU.size() == texV.size()) {
        for (size_t i = 0; i < texU.size(); i++) {
            scene::VectorFloatPair uv;
            uv.x = texU[i];
            uv.y = texV[i];
            texCoordList.push_back(uv);
        }
        verbose("[+] Extracted " + std::to_string(texU.size()) + " texture coordinates from PLY file");
    }

    // Get faces from the PLY file. The underlying library (`happly`) already
    // handles both the standard `vertex_indices` and the common variant
    // `vertex_index` property names for face definitions.
    std::vector<std::vector<unsigned long>> faces = plyData.getFaceIndices();

    // Triangulate faces in case some polygons have more than 3 vertices
    // (e.g., quads in cube meshes). We use a simple fan triangulation:
    //   (v0, v1, v2), (v0, v2, v3), ...
    result.vertexFaces.reserve(faces.size()); // lower bound; may grow if we split quads/ngons
    result.texCoordFaces.reserve(faces.size());
    
    for (const auto& face : faces) {
        if (face.size() < 3) {
            // Degenerate face; skip it.
            continue;
        }

        for (size_t k = 1; k + 1 < face.size(); ++k) {
            scene::VectorIntTriplet vertexTri;
            // PLY indices are 0-based; shift them by baseVertexIndex so they refer to
            // the vertices we append below.
            vertexTri.x = static_cast<int>(baseVertexIndex + face[0]);
            vertexTri.y = static_cast<int>(baseVertexIndex + face[k]);
            vertexTri.z = static_cast<int>(baseVertexIndex + face[k + 1]);
            result.vertexFaces.push_back(vertexTri);
            
            // Texture coordinate indices: if we found texture coordinates, use them
            // Otherwise, use the same indices as vertices (fallback)
            if (hasTexCoords) {
                scene::VectorIntTriplet texTri;
                texTri.x = static_cast<int>(baseTexCoordIndex + face[0]);
                texTri.y = static_cast<int>(baseTexCoordIndex + face[k]);
                texTri.z = static_cast<int>(baseTexCoordIndex + face[k + 1]);
                result.texCoordFaces.push_back(texTri);
            } else {
                // Fallback: use vertex indices for texture coordinates
                result.texCoordFaces.push_back(vertexTri);
            }
        }
    }

    // Extract and add vertices
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
    
    return result;
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
    // First check if it's a direct mesh
    for (const auto& mesh : meshes) {
        if (mesh._id == id) {
            return &mesh;
        }
    }
    // If it's an instance, recursively follow the chain to find the actual base mesh
    for (const auto& instance : meshInstances) {
        if (instance._id == id) {
            // Recursively follow baseMeshId to get the actual mesh
            return findMeshOrInstanceById(instance.baseMeshId);
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
        
        verbose("[+] Mesh instance " + std::to_string(instance._id) + " references base " + std::to_string(instance.baseMeshId) + 
                " -> baseMesh=" + (instance.baseMesh ? "OK" : "NULL") + 
                " baseMeshIndex=" + std::to_string(instance.baseMeshIndex));
        
        if (instance.baseMesh) {
            // Inherit material from base mesh if not specified
            if (instance.material == nullptr) {
                instance.material = instance.baseMesh->material;
                verbose("[+] Mesh instance " + std::to_string(instance._id) + " inherited material from base mesh");
            } else {
                verbose("[+] Mesh instance " + std::to_string(instance._id) + " has explicit material " + std::to_string(instance.material->_id));
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
            // Derive a zoom factor from any Scaling transformations referenced by this camera.
            // This is required for camera-zoom scenes where the zoom is encoded as a scale
            // (e.g. via s3) rather than by directly animating the near plane / FOV.
            double zoomFactor = 1.0;
            for (const auto& ref : camera.transformations) {
                if (ref.type == 's') {
                    auto it = scalingIdToIndex.find(ref.id);
                    if (it != scalingIdToIndex.end()) {
                        const Scaling& s = scalings[it->second];
                        // Use the Y scale component as the zoom driver (matches provided scenes).
                        if (s.data.y > 0.0) {
                            zoomFactor *= s.data.y;
                        }
                    }
                }
            }

            Matrix4x4 cameraTransform = buildObjectTransformMatrix(*this, camera.transformations);
            camera.position = transformPoint(cameraTransform, camera.position);
            camera.gaze = normalize(transformDirection(cameraTransform, camera.gaze));
            camera.up = normalize(transformDirection(cameraTransform, camera.up));

            // Apply zoom by shrinking / expanding the camera's near-plane extents.
            // Larger zoomFactor => smaller near-plane window => narrower FOV => zoom in.
            if (zoomFactor != 1.0) {
                camera.nearPlane.x /= zoomFactor;
                camera.nearPlane.y /= zoomFactor;
                camera.nearPlane.z /= zoomFactor;
                camera.nearPlane.w /= zoomFactor;
                verbose("[+] Applied camera zoom factor " + std::to_string(zoomFactor) +
                        " to near plane for camera " + std::to_string(camera._id));
            }

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
    
    for (auto& areaLight : areaLights) {
        if (!areaLight.transformations.empty()) {
            Matrix4x4 lightTransform = buildObjectTransformMatrix(*this, areaLight.transformations);
            areaLight.position = transformPoint(lightTransform, areaLight.position);
            areaLight.normal = normalize(transformDirection(lightTransform, areaLight.normal));
            verbose("[+] Applied transformation to area light " + std::to_string(areaLight._id));
        }
    }
    
    for (auto& light : directionalLights) {
        if (!light.transformations.empty()) {
            Matrix4x4 lightTransform = buildObjectTransformMatrix(*this, light.transformations);
            light.direction = normalize(transformDirection(lightTransform, light.direction));
            verbose("[+] Applied transformation to directional light " + std::to_string(light._id));
        }
    }
    
    for (auto& light : spotLights) {
        if (!light.transformations.empty()) {
            Matrix4x4 lightTransform = buildObjectTransformMatrix(*this, light.transformations);
            light.position = transformPoint(lightTransform, light.position);
            light.direction = normalize(transformDirection(lightTransform, light.direction));
            verbose("[+] Applied transformation to spot light " + std::to_string(light._id));
        }
    }
    
    // SphericalDirectionalLight doesn't need transformations (uses image lookup)
    
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

scene::Image::~Image() {
    if (data) {
        stbi_image_free(data);
        data = nullptr;
    }
    if (hdrData) {
        // For EXR, use free (tinyexr uses malloc)
        // For HDR loaded with stbi_loadf, use stbi_image_free
        std::string lowerFilename = filename;
        std::transform(lowerFilename.begin(), lowerFilename.end(), lowerFilename.begin(), ::tolower);
        if (lowerFilename.length() >= 4 && lowerFilename.substr(lowerFilename.length() - 4) == ".exr") {
            free(hdrData);
        } else {
            stbi_image_free(hdrData);
        }
        hdrData = nullptr;
    }
}

scene::Image scene::parseImage(const json& imageData) {
    Image newImage;
    try {
        newImage._id = parseSingleValue<unsigned int>(imageData["_id"]);
        newImage.filename = imageData["_data"].get<std::string>();
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to parse Image: " + std::string(e.what()));
    }
    return newImage;
}

scene::TextureMap scene::parseTextureMap(const json& textureMapData) {
    TextureMap newTextureMap;
    auto normalizeKey = [](std::string s) {
        // trim
        const char* ws = " \t\n\r";
        size_t start = s.find_first_not_of(ws);
        size_t end = s.find_last_not_of(ws);
        if (start == std::string::npos) return std::string();
        s = s.substr(start, end - start + 1);
        // lowercase
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return (char)std::tolower(c); });
        return s;
    };
    auto parseUInt = [&](const json& v) -> unsigned int {
        if (v.is_string()) return parseSingleValue<unsigned int>(v.get<std::string>());
        if (v.is_number_unsigned()) return v.get<unsigned int>();
        if (v.is_number_integer()) return (unsigned int)std::max<long long>(0, v.get<long long>());
        throw std::runtime_error("Expected unsigned int (string/number)");
    };
    auto parseDouble = [&](const json& v) -> double {
        if (v.is_string()) return parseSingleValue<double>(v.get<std::string>());
        if (v.is_number_float()) return v.get<double>();
        if (v.is_number_integer()) return (double)v.get<long long>();
        if (v.is_number_unsigned()) return (double)v.get<unsigned long long>();
        throw std::runtime_error("Expected double (string/number)");
    };
    try {
        newTextureMap._id = parseUInt(textureMapData["_id"]);
        
        if (textureMapData.contains("_type") && !textureMapData["_type"].is_null()) {
            newTextureMap.type = normalizeKey(textureMapData["_type"].get<std::string>());
        } else {
            throw std::runtime_error("TextureMap missing _type");
        }
        
        if (newTextureMap.type == "image") {
            if (textureMapData.contains("ImageId") && !textureMapData["ImageId"].is_null()) {
                newTextureMap.imageId = parseUInt(textureMapData["ImageId"]);
            } else {
                throw std::runtime_error("Image texture missing ImageId");
            }
        }
        
        if (textureMapData.contains("DecalMode") && !textureMapData["DecalMode"].is_null()) {
            std::string decalModeStr = normalizeKey(textureMapData["DecalMode"].get<std::string>());
            if (decalModeStr == "replace_kd") {
                newTextureMap.decalMode = DecalMode::ReplaceKd;
            } else if (decalModeStr == "blend_kd") {
                newTextureMap.decalMode = DecalMode::BlendKd;
            } else if (decalModeStr == "replace_ks") {
                newTextureMap.decalMode = DecalMode::ReplaceKs;
            } else if (decalModeStr == "replace_background") {
                newTextureMap.decalMode = DecalMode::ReplaceBackground;
            } else if (decalModeStr == "replace_normal") {
                newTextureMap.decalMode = DecalMode::ReplaceNormal;
            } else if (decalModeStr == "bump_normal") {
                newTextureMap.decalMode = DecalMode::BumpNormal;
            } else if (decalModeStr == "replace_all") {
                newTextureMap.decalMode = DecalMode::ReplaceAll;
            }
        }
        
        if (textureMapData.contains("Interpolation") && !textureMapData["Interpolation"].is_null()) {
            std::string interpStr = normalizeKey(textureMapData["Interpolation"].get<std::string>());
            if (interpStr == "nearest") {
                newTextureMap.interpolation = InterpolationMode::Nearest;
            } else if (interpStr == "bilinear") {
                newTextureMap.interpolation = InterpolationMode::Bilinear;
            } else if (interpStr == "trilinear") {
                newTextureMap.interpolation = InterpolationMode::Trilinear;
            }
        }
        
        if (textureMapData.contains("BumpFactor") && !textureMapData["BumpFactor"].is_null()) {
            newTextureMap.bumpFactor = parseDouble(textureMapData["BumpFactor"]);
        }
        
        if (textureMapData.contains("NoiseScale") && !textureMapData["NoiseScale"].is_null()) {
            newTextureMap.noiseScale = parseDouble(textureMapData["NoiseScale"]);
        }
        
        if (textureMapData.contains("NoiseConversion") && !textureMapData["NoiseConversion"].is_null()) {
            std::string noiseConvStr = textureMapData["NoiseConversion"].get<std::string>();
            if (noiseConvStr == "absval") {
                newTextureMap.noiseConversion = NoiseConversion::AbsVal;
            } else if (noiseConvStr == "linear") {
                newTextureMap.noiseConversion = NoiseConversion::Linear;
            }
        }
        
        if (textureMapData.contains("NumOctaves") && !textureMapData["NumOctaves"].is_null()) {
            if (textureMapData["NumOctaves"].is_string()) {
                newTextureMap.numOctaves = parseSingleValue<int>(textureMapData["NumOctaves"].get<std::string>());
            } else if (textureMapData["NumOctaves"].is_number_integer()) {
                newTextureMap.numOctaves = textureMapData["NumOctaves"].get<int>();
            } else if (textureMapData["NumOctaves"].is_number_unsigned()) {
                newTextureMap.numOctaves = (int)textureMapData["NumOctaves"].get<unsigned int>();
            }
        }
        
        if (textureMapData.contains("Normalizer") && !textureMapData["Normalizer"].is_null()) {
            newTextureMap.normalizer = parseDouble(textureMapData["Normalizer"]);
            verbose("[+] TextureMap Normalizer parsed: " + std::to_string(newTextureMap.normalizer));
        }
        
        if (textureMapData.contains("Scale") && !textureMapData["Scale"].is_null()) {
            newTextureMap.scale = parseDouble(textureMapData["Scale"]);
        }
        
        if (textureMapData.contains("Offset") && !textureMapData["Offset"].is_null()) {
            newTextureMap.offset = parseTriplet<VectorFloatTriplet>(textureMapData["Offset"]);
        }
        
        if (textureMapData.contains("BlackColor") && !textureMapData["BlackColor"].is_null()) {
            newTextureMap.blackColor = parseTriplet<VectorFloatTriplet>(textureMapData["BlackColor"]);
        }
        
        if (textureMapData.contains("WhiteColor") && !textureMapData["WhiteColor"].is_null()) {
            newTextureMap.whiteColor = parseTriplet<VectorFloatTriplet>(textureMapData["WhiteColor"]);
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to parse TextureMap: " + std::string(e.what()));
    }
    return newTextureMap;
}

std::vector<scene::VectorFloatPair> scene::parseTexCoordData(const json& texCoordData) {
    std::vector<VectorFloatPair> texCoords;
    try {
        if (texCoordData.contains("_type") && !texCoordData["_type"].is_null()) {
            std::string type = texCoordData["_type"].get<std::string>();
            if (type != "uv") {
                throw std::runtime_error("Unsupported TexCoordData type: " + type);
            }
        }
        
        if (texCoordData.contains("_data") && !texCoordData["_data"].is_null()) {
            std::stringstream stream(texCoordData["_data"].get<std::string>());
            VectorFloatPair uv;
            while (stream >> uv.x >> uv.y) {
                texCoords.push_back(uv);
            }
        } else {
            throw std::runtime_error("TexCoordData missing _data");
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to parse TexCoordData: " + std::string(e.what()));
    }
    return texCoords;
}