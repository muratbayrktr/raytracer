#ifndef __SCENE__
#define __SCENE__

#include "scene.h"
#include "json.hpp"
#include <string>
#include <vector>

namespace scene {
    using json = nlohmann::json;

    /*
    * Atomic vector types
    * I am defining these structs to have a simple and faster access to the data.
    */
   struct VectorIntPair {
        int x, y;
   };

    struct VectorFloatTriplet {
        float x, y, z;
    };

    struct VectorIntTriplet {
        int x, y, z;
    };

    struct VectorFloatQuad {
        float x, y, z, w;
    };

    struct VectorIntQuad {
        int x, y, z, w;
    };

    struct VertexData {
        std::string _type;
        std::vector<float> data; // arbitrary number of floats, dataSize for fast access
        unsigned int dataSize;
    };

    /*
    * Scene objects
    * More sophisticated structs that define the scene objects.
    */
    struct Camera {
        unsigned int _id;
        VectorFloatTriplet position;
        VectorFloatTriplet gaze;
        VectorFloatTriplet up;
        VectorFloatQuad nearPlane;
        float nearDistance;
        VectorIntPair imageResolution;
        std::string imageName;
    };

    struct PointLight {
        unsigned int _id;
        VectorFloatTriplet position;
        VectorFloatTriplet intensity;
    };

    // Tbh no need for defining ambient light as a struct, but I did it for consistency.
    struct AmbientLight {
        VectorFloatTriplet intensity;
    };

    struct Material {
        unsigned int _id;
        VectorFloatTriplet ambientReflectance;
        VectorFloatTriplet diffuseReflectance;
        VectorFloatTriplet specularReflectance;
        float phongExponent;
    };

    struct Mesh {
        unsigned int _id;
        Material* material; // Define materials once, use pointer to access
        VertexData faces;
    };

    struct Triangle {
        unsigned int _id;
        Material* material;
        VectorFloatTriplet indices;
    };

    struct Sphere {
        unsigned int _id;
        Material* material;
        float center;
        float radius;
    };

    struct Plane {
        unsigned int _id;
        Material* material;
        float point;
        VectorFloatTriplet normal;
    };

    struct Scene {
        VectorFloatTriplet backgroundColor;
        float shadowRayEpsilon;
        float intersectionTestEpsilon;
        std::vector<Camera> cameras;
        AmbientLight ambientLight;
        std::vector<PointLight> pointLights;
        std::vector<Material> materials;
        VertexData vertexData;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;
        std::vector<Plane> planes;

        void loadSceneFromFile(const std::string& filename);
    };

    Camera parseCamera(const json& cameraData);
    PointLight parsePointLight(const json& pointLightData);
    Material parseMaterial(const json& materialData);
    
}

#endif