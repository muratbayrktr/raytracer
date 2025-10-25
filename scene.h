#ifndef __SCENE__
#define __SCENE__

#include "scene.h"
#include "json.hpp"
#include <string>
#include <vector>
#include <map>

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

    struct Object {
        unsigned int _id;
        Material* material;
    };

    struct Mesh : public Object {
        std::vector<VectorIntTriplet> faces;
    };

    struct Triangle : public Object {
        VectorIntTriplet indices;
    };

    struct Sphere : public Object {
        unsigned int center;
        float radius;
    };

    struct Plane : public Object {
        unsigned int point;
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
        std::map<unsigned int, size_t> materialIdToIndex;
        std::vector<VectorFloatTriplet> vertices;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;
        std::vector<Plane> planes;

        void loadSceneFromFile(const std::string& filename);
        
        Material* getMaterialById(unsigned int id);
        
        template<typename T> 
        std::vector<T> parseObjects(const json& objectsData);

        template<typename T>
        void parseSpecificAttributes(T& object, const json& objectData);

        void getSummary();
    };

    Camera parseCamera(const json& cameraData);
    PointLight parsePointLight(const json& pointLightData);
    Material parseMaterial(const json& materialData);
    std::vector<VectorFloatTriplet> parseVertex(const json& vertexData);
    std::vector<VectorIntTriplet> parseFaces(const json& facesData);
    Object parseObject(const json& objectData);

    struct Ray {
        VectorFloatTriplet origin;
        VectorFloatTriplet direction;
        int depth;
    };

    struct Intersection {
        bool hit;
        float distance;
        VectorFloatTriplet point;
        VectorFloatTriplet normal;
        Material* material;
    };
}

#endif