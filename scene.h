#ifndef __SCENE__
#define __SCENE__

#include "scene.h"
#include "json.hpp"
#include <string>
#include <vector>
#include <map>
#include <vector>

namespace scene {
    using json = nlohmann::json;
    

    class MeshBVH;

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
        bool isMirror = false;
        VectorFloatTriplet mirrorReflectance;
        std::string type = "";
        float refractionIndex = 1.0f;
        float absorptionIndex = 0.0f;
        VectorFloatTriplet absorptionCoefficient = {0, 0, 0};
    };

    struct Object {
        unsigned int _id;
        Material* material;
    };

    struct Mesh : public Object {
        char shadingMode;
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
        int maxRecursionDepth = 0;
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

        // Additional precomputed data
        std::vector<VectorFloatTriplet> triangleNormals;
        std::vector<std::vector<VectorFloatTriplet>> meshVertexNormals;
        std::vector<std::vector<float>> cameraTriangleDeterminant;
        std::vector<std::vector<std::vector<float>>> cameraMeshDeterminant;

        std::vector<MeshBVH*> meshBVHs;

        // Informational variables
        unsigned char currentCameraIndex;
        std::string baseDirectory;
        
        bool enableBackFaceCulling = true;
        void loadSceneFromFile(const std::string& filename);
        
        Material* getMaterialById(unsigned int id);
        
        template<typename T> 
        std::vector<T> parseObjects(const json& objectsData);

        template<typename T>
        void parseSpecificAttributes(T& object, const json& objectData);

        std::vector<VectorIntTriplet> parseFaces(const json& facesData);
        
        void buildBVH();
        
        void getSummary();
        
        void writePPM(const std::string& filename, unsigned char* image, int width, int height);
    };

    Camera parseCamera(const json& cameraData);
    PointLight parsePointLight(const json& pointLightData);
    Material parseMaterial(const json& materialData);
    std::vector<VectorFloatTriplet> parseVertex(const json& vertexData);
    
    std::vector<VectorIntTriplet> parsePLYFile(const std::string& plyFile, 
                                                std::vector<VectorFloatTriplet>& vertexList);
    
    Object parseObject(const json& objectData);

    struct Ray {
        VectorFloatTriplet origin;
        VectorFloatTriplet direction;
        int depth;
        bool shadowRay;
        bool reflectionRay;
        bool refractionRay;
    };

    struct Intersection {
        bool hit = false;
        float distance = 0.0f;
        VectorFloatTriplet point;
        VectorFloatTriplet geometricNormal;
        VectorFloatTriplet shadingNormal;
        

        enum class Kind { None, Plane, Sphere, Triangle, Mesh } kind = Kind::None;
        int containerIndex = -1;
        int faceIndex = -1;
        
        float beta = 0.0f;
        float gamma = 0.0f;
        
        Material* material = nullptr;
    };

    struct Args {
        std::string sceneFile;
        bool isMultiThreaded;
        bool useBVH;
        bool enableBackFaceCulling;
        Args() : isMultiThreaded(true), useBVH(false), enableBackFaceCulling(true) {}
    };
}

#endif