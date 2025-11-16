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
        double x, y, z;
    };

    struct VectorIntTriplet {
        int x, y, z;
    };

    struct VectorFloatQuad {
        double x, y, z, w;
    };

    struct VectorIntQuad {
        int x, y, z, w;
    };

    /*
    * Transformation structs
    */
    struct Translation {
        unsigned int _id;
        VectorFloatTriplet data; // dx, dy, dz
    };

    struct Scaling {
        unsigned int _id;
        VectorFloatTriplet data; // sx, sy, sz
    };

    struct Rotation {
        unsigned int _id;
        double angle; // in degrees
        VectorFloatTriplet axis; // x, y, z
    };

    struct Composite {
        unsigned int _id;
        double data[16]; // 4x4 matrix in row-major order
    };

    struct TransformationRef {
        char type; // 't' = translation, 's' = scaling, 'r' = rotation, 'c' = composite
        unsigned int id;
    };

    struct Matrix4x4 {
        double m[16];
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
        double nearDistance;
        VectorIntPair imageResolution;
        std::string imageName;
        std::vector<TransformationRef> transformations;
    };

    struct PointLight {
        unsigned int _id;
        VectorFloatTriplet position;
        VectorFloatTriplet intensity;
        std::vector<TransformationRef> transformations;
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
        double phongExponent;
        bool isMirror = false;
        VectorFloatTriplet mirrorReflectance;
        std::string type = "";
        double refractionIndex = 1.0;
        double absorptionIndex = 0.0;
        VectorFloatTriplet absorptionCoefficient = {0, 0, 0};
    };
    
    // Forward declare AABB
    struct AABB;
    
    struct Object {
        unsigned int _id;
        Material* material;
        std::vector<TransformationRef> transformations;
        Matrix4x4* transformMatrix = nullptr;
        Matrix4x4* inverseTransformMatrix = nullptr;
        Matrix4x4* normalMatrix = nullptr;
        bool hasTransformation = false;
        bool hasNegativeScale = false;  // True if transformation includes reflection (negative scale)
        AABB* worldSpaceBounds = nullptr;  // World-space bounding box for transformed objects
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
        double radius;
    };

    struct Plane : public Object {
        unsigned int point;
        VectorFloatTriplet normal;
    };

    struct MeshInstance : public Object {
        unsigned int baseMeshId;
        bool resetTransform = false;
        const Mesh* baseMesh = nullptr;
        int baseMeshIndex = -1;
    };

    struct Scene {
        VectorFloatTriplet backgroundColor;
        double shadowRayEpsilon;
        double intersectionTestEpsilon;
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
        std::vector<MeshInstance> meshInstances;
        
        // Transformation storage
        std::vector<Translation> translations;
        std::vector<Scaling> scalings;
        std::vector<Rotation> rotations;
        std::vector<Composite> composites;
        std::map<unsigned int, size_t> translationIdToIndex;
        std::map<unsigned int, size_t> scalingIdToIndex;
        std::map<unsigned int, size_t> rotationIdToIndex;
        std::map<unsigned int, size_t> compositeIdToIndex;

        // Additional precomputed data
        std::vector<VectorFloatTriplet> triangleNormals;
        std::vector<std::vector<VectorFloatTriplet>> meshVertexNormals;
        std::vector<std::vector<double>> cameraTriangleDeterminant;
        std::vector<std::vector<std::vector<double>>> cameraMeshDeterminant;

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
        
        void precomputeTransformations();
        int findBaseMeshIndex(unsigned int meshOrInstanceId) const;
        const Mesh* findMeshOrInstanceById(unsigned int id) const;
    };

    // Helper function to parse transformation string into TransformationRef vector
    std::vector<TransformationRef> parseTransformationString(const std::string& transformStr);
    
    Camera parseCamera(const json& cameraData);
    PointLight parsePointLight(const json& pointLightData);
    Material parseMaterial(const json& materialData);
    std::vector<VectorFloatTriplet> parseVertex(const json& vertexData);
    
    std::vector<VectorIntTriplet> parsePLYFile(const std::string& plyFile, 
                                                std::vector<VectorFloatTriplet>& vertexList);
    
    Object parseObject(const json& objectData);
    
    // Transformation parsing functions
    Translation parseTranslation(const json& translationData);
    Scaling parseScaling(const json& scalingData);
    Rotation parseRotation(const json& rotationData);
    Composite parseComposite(const json& compositeData);

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
        double distance = 0.0;
        VectorFloatTriplet point;
        VectorFloatTriplet geometricNormal;
        VectorFloatTriplet shadingNormal;
        

        enum class Kind { None, Plane, Sphere, Triangle, Mesh } kind = Kind::None;
        int containerIndex = -1;
        int faceIndex = -1;
        
        double beta = 0.0;
        double gamma = 0.0;
        
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