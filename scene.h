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

   struct VectorFloatPair {
        double x, y;
   };

   struct VectorFloatTriplet {
        double x, y, z;
   };

   // 5D float vector for precomputed sample dimensions (jitter, time, lens, etc.)
   struct VectorFloatPenta {
        double x, y, z, w, v;
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
    struct TonemapSettings {
        std::string tmo;  // "Photographic", "Filmic", "ACES"
        std::string tmoOptions;  // "key burnOutPercent" as string
        double saturation;
        double gamma;
        std::string extension;  // Extension for output filename
    };

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
        int numSamples = 1;  // Perfect square (1, 4, 16, etc.) for jittered sampling
        double apertureSize = 0.0;  // 0 means no depth-of-field
        double focusDistance = 0.0;  // Focus distance for depth-of-field
        std::vector<TonemapSettings> tonemapSettings;  // Tone mapping settings (can be multiple)
        // Per-pixel per-sample precomputed random dimensions:
        //   x, y : subpixel jitter in [0,1)
        //   z    : motion blur time in [0,1)
        //   w, v : extra random dims reused for lens/roughness/area lights
        VectorFloatPenta* samples = nullptr;
        
        // Path tracing settings
        std::string renderer = "";  // "PathTracing" or "" (default shading)
        bool importanceSampling = false;
        bool nextEventEstimation = false;
        std::string misHeuristic = "";  // "balance", "power", "01"
        bool russianRoulette = false;
        int maxRecursionDepth = 6;
        int minRecursionDepth = 0;
        int splittingFactor = 1;
        double sampleMaxVal = 0.0;  // 0 = no clamping
    };

    struct PointLight {
        unsigned int _id;
        VectorFloatTriplet position;
        VectorFloatTriplet intensity;
        std::vector<TransformationRef> transformations;
    };

    struct AreaLight {
        unsigned int _id;
        VectorFloatTriplet position;  // Center point of the area light
        VectorFloatTriplet normal;    // Surface normal of the light
        double size;                  // Edge length of the square area light
        VectorFloatTriplet radiance;   // Radiance of the light
        std::vector<TransformationRef> transformations;
    };

    struct DirectionalLight {
        unsigned int _id;
        VectorFloatTriplet direction;
        VectorFloatTriplet radiance;
        std::vector<TransformationRef> transformations;
    };

    struct SpotLight {
        unsigned int _id;
        VectorFloatTriplet position;
        VectorFloatTriplet direction;
        VectorFloatTriplet intensity;
        double coverageAngle;  // in degrees
        double falloffAngle;  // in degrees
        std::vector<TransformationRef> transformations;
    };

    struct SphericalDirectionalLight {
        unsigned int _id;
        unsigned int imageId;
        std::string type;  // "latlong" or "probe"
        std::string sampler;  // "uniform" or "cosine"
        std::vector<TransformationRef> transformations;
    };

    // Tbh no need for defining ambient light as a struct, but I did it for consistency.
    struct AmbientLight {
        VectorFloatTriplet intensity;
    };

    enum class BRDFType {
        OriginalBlinnPhong,
        OriginalPhong,
        ModifiedBlinnPhong,
        ModifiedPhong,
        TorranceSparrow
    };

    struct BRDF {
        unsigned int _id;
        BRDFType type;
        double exponent = 1.0;
        bool normalized = false;
        bool kdFresnel = false;  // For TorranceSparrow: use (1-F)*kd/pi instead of kd/pi
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
        double roughness = 0.0;  // Roughness for mirrors, conductors, and dielectrics
        unsigned int brdfId = 0;  // 0 = use default (OriginalBlinnPhong)
    };

    struct Image {
        unsigned int _id;
        std::string filename;
        unsigned char* data = nullptr;  // For LDR images
        float* hdrData = nullptr;  // For HDR images (EXR/HDR)
        int width = 0;
        int height = 0;
        int channels = 0;
        bool isHDR = false;  // True if loaded as HDR
        ~Image();
    };

    enum class InterpolationMode {
        Nearest,
        Bilinear,
        Trilinear
    };

    enum class DecalMode {
        ReplaceKd,
        BlendKd,
        ReplaceKs,
        ReplaceBackground,
        ReplaceNormal,
        BumpNormal,
        ReplaceAll
    };

    enum class NoiseConversion {
        AbsVal,
        Linear
    };

    struct TextureMap {
        unsigned int _id;
        std::string type;  // "image", "perlin", "checkerboard"
        unsigned int imageId = 0;  // Only for image type
        DecalMode decalMode = DecalMode::ReplaceKd;
        InterpolationMode interpolation = InterpolationMode::Nearest;
        double bumpFactor = 0.01;  // Reasonable default (homework scenes use 0.01)
        double noiseScale = 1.0;
        NoiseConversion noiseConversion = NoiseConversion::Linear;
        int numOctaves = 1;
        double normalizer = 0.0;  // 0 = not set (use 255), otherwise divide pixel by this value
        
        // Checkerboard parameters
        double scale = 1.0;
        VectorFloatTriplet offset = {0, 0, 0};
        VectorFloatTriplet blackColor = {0, 0, 0};
        VectorFloatTriplet whiteColor = {1, 1, 1};
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
        VectorFloatTriplet motionBlur = {0, 0, 0};  // Displacement vector for motion blur (translational only)
        bool hasMotionBlur = false;
        std::vector<unsigned int> textureIds;  // TextureMap IDs
    };

    struct Mesh : public Object {
        char shadingMode;
        std::vector<VectorIntTriplet> faces;  // Vertex indices
        std::vector<VectorIntTriplet> texCoordIndices;  // Texture coordinate indices (separate from vertex indices)
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

    struct LightSphere : public Object {
        unsigned int center;
        double radius;
        VectorFloatTriplet radiance;
    };

    struct LightMesh : public Object {
        char shadingMode = 'f';
        std::vector<VectorIntTriplet> faces;  // Vertex indices
        std::vector<VectorIntTriplet> texCoordIndices;  // Texture coordinate indices
        VectorFloatTriplet radiance;
        double totalArea = 0.0;  // Precomputed for sampling
        std::vector<double> cdfTriangleAreas;  // For importance sampling (cumulative distribution function)
    };

    struct Scene {
        VectorFloatTriplet backgroundColor;
        double shadowRayEpsilon;
        double intersectionTestEpsilon;
        int maxRecursionDepth = 5;
        std::vector<Camera> cameras;
        AmbientLight ambientLight;
        std::vector<PointLight> pointLights;
        std::vector<AreaLight> areaLights;
        std::vector<DirectionalLight> directionalLights;
        std::vector<SpotLight> spotLights;
        std::vector<SphericalDirectionalLight> sphericalDirectionalLights;
        std::vector<Material> materials;
        std::map<unsigned int, size_t> materialIdToIndex;
        std::vector<BRDF> brdfs;
        std::map<unsigned int, size_t> brdfIdToIndex;
        std::vector<VectorFloatTriplet> vertices;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;
        std::vector<Plane> planes;
        std::vector<MeshInstance> meshInstances;
        std::vector<LightSphere> lightSpheres;
        std::vector<LightMesh> lightMeshes;
        
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

        // Texture storage
        std::vector<Image> images;
        std::map<unsigned int, size_t> imageIdToIndex;
        std::vector<TextureMap> textureMaps;
        std::map<unsigned int, size_t> textureMapIdToIndex;
        std::vector<VectorFloatPair> texCoords;  // Texture coordinates (UV pairs)
        unsigned int backgroundTextureId = 0;  // TextureMap id used for background (ReplaceBackground)

        // Informational variables
        unsigned char currentCameraIndex;
        std::string baseDirectory;
        
        bool enableBackFaceCulling = true;
        void loadSceneFromFile(const std::string& filename);
        
        Material* getMaterialById(unsigned int id);
        const BRDF* getBRDFById(unsigned int id) const;
        const Image* getImageById(unsigned int id) const;
        const TextureMap* getTextureMapById(unsigned int id) const;
        
        template<typename T> 
        std::vector<T> parseObjects(const json& objectsData);

        template<typename T>
        void parseSpecificAttributes(T& object, const json& objectData);

        struct FaceParseResult {
            std::vector<VectorIntTriplet> vertexFaces;
            std::vector<VectorIntTriplet> texCoordFaces;
        };
        
        FaceParseResult parseFacesWithOffsets(const json& facesData);
        std::vector<VectorIntTriplet> parseFaces(const json& facesData);
        
        FaceParseResult parsePLYFile(const std::string& plyFile, 
                                      std::vector<VectorFloatTriplet>& vertexList,
                                      std::vector<VectorFloatPair>& texCoordList);
        
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
    AreaLight parseAreaLight(const json& areaLightData);
    DirectionalLight parseDirectionalLight(const json& directionalLightData);
    SpotLight parseSpotLight(const json& spotLightData);
    SphericalDirectionalLight parseSphericalDirectionalLight(const json& sphericalDirectionalLightData);
    Material parseMaterial(const json& materialData);
    BRDF parseBRDF(const json& brdfData, BRDFType type);
    std::vector<VectorFloatTriplet> parseVertex(const json& vertexData);
    
    Object parseObject(const json& objectData);
    
    // Transformation parsing functions
    Translation parseTranslation(const json& translationData);
    Scaling parseScaling(const json& scalingData);
    Rotation parseRotation(const json& rotationData);
    Composite parseComposite(const json& compositeData);
    
    // Texture parsing functions
    Image parseImage(const json& imageData);
    TextureMap parseTextureMap(const json& textureMapData);
    std::vector<VectorFloatPair> parseTexCoordData(const json& texCoordData);

    struct Ray {
        VectorFloatTriplet origin;
        VectorFloatTriplet direction;
        int depth;
        bool shadowRay;
        bool reflectionRay;
        bool refractionRay;
        double time;      // Time value for motion blur (0.0 to 1.0)
        double random1;   // Precomputed random in [0,1] (lens/roughness/area light)
        double random2;   // Precomputed random in [0,1] (lens/roughness/area light)
        
        // Default constructor
        Ray()
            : origin({0, 0, 0}),
              direction({0, 0, 0}),
              depth(0),
              shadowRay(false),
              reflectionRay(false),
              refractionRay(false),
              time(0.0),
              random1(0.0),
              random2(0.0) {}
        
        // Constructor with all parameters
        Ray(const VectorFloatTriplet& o,
            const VectorFloatTriplet& d,
            int dep, 
            bool shadow,
            bool reflection,
            bool refraction,
            double t,
            double r1 = 0.0,
            double r2 = 0.0)
            : origin(o),
              direction(d),
              depth(dep),
              shadowRay(shadow), 
              reflectionRay(reflection),
              refractionRay(refraction),
              time(t),
              random1(r1),
              random2(r2) {}
    };

    struct Intersection {
        bool hit = false;
        double distance = 0.0;
        VectorFloatTriplet point;
        VectorFloatTriplet geometricNormal;
        VectorFloatTriplet shadingNormal;
        

        enum class Kind { None, Plane, Sphere, Triangle, Mesh, AreaLight, LightSphere, LightMesh } kind = Kind::None;
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
        bool iterativeSampling;
        bool useGUI;
        Args()
            : isMultiThreaded(true),
              useBVH(false),
              enableBackFaceCulling(true),
              iterativeSampling(false),
              useGUI(false) {}
    };
}

#endif