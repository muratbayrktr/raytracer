#ifndef __UTILS__
#define __UTILS__

#include "scene.h"

using namespace std;
using namespace scene;

Matrix4x4 buildTranslationMatrix(const Translation& t);
Matrix4x4 buildScalingMatrix(const Scaling& s);
Matrix4x4 buildRotationMatrix(const Rotation& r);
Matrix4x4 buildCompositeMatrix(const Composite& c);
Matrix4x4 identityMatrix();
Matrix4x4 multiplyMatrices(const Matrix4x4& a, const Matrix4x4& b);
Matrix4x4 invertMatrix(const Matrix4x4& m);
Matrix4x4 transposeMatrix(const Matrix4x4& m);

VectorFloatTriplet transformPoint(const Matrix4x4& m, const VectorFloatTriplet& p);
VectorFloatTriplet transformDirection(const Matrix4x4& m, const VectorFloatTriplet& d);
VectorFloatTriplet transformNormal(const Matrix4x4& invTranspose, const VectorFloatTriplet& n);

Matrix4x4 buildObjectTransformMatrix(const Scene& scene, const std::vector<TransformationRef>& transforms);
bool hasNegativeScale(const Matrix4x4& m);

void clamp(VectorFloatTriplet& color, int min, int max);

/* Ray Functions */
Ray castRay(const Camera& camera,
            double x,
            double y,
            int width,
            int height,
            double time = 0.0,
            double random1 = 0.0,
            double random2 = 0.0);

/* Intersection Functions */
Intersection intersect(const Scene& scene, Ray& ray); 

bool rayHitsPlane(Ray& ray, const Plane& plane, const vector<VectorFloatTriplet>& vertices, double& t_min, Intersection& intersection, int planeIndex, double minDistance = 0.0);
bool rayHitsSphere(Ray& ray, const Sphere& sphere, const vector<VectorFloatTriplet>& vertices, double& t_min, Intersection& intersection, int sphereIndex, double minDistance = 0.0);
bool rayHitsTriangle(Ray& ray, const VectorIntTriplet& face, const vector<VectorFloatTriplet>& vertices, double& t_min, Intersection& intersection, double intersectionTestEpsilon, double determinantT, Material* material, bool enableBackFaceCulling, int containerIndex, int faceIndex, double minDistance = 0.0);
bool rayHitsMesh(Ray& ray, const Mesh& mesh, const vector<VectorFloatTriplet>& vertices, const vector<double>& determinants, double& t_min, Intersection& intersection, double intersectionTestEpsilon, scene::MeshBVH* bvh, bool enableBackFaceCulling, int meshIndex, const Matrix4x4* transformMatrix = nullptr, const Matrix4x4* inverseTransformMatrix = nullptr, const Matrix4x4* normalMatrix = nullptr, const Scene* scene = nullptr, Material* materialOverride = nullptr, const scene::AABB* worldSpaceBoundsOverride = nullptr, double minDistance = 0.0);
bool rayHitsGaussianField(Ray& ray, const GaussianField& field, double& t_min, Intersection& intersection, int fieldIndex, double minDistance = 0.0);
void rayMarchGaussianVolume(Ray& ray, const GaussianField& field, double t_min, double minDistance, VectorFloatTriplet& accumulatedColor, double& accumulatedOpacity);

/* Pixel Color Functions */
VectorFloatTriplet computePixelColor(const Scene& scene, Ray& ray, const Intersection& intersection);

/* Shading Functions */
VectorFloatTriplet computeShading(const Scene& scene, Ray& ray, const Intersection& intersection);

/* Reflection Functions */
Ray reflect(Ray& ray, const VectorFloatTriplet normal, VectorFloatTriplet point, const double shadowRayEpsilon);

/* Fresnel Functions */
double fresnelConductor(double cosTheta, double n, double k);
double fresnelDielectric(double cosTheta, double n1, double n2);

/* Refraction Functions */
Ray refract(Ray& ray, const VectorFloatTriplet normal, double n1, double n2, VectorFloatTriplet point, const double shadowRayEpsilon, bool& totalInternalReflection);

/* Shadow */
bool isInShadow(const Scene& scene, Ray& ray, const PointLight& light, const Intersection& intersection);
double gaussianFieldTransmittance(const Scene& scene, const Ray& shadowRay, double maxDist);

/* Utility Functions */
void orthonormalBasis(const VectorFloatTriplet& n, VectorFloatTriplet& u, VectorFloatTriplet& v);

void printPerfStats();
void printPerfStatsInline();

#endif