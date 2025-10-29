#ifndef __UTILS__
#define __UTILS__

#include "scene.h"

using namespace std;
using namespace scene;

void clamp(VectorFloatTriplet& color, int min, int max);

/* Ray Functions */
Ray castRay(const Camera& camera, int x, int y, int width, int height);

/* Intersection Functions */
Intersection intersect(const Scene& scene, Ray& ray); 

bool rayHitsPlane(Ray& ray, const Plane& plane, const vector<VectorFloatTriplet>& vertices, float& t_min, Intersection& intersection);
bool rayHitsSphere(Ray& ray, const Sphere& sphere, const vector<VectorFloatTriplet>& vertices, float& t_min, Intersection& intersection);
bool rayHitsTriangle(Ray& ray, const VectorIntTriplet& face, const VectorFloatTriplet& triangleNormal, const vector<VectorFloatTriplet>& vertices, float& t_min, Intersection& intersection, float intersectionTestEpsilon, float determinantT, Material* material);
bool rayHitsMesh(Ray& ray, const Mesh& mesh, const vector<VectorFloatTriplet>& normals, const vector<VectorFloatTriplet>& vertices, const vector<float>& determinants, float& t_min, Intersection& intersection, float intersectionTestEpsilon, scene::MeshBVH* bvh = nullptr);

/* Pixel Color Functions */
VectorFloatTriplet computePixelColor(const Scene& scene, Ray& ray, const Intersection& intersection);

/* Shading Functions */
VectorFloatTriplet computeShading(const Scene& scene, Ray& ray, const Intersection& intersection);

/* Reflection Functions */
Ray reflect(Ray& ray, const VectorFloatTriplet normal, VectorFloatTriplet point, const float shadowRayEpsilon);

/* Shadow */
bool isInShadow(const Scene& scene, Ray& ray, const PointLight& light, const Intersection& intersection);

#endif