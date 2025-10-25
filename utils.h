#ifndef __UTILS__
#define __UTILS__

#include "scene.h"

using namespace std;
using namespace scene;


/* Precompute Functions */
void precomputeMeshNormals(const vector<Mesh>& meshes, vector<vector<VectorFloatTriplet>>& meshNormals, const vector<VectorFloatTriplet>& vertices);
void precomputeTriangleNormals(const vector<Triangle>& triangles, vector<VectorFloatTriplet>& triangleNormals, const vector<VectorFloatTriplet>& vertices);



/* Write Functions */
void writePPM(const string& filename, unsigned char* image, int width, int height);


/* Ray Functions */
Ray castRay(const Camera& camera, int x, int y, int width, int height);

Intersection intersect(const Scene& scene, const Ray& ray, const vector<vector<VectorFloatTriplet>>& meshNormals, const vector<VectorFloatTriplet>& triangleNormals);

VectorFloatTriplet computePixelColor(const Scene& scene, const Ray& ray, const Intersection& intersection, const vector<vector<VectorFloatTriplet>>& meshNormals, const vector<VectorFloatTriplet>& triangleNormals);

VectorFloatTriplet computeShading(const Scene& scene, const Ray& ray, const Intersection& intersection, const vector<vector<VectorFloatTriplet>>& meshNormals, const vector<VectorFloatTriplet>& triangleNormals);

/* Operator Overloads */
/* VectorFloatTriplet */
VectorFloatTriplet crossProduct(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
float dotProduct(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet normalize(const VectorFloatTriplet& a);
VectorFloatTriplet operator-(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet operator+(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet operator*(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet& operator+=(VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet& operator-=(VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet& operator*=(VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet& operator/=(VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet operator-(const VectorFloatTriplet& a);
VectorFloatTriplet operator*(const VectorFloatTriplet& a, const float& b);
VectorFloatTriplet operator*(const float& a, const VectorFloatTriplet& b);

/* VectorIntTriplet */
VectorIntTriplet crossProduct(const VectorIntTriplet& a, const VectorIntTriplet& b);
int dotProduct(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet normalize(const VectorIntTriplet& a);
VectorIntTriplet operator-(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet operator+(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet operator*(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet& operator+=(VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet& operator-=(VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet& operator*=(VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet& operator/=(VectorIntTriplet& a, const VectorIntTriplet& b);


#endif