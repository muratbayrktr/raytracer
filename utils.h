#ifndef __UTILS__
#define __UTILS__

#include "scene.h"

using namespace std;
using namespace scene;


/* Precompute Functions */
void precomputeMeshNormals(const vector<Mesh>& meshes, vector<VectorFloatTriplet>& meshNormals, const vector<VectorFloatTriplet>& vertices);
void precomputeTriangleNormals(const vector<Triangle>& triangles, vector<VectorFloatTriplet>& triangleNormals, const vector<VectorFloatTriplet>& vertices);



/* Write Functions */
void writePPM(const string& filename, const vector<unsigned char>& image, int width, int height);


/* Operator Overloads */
/* VectorFloatTriplet */
VectorFloatTriplet crossProduct(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
float dotProduct(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet normalize(const VectorFloatTriplet& a);
VectorFloatTriplet operator-(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet operator+(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet operator*(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet operator+=(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet operator-=(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet operator*=(const VectorFloatTriplet& a, const VectorFloatTriplet& b);
VectorFloatTriplet operator/=(const VectorFloatTriplet& a, const VectorFloatTriplet& b);

/* VectorIntTriplet */
VectorIntTriplet crossProduct(const VectorIntTriplet& a, const VectorIntTriplet& b);
int dotProduct(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet normalize(const VectorIntTriplet& a);
VectorIntTriplet operator-(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet operator+(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet operator*(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet operator+=(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet operator-=(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet operator*=(const VectorIntTriplet& a, const VectorIntTriplet& b);
VectorIntTriplet operator/=(const VectorIntTriplet& a, const VectorIntTriplet& b);


#endif