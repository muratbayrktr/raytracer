#ifndef __PRECOMPUTE__
#define __PRECOMPUTE__

#include "scene.h"
#include "utils.h"

using namespace std;
using namespace scene;

/* Precompute Functions */
void precomputeMeshNormals(const vector<Mesh>& meshes, vector<vector<VectorFloatTriplet>>& meshVertexNormals, const vector<VectorFloatTriplet>& vertices);
void precomputeTriangleNormals(const vector<Triangle>& triangles, vector<VectorFloatTriplet>& triangleNormals, const vector<VectorFloatTriplet>& vertices);
void precomputeCameraTriangleDeterminant(const Scene& scene, vector<vector<float>>& cameraTriangleDeterminant);
void precomputeCameraMeshDeterminant(const Scene& scene, vector<vector<vector<float>>>& cameraMeshDeterminant);

#endif