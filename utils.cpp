#include <math.h>
#include <fstream>
#include <stdexcept>
#include "scene.h"
#include "utils.h"

using namespace std;
using namespace scene;


void precomputeMeshNormals(
    const vector<Mesh>& meshes, 
    vector<VectorFloatTriplet>& meshNormals, 
    const vector<VectorFloatTriplet>& vertices
) {
    meshNormals.reserve(meshes.size());
    for (const Mesh& mesh : meshes) {
        VectorFloatTriplet normal = {0, 0, 0};
        for (const VectorIntTriplet& face : mesh.faces) {
            VectorFloatTriplet v0 = vertices[face.x];
            VectorFloatTriplet v1 = vertices[face.y];
            VectorFloatTriplet v2 = vertices[face.z];
            VectorFloatTriplet normal = crossProduct(v1 - v0, v2 - v0);
            normal = normalize(normal);
            normal += normal;
        }
        meshNormals.push_back(normal);
    }
    return;
}

void precomputeTriangleNormals(
    const vector<Triangle>& triangles, 
    vector<VectorFloatTriplet>& triangleNormals, 
    const vector<VectorFloatTriplet>& vertices
) {
    triangleNormals.reserve(triangles.size());
    for (const Triangle& triangle : triangles) {
        VectorFloatTriplet v0 = vertices[triangle.indices.x];
        VectorFloatTriplet v1 = vertices[triangle.indices.y];
        VectorFloatTriplet v2 = vertices[triangle.indices.z];
        VectorFloatTriplet v1_v0 = {v1.x - v0.x, v1.y - v0.y, v1.z - v0.z};
        VectorFloatTriplet v2_v0 = {v2.x - v0.x, v2.y - v0.y, v2.z - v0.z};
        VectorFloatTriplet normal = crossProduct(v1_v0, v2_v0);
        normal = normalize(normal);
        triangleNormals.push_back(normal);
    }
    return;
}

void writePPM(const string& filename, const vector<unsigned char>& image, int width, int height) {
    FILE *outfile;
    if ((outfile = fopen(filename.c_str(), "w")) == NULL) {
        throw runtime_error("Error: The ppm file cannot be opened for writing.");
    }
    fprintf(outfile, "P3\n%d %d\n255\n", width, height);
    for (int i = 0; i < width * height * 3; i++) {
        fprintf(outfile, "%d ", image[i]);
    }
    fclose(outfile);
}

/* Operator Overloads */
/* VectorFloatTriplet */
VectorFloatTriplet crossProduct(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; 
}

float dotProduct(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return a.x * b.x + a.y * b.y + a.z * b.z; 
}

VectorFloatTriplet normalize(const VectorFloatTriplet& a) { 
    float length = sqrt(a.x * a.x + a.y * a.y + a.z * a.z); 
    return VectorFloatTriplet{a.x / length, a.y / length, a.z / length}; 
}

VectorFloatTriplet operator-(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.x - b.x, a.y - b.y, a.z - b.z}; 
}

VectorFloatTriplet operator+(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.x + b.x, a.y + b.y, a.z + b.z}; 
}

VectorFloatTriplet operator*(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.x * b.x, a.y * b.y, a.z * b.z}; 
}

VectorFloatTriplet operator+=(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.x + b.x, a.y + b.y, a.z + b.z}; 
}

VectorFloatTriplet operator-=(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.x - b.x, a.y - b.y, a.z - b.z}; 
}

VectorFloatTriplet operator*=(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.x * b.x, a.y * b.y, a.z * b.z}; 
}

VectorFloatTriplet operator/=(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.x / b.x, a.y / b.y, a.z / b.z}; 
}

/* VectorIntTriplet */
VectorIntTriplet crossProduct(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; 
}

int dotProduct(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return a.x * b.x + a.y * b.y + a.z * b.z; 
}

VectorIntTriplet normalize(const VectorIntTriplet& a) { 
    int length = sqrt(a.x * a.x + a.y * a.y + a.z * a.z); 
    return VectorIntTriplet{a.x / length, a.y / length, a.z / length}; 
}

VectorIntTriplet operator-(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.x - b.x, a.y - b.y, a.z - b.z}; 
}

VectorIntTriplet operator+(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.x + b.x, a.y + b.y, a.z + b.z}; 
}

VectorIntTriplet operator*(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.x * b.x, a.y * b.y, a.z * b.z}; 
}

VectorIntTriplet operator+=(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.x + b.x, a.y + b.y, a.z + b.z}; 
}

VectorIntTriplet operator-=(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.x - b.x, a.y - b.y, a.z - b.z}; 
}

VectorIntTriplet operator*=(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.x * b.x, a.y * b.y, a.z * b.z}; 
}

VectorIntTriplet operator/=(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.x / b.x, a.y / b.y, a.z / b.z}; 
}