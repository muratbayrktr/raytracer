#include "precompute.h"
#include "scene.h"
#include "overloads.h"

using namespace std;
using namespace scene;


void precomputeMeshNormals(
    const vector<Mesh>& meshes, 
    vector<vector<VectorFloatTriplet>>& meshVertexNormals,  // per-vertex normals
    const vector<VectorFloatTriplet>& vertices
) {
    meshVertexNormals.clear();
    meshVertexNormals.reserve(meshes.size());

    auto safeNormalize = [](const VectorFloatTriplet& v) -> VectorFloatTriplet {
        double len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
        if (len > 0.0) return v * (1.0 / len);
        return VectorFloatTriplet{0.0, 0.0, 1.0};
    };
    for (size_t meshIdx = 0; meshIdx < meshes.size(); ++meshIdx) {
        const Mesh& mesh = meshes[meshIdx];
        vector<VectorFloatTriplet> vertexAccum(vertices.size(), VectorFloatTriplet{0,0,0});
        for (const VectorIntTriplet& face : mesh.faces) {
            const int i0 = face.x, i1 = face.y, i2 = face.z;

            const VectorFloatTriplet& v0 = vertices[i0];
            const VectorFloatTriplet& v1 = vertices[i1];
            const VectorFloatTriplet& v2 = vertices[i2];
            const VectorFloatTriplet e1 = v1 - v0;
            const VectorFloatTriplet e2 = v2 - v0;
            const VectorFloatTriplet n_unnormalized = crossProduct(e1, e2);
            vertexAccum[i0] = vertexAccum[i0] + n_unnormalized;
            vertexAccum[i1] = vertexAccum[i1] + n_unnormalized;
            vertexAccum[i2] = vertexAccum[i2] + n_unnormalized;
        }
        vector<VectorFloatTriplet> vertexNormals(vertices.size());
        for (size_t i = 0; i < vertices.size(); ++i) {
            vertexNormals[i] = safeNormalize(vertexAccum[i]);
        }
        meshVertexNormals.push_back(std::move(vertexNormals));
    }
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
        VectorFloatTriplet v1_v0 = v1 - v0;
        VectorFloatTriplet v2_v0 = v2 - v0;
        VectorFloatTriplet normal = crossProduct(v1_v0, v2_v0);
        normal = normalize(normal);
        triangleNormals.push_back(normal);
    }
    return;
}

void precomputeCameraTriangleDeterminant(
    const Scene& scene, 
    vector<vector<double>>& cameraTriangleDeterminant
) {
    cameraTriangleDeterminant.reserve(scene.cameras.size());
    for (const Camera& camera : scene.cameras) {
        cameraTriangleDeterminant.push_back(vector<double>());
        cameraTriangleDeterminant.back().reserve(scene.triangles.size());
        for (const Triangle& triangle : scene.triangles) {
            VectorFloatTriplet o = camera.position;
            VectorFloatTriplet a = scene.vertices[triangle.indices.x];
            VectorFloatTriplet b = scene.vertices[triangle.indices.y];
            VectorFloatTriplet c = scene.vertices[triangle.indices.z];
            double ax = a.x, ay = a.y, az = a.z;
            double bx = b.x, by = b.y, bz = b.z;
            double cx = c.x, cy = c.y, cz = c.z;
            double ox = o.x, oy = o.y, oz = o.z;
            const double e1x = bx - ax, e1y = by - ay, e1z = bz - az;
            const double e2x = cx - ax, e2y = cy - ay, e2z = cz - az;
            const double rx  = ox - ax, ry  = oy - ay, rz  = oz - az;
            cameraTriangleDeterminant.back().push_back(e1x * (e2y * rz - e2z * ry)
                                                    - e1y * (e2x * rz - e2z * rx)
                                                    + e1z * (e2x * ry - e2y * rx));
        }
    }
}

void precomputeCameraMeshDeterminant(
    const Scene& scene, 
    vector<vector<vector<double>>>& cameraMeshDeterminant
) {
    cameraMeshDeterminant.reserve(scene.cameras.size());
    for (const Camera& camera : scene.cameras) {
        cameraMeshDeterminant.push_back(vector<vector<double>>());
        cameraMeshDeterminant.back().reserve(scene.meshes.size());
        for (const Mesh& mesh : scene.meshes) {
            cameraMeshDeterminant.back().emplace_back(vector<double>());
            cameraMeshDeterminant.back().back().reserve(mesh.faces.size());
            for (const VectorIntTriplet& face : mesh.faces) {
                VectorFloatTriplet v0 = scene.vertices[face.x];
                VectorFloatTriplet v1 = scene.vertices[face.y];
                VectorFloatTriplet v2 = scene.vertices[face.z];
                double ax = v0.x, ay = v0.y, az = v0.z;
                double bx = v1.x, by = v1.y, bz = v1.z;
                double cx = v2.x, cy = v2.y, cz = v2.z;
                double ox = camera.position.x, oy = camera.position.y, oz = camera.position.z;
                const double e1x = bx - ax, e1y = by - ay, e1z = bz - az;
                const double e2x = cx - ax, e2y = cy - ay, e2z = cz - az;
                const double rx  = ox - ax, ry  = oy - ay, rz  = oz - az;
                cameraMeshDeterminant.back().back().push_back(e1x * (e2y * rz - e2z * ry)
                                                                - e1y * (e2x * rz - e2z * rx)
                                                                + e1z * (e2x * ry - e2y * rx));
            }
        }
    }
    return;
}