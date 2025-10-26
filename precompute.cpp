#include "precompute.h"
#include "scene.h"
#include "overloads.h"

using namespace std;
using namespace scene;


void precomputeMeshNormals(
    const vector<Mesh>& meshes, 
    vector<vector<VectorFloatTriplet>>& meshNormals, 
    const vector<VectorFloatTriplet>& vertices
) {
    meshNormals.reserve(meshes.size());
    for (const Mesh& mesh : meshes) {
        meshNormals.push_back(vector<VectorFloatTriplet>());
        for (const VectorIntTriplet& face : mesh.faces) {
            VectorFloatTriplet v0 = vertices[face.x];
            VectorFloatTriplet v1 = vertices[face.y];
            VectorFloatTriplet v2 = vertices[face.z];
            VectorFloatTriplet v1_v0 = v1 - v0;
            VectorFloatTriplet v2_v0 = v2 - v0;
            VectorFloatTriplet normal = crossProduct(v1_v0, v2_v0);
            normal = normalize(normal);
            meshNormals.back().push_back(normal);
        }
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
    vector<vector<float>>& cameraTriangleDeterminant
) {
    cameraTriangleDeterminant.reserve(scene.cameras.size());
    for (const Camera& camera : scene.cameras) {
        cameraTriangleDeterminant.push_back(vector<float>());
        cameraTriangleDeterminant.back().reserve(scene.triangles.size());
        for (const Triangle& triangle : scene.triangles) {
            VectorFloatTriplet o = camera.position;
            VectorFloatTriplet a = scene.vertices[triangle.indices.x];
            VectorFloatTriplet b = scene.vertices[triangle.indices.y];
            VectorFloatTriplet c = scene.vertices[triangle.indices.z];
            float ax = a.x, ay = a.y, az = a.z;
            float bx = b.x, by = b.y, bz = b.z;
            float cx = c.x, cy = c.y, cz = c.z;
            float ox = o.x, oy = o.y, oz = o.z;
            const float e1x = bx - ax, e1y = by - ay, e1z = bz - az;
            const float e2x = cx - ax, e2y = cy - ay, e2z = cz - az;
            const float rx  = ox - ax, ry  = oy - ay, rz  = oz - az;
            cameraTriangleDeterminant.back().push_back(e1x * (e2y * rz - e2z * ry)
                                                    - e1y * (e2x * rz - e2z * rx)
                                                    + e1z * (e2x * ry - e2y * rx));
        }
    }
}

void precomputeCameraMeshDeterminant(
    const Scene& scene, 
    vector<vector<vector<float>>>& cameraMeshDeterminant
) {
    cameraMeshDeterminant.reserve(scene.cameras.size());
    for (const Camera& camera : scene.cameras) {
        cameraMeshDeterminant.push_back(vector<vector<float>>());
        cameraMeshDeterminant.back().reserve(scene.meshes.size());
        for (const Mesh& mesh : scene.meshes) {
            cameraMeshDeterminant.back().emplace_back(vector<float>());
            cameraMeshDeterminant.back().back().reserve(mesh.faces.size());
            for (const VectorIntTriplet& face : mesh.faces) {
                VectorFloatTriplet v0 = scene.vertices[face.x];
                VectorFloatTriplet v1 = scene.vertices[face.y];
                VectorFloatTriplet v2 = scene.vertices[face.z];
                float ax = v0.x, ay = v0.y, az = v0.z;
                float bx = v1.x, by = v1.y, bz = v1.z;
                float cx = v2.x, cy = v2.y, cz = v2.z;
                float ox = camera.position.x, oy = camera.position.y, oz = camera.position.z;
                const float e1x = bx - ax, e1y = by - ay, e1z = bz - az;
                const float e2x = cx - ax, e2y = cy - ay, e2z = cz - az;
                const float rx  = ox - ax, ry  = oy - ay, rz  = oz - az;
                cameraMeshDeterminant.back().back().push_back(e1x * (e2y * rz - e2z * ry)
                                                                - e1y * (e2x * rz - e2z * rx)
                                                                + e1z * (e2x * ry - e2y * rx));
            }
        }
    }
    return;
}