#include "bvh.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace scene {

namespace {

inline float getAxisValue(const VectorFloatTriplet& v, int axis) {
    switch (axis) {
        case 0: return v.x;
        case 1: return v.y;
        default: return v.z;
    }
}

} // anonymous namespace

bool AABB::intersect(const Ray& ray, float t_min, float t_max) const {
    float t0 = t_min;
    float t1 = t_max;

    const VectorFloatTriplet& origin = ray.origin;
    const VectorFloatTriplet& direction = ray.direction;

    for (int axis = 0; axis < 3; ++axis) {
        const float dirComponent = getAxisValue(direction, axis);
        const float originComponent = getAxisValue(origin, axis);
        const float slabMin = getAxisValue(min, axis);
        const float slabMax = getAxisValue(max, axis);

        if (std::fabs(dirComponent) < 1e-9f) {
            if (originComponent < slabMin || originComponent > slabMax) {
                return false;
            }
            continue;
        }

        float invD = 1.0f / dirComponent;
        float tNear = (slabMin - originComponent) * invD;
        float tFar = (slabMax - originComponent) * invD;
        if (tNear > tFar) std::swap(tNear, tFar);

        t0 = std::max(t0, tNear);
        t1 = std::min(t1, tFar);
        if (t1 < t0) {
            return false;
        }
    }

    return t1 >= t0;
}

void MeshBVH::build(const Mesh& mesh, const std::vector<VectorFloatTriplet>& vertices) {
    nodes.clear();
    faceIndices.clear();
    faceCentroids.clear();

    const size_t faceCount = mesh.faces.size();
    if (faceCount == 0) {
        return;
    }

    nodes.reserve(faceCount * 2);
    faceIndices.resize(faceCount);
    std::iota(faceIndices.begin(), faceIndices.end(), 0);

    faceCentroids.resize(faceCount);
    for (size_t i = 0; i < faceCount; ++i) {
        faceCentroids[i] = computeFaceCentroid(mesh.faces[i], vertices);
    }

    buildRecursive(mesh, vertices, 0, static_cast<int>(faceCount), 0);

    std::vector<VectorFloatTriplet>().swap(faceCentroids);
}

bool MeshBVH::traverse(
    Ray& ray,
    const Mesh& mesh,
    const std::vector<VectorFloatTriplet>& vertices,
    const std::vector<float>& determinants,
    float& t_min,
    Intersection& intersection,
    float intersectionTestEpsilon,
    bool enableBackFaceCulling,
    int meshIndex
) const {
    if (nodes.empty()) {
        return false;
    }

    bool hit = false;
    std::vector<int> stack;
    stack.reserve(64);
    stack.push_back(0);

    auto computeDeterminantT = [&](const VectorIntTriplet& face) -> float {
        const VectorFloatTriplet& a = vertices[face.x];
        const VectorFloatTriplet& b = vertices[face.y];
        const VectorFloatTriplet& c = vertices[face.z];
        const VectorFloatTriplet& o = ray.origin;

        const float e1x = b.x - a.x;
        const float e1y = b.y - a.y;
        const float e1z = b.z - a.z;
        const float e2x = c.x - a.x;
        const float e2y = c.y - a.y;
        const float e2z = c.z - a.z;
        const float rx = o.x - a.x;
        const float ry = o.y - a.y;
        const float rz = o.z - a.z;

        return e1x * (e2y * rz - e2z * ry)
             - e1y * (e2x * rz - e2z * rx)
             + e1z * (e2x * ry - e2y * rx);
    };

    while (!stack.empty()) {
        int nodeIndex = stack.back();
        stack.pop_back();

        const BVHNode& node = nodes[nodeIndex];
        if (!node.bounds.intersect(ray, 0.0f, t_min)) {
            continue;
        }

        if (node.isLeaf()) {
            for (int i = 0; i < node.faceCount; ++i) {
                const int faceIndex = faceIndices[node.firstFace + i];
                if (faceIndex < 0 || faceIndex >= static_cast<int>(mesh.faces.size())) {
                    continue;
                }

                const VectorIntTriplet& face = mesh.faces[faceIndex];
                float determinantT = 0.0f;
                if (faceIndex < static_cast<int>(determinants.size())) {
                    determinantT = determinants[faceIndex];
                } else {
                    determinantT = computeDeterminantT(face);
                }

                if (rayHitsTriangle(
                        ray,
                        face,
                        vertices,
                        t_min,
                        intersection,
                        intersectionTestEpsilon,
                        determinantT,
                        mesh.material,
                        enableBackFaceCulling,
                        meshIndex,
                        faceIndex)) {
                    hit = true;
                    intersection.kind = Intersection::Kind::Mesh;
                }
            }
        } else {
            const int leftChild = node.leftChild;
            const int rightChild = node.rightChild;

            if (leftChild < 0 || rightChild < 0) {
                continue;
            }

            const float dirComponent = getAxisValue(ray.direction, node.splitAxis);
            if (dirComponent >= 0.0f) {
                stack.push_back(rightChild);
                stack.push_back(leftChild);
            } else {
                stack.push_back(leftChild);
                stack.push_back(rightChild);
            }
        }
    }

    return hit;
}

int MeshBVH::buildRecursive(
    const Mesh& mesh,
    const std::vector<VectorFloatTriplet>& vertices,
    int start,
    int end,
    int depth
) {
    if (start >= end) {
        return -1;
    }

    const int nodeIndex = static_cast<int>(nodes.size());
    nodes.emplace_back();

    BVHNode& node = nodes[nodeIndex];
    node.bounds = computeBounds(mesh, vertices, start, end);
    node.firstFace = start;
    node.faceCount = end - start;
    node.leftChild = -1;
    node.rightChild = -1;
    node.splitAxis = 0;

    const int faceCount = end - start;
    if (faceCount <= MAX_FACES_PER_LEAF || depth >= MAX_BUILD_DEPTH) {
        return nodeIndex;
    }

    AABB centroidBounds;
    for (int i = start; i < end; ++i) {
        const int faceIdx = faceIndices[i];
        centroidBounds.expand(faceCentroids[faceIdx]);
    }

    const VectorFloatTriplet extent{
        centroidBounds.max.x - centroidBounds.min.x,
        centroidBounds.max.y - centroidBounds.min.y,
        centroidBounds.max.z - centroidBounds.min.z
    };

    int splitAxis = 0;
    float maxExtent = extent.x;
    if (extent.y > maxExtent) {
        splitAxis = 1;
        maxExtent = extent.y;
    }
    if (extent.z > maxExtent) {
        splitAxis = 2;
        maxExtent = extent.z;
    }

    if (maxExtent < 1e-5f) {
        return nodeIndex;
    }

    const int mid = start + faceCount / 2;
    auto comparator = [&](int lhs, int rhs) {
        return getAxisValue(faceCentroids[lhs], splitAxis) < getAxisValue(faceCentroids[rhs], splitAxis);
    };
    std::nth_element(faceIndices.begin() + start, faceIndices.begin() + mid, faceIndices.begin() + end, comparator);

    node.splitAxis = splitAxis;
    node.firstFace = -1;
    node.faceCount = 0;

    const int leftChild = buildRecursive(mesh, vertices, start, mid, depth + 1);
    const int rightChild = buildRecursive(mesh, vertices, mid, end, depth + 1);

    node.leftChild = leftChild;
    node.rightChild = rightChild;

    return nodeIndex;
}

AABB MeshBVH::computeBounds(
    const Mesh& mesh,
    const std::vector<VectorFloatTriplet>& vertices,
    int start,
    int end
) const {
    AABB bounds;
    for (int i = start; i < end; ++i) {
        const int faceIndex = faceIndices[i];
        if (faceIndex < 0 || faceIndex >= static_cast<int>(mesh.faces.size())) {
            continue;
        }

        const VectorIntTriplet& face = mesh.faces[faceIndex];
        bounds.expand(vertices[face.x]);
        bounds.expand(vertices[face.y]);
        bounds.expand(vertices[face.z]);
    }
    return bounds;
}

VectorFloatTriplet MeshBVH::computeFaceCentroid(
    const VectorIntTriplet& face,
    const std::vector<VectorFloatTriplet>& vertices
) const {
    const VectorFloatTriplet& v0 = vertices[face.x];
    const VectorFloatTriplet& v1 = vertices[face.y];
    const VectorFloatTriplet& v2 = vertices[face.z];
    const float inv3 = 1.0f / 3.0f;

    return VectorFloatTriplet{
        (v0.x + v1.x + v2.x) * inv3,
        (v0.y + v1.y + v2.y) * inv3,
        (v0.z + v1.z + v2.z) * inv3
    };
}

} // namespace scene

