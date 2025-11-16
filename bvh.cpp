#include "bvh.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace scene {

namespace {

inline double getAxisValue(const VectorFloatTriplet& v, int axis) {
    switch (axis) {
        case 0: return v.x;
        case 1: return v.y;
        default: return v.z;
    }
}

} // anonymous namespace

bool AABB::intersect(const Ray& ray, double t_min, double t_max) const {
    double t0 = t_min;
    double t1 = t_max;

    const VectorFloatTriplet& origin = ray.origin;
    const VectorFloatTriplet& direction = ray.direction;

    for (int axis = 0; axis < 3; ++axis) {
        const double dirComponent = getAxisValue(direction, axis);
        const double originComponent = getAxisValue(origin, axis);
        const double slabMin = getAxisValue(min, axis);
        const double slabMax = getAxisValue(max, axis);

        if (std::fabs(dirComponent) < 1e-9f) {
            if (originComponent < slabMin || originComponent > slabMax) {
                return false;
            }
            continue;
        }

        double invD = 1.0 / dirComponent;
        double tNear = (slabMin - originComponent) * invD;
        double tFar = (slabMax - originComponent) * invD;
        if (tNear > tFar) std::swap(tNear, tFar);

        t0 = std::max(t0, tNear);
        t1 = std::min(t1, tFar);
        if (t1 < t0) {
            return false;
        }
    }

    return t1 >= t0;
}

AABB AABB::transform(const Matrix4x4& mat) const {
    // Transform all 8 corners of the bounding box
    VectorFloatTriplet corners[8] = {
        {min.x, min.y, min.z},
        {min.x, min.y, max.z},
        {min.x, max.y, min.z},
        {min.x, max.y, max.z},
        {max.x, min.y, min.z},
        {max.x, min.y, max.z},
        {max.x, max.y, min.z},
        {max.x, max.y, max.z}
    };
    
    AABB result;
    for (int i = 0; i < 8; i++) {
        VectorFloatTriplet transformed = transformPoint(mat, corners[i]);
        result.expand(transformed);
    }
    
    return result;
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
    const std::vector<double>& determinants,
    double& t_min,
    Intersection& intersection,
    double intersectionTestEpsilon,
    bool enableBackFaceCulling,
    int meshIndex,
    Material* materialOverride
) const {
    if (nodes.empty()) {
        return false;
    }
    
    // Use material override if provided (for instances), otherwise use mesh's material
    Material* materialToUse = materialOverride ? materialOverride : mesh.material;

    bool hit = false;
    int stack[64];
    int stackSize = 0;
    stack[stackSize++] = 0;

    while (stackSize > 0) {
        int nodeIndex = stack[--stackSize];

        const BVHNode& node = nodes[nodeIndex];
        if (!node.bounds.intersect(ray, 0.0, t_min)) {
            continue;
        }

        if (node.isLeaf()) {
            for (int i = 0; i < node.faceCount; ++i) {
                const int faceIndex = faceIndices[node.firstFace + i];
                if (faceIndex < 0 || faceIndex >= static_cast<int>(mesh.faces.size())) {
                    continue;
                }

                const VectorIntTriplet& face = mesh.faces[faceIndex];
                double determinantT = 0.0;
                if (faceIndex < static_cast<int>(determinants.size()) && determinants[faceIndex] != 0.0) {
                    determinantT = determinants[faceIndex];
                }

                if (rayHitsTriangle(
                        ray,
                        face,
                        vertices,
                        t_min,
                        intersection,
                        intersectionTestEpsilon,
                        determinantT,
                        materialToUse,
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

            const double dirComponent = getAxisValue(ray.direction, node.splitAxis);
            if (dirComponent >= 0.0) {
                stack[stackSize++] = rightChild;
                stack[stackSize++] = leftChild;
            } else {
                stack[stackSize++] = leftChild;
                stack[stackSize++] = rightChild;
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
    double maxExtent = extent.x;
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
    const double inv3 = 1.0 / 3.0;

    return VectorFloatTriplet{
        (v0.x + v1.x + v2.x) * inv3,
        (v0.y + v1.y + v2.y) * inv3,
        (v0.z + v1.z + v2.z) * inv3
    };
}

AABB computeMeshAABB(const Mesh& mesh, const std::vector<VectorFloatTriplet>& vertices) {
    AABB bounds;
    for (const auto& face : mesh.faces) {
        bounds.expand(vertices[face.x]);
        bounds.expand(vertices[face.y]);
        bounds.expand(vertices[face.z]);
    }
    return bounds;
}

} // namespace scene

