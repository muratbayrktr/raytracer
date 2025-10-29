#include "bvh.h"
#include "utils.h"
#include <iostream>
#include <limits>
#include <algorithm>

using namespace scene;

bool AABB::intersect(const Ray& ray, float tMin, float tMax) const {
    /*
    From the lecture slides:
        Repeat this process for each slab
        Find the maximum t1 (entrance parameter)
        Find the minimum t2 (exit parameter)
        If t1 > t2 then the ray does not intersect
        If t1 < t2 then the ray enters at r(t1) and exits at r(t2)
    */
    float t_min = tMin;
    float t_max = tMax;
    
    for (int axis = 0; axis < 3; axis++) {
        float origin, direction, bmin, bmax;
        
        if (axis == 0) {
            origin = ray.origin.x;
            direction = ray.direction.x;
            bmin = min.x;
            bmax = max.x;
        } else if (axis == 1) {
            origin = ray.origin.y;
            direction = ray.direction.y;
            bmin = min.y;
            bmax = max.y;
        } else {
            origin = ray.origin.z;
            direction = ray.direction.z;
            bmin = min.z;
            bmax = max.z;
        }
        
        float invD = 1.0f / direction;
        float t0 = (bmin - origin) * invD;
        float t1 = (bmax - origin) * invD;
        
        if (invD < 0.0f) {
            std::swap(t0, t1);
        }
        
        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;
        
        if (t_max <= t_min) {
            return false;
        }
    }
    
    return true;
}

VectorFloatTriplet MeshBVH::computeFaceCentroid(
    const VectorIntTriplet& face,
    const std::vector<VectorFloatTriplet>& vertices
) const {
    VectorFloatTriplet a = vertices[face.x];
    VectorFloatTriplet b = vertices[face.y];
    VectorFloatTriplet c = vertices[face.z];
    
    return VectorFloatTriplet{
        (a.x + b.x + c.x) / 3.0f,
        (a.y + b.y + c.y) / 3.0f,
        (a.z + b.z + c.z) / 3.0f
    };
}

AABB MeshBVH::computeBounds(
    const Mesh& mesh,
    const std::vector<VectorFloatTriplet>& vertices,
    int start,
    int end
) const {
    AABB bounds;
    
    for (int i = start; i < end; i++) {
        const VectorIntTriplet& face = mesh.faces[faceIndices[i]];
        bounds.expand(vertices[face.x]);
        bounds.expand(vertices[face.y]);
        bounds.expand(vertices[face.z]);
    }
    
    return bounds;
}

void MeshBVH::build(const Mesh& mesh, const std::vector<VectorFloatTriplet>& vertices) {
    if (mesh.faces.empty()) {
        return;
    }
    
    faceIndices.resize(mesh.faces.size());
    for (size_t i = 0; i < mesh.faces.size(); i++) {
        faceIndices[i] = i;
    }
    
    nodes.reserve(mesh.faces.size() * 2);
    
    buildRecursive(mesh, vertices, 0, mesh.faces.size(), 0);
    
    std::cout << "[BVH] Built BVH with " << nodes.size() << " nodes for " 
              << mesh.faces.size() << " faces" << std::endl;
}

int MeshBVH::buildRecursive(
    const Mesh& mesh,
    const std::vector<VectorFloatTriplet>& vertices,
    int start,
    int end,
    int depth
) {
    /*
        I will apply the algorithm verbatim from the lecture slides:
            Step 1: Compute the global bounding box of all objects
            Step 2: Decide on a splitting axis
                May start with the largest extent (assume x) and progress as y, z, back to x and so on
            Step 3: Decide on a splitting point on the axis
                May choose the middle for simplicity
                Or sort the objects and choose the median object’s centroid
                Equal count heuristic
                Choose the point based on bounding box areas and object counts
                Surface area heuristic (SAH)
            Step 4: Assign objects to each part
                Can be done by swapping the elements in the container (e.g. vector)
    */
    int nodeIndex = nodes.size();
    nodes.push_back(BVHNode());
    BVHNode& node = nodes[nodeIndex];
    node.bounds = computeBounds(mesh, vertices, start, end);
    
    int faceCount = end - start;

    const int MAX_LEAF_SIZE = 4;
    const int MAX_DEPTH = 32;
    
    if (faceCount <= MAX_LEAF_SIZE || depth >= MAX_DEPTH) {
        node.firstFace = start;
        node.faceCount = faceCount;
        return nodeIndex;
    }

    int bestAxis = 0;
    int bestSplit = start + faceCount / 2;
    float bestCost = std::numeric_limits<float>::max();
    

    for (int axis = 0; axis < 3; axis++) {
        std::sort(faceIndices.begin() + start, faceIndices.begin() + end,
            [&](int a, int b) {
                VectorFloatTriplet centroidA = computeFaceCentroid(mesh.faces[a], vertices);
                VectorFloatTriplet centroidB = computeFaceCentroid(mesh.faces[b], vertices);

                if (axis == 0) return centroidA.x < centroidB.x;
                else if (axis == 1) return centroidA.y < centroidB.y;
                else return centroidA.z < centroidB.z;
            }
        );

        const int NUM_BUCKETS = 8;
        for (int i = 1; i < NUM_BUCKETS; i++) {
            int splitPos = start + (faceCount * i) / NUM_BUCKETS;


            AABB leftBounds = computeBounds(mesh, vertices, start, splitPos);
            AABB rightBounds = computeBounds(mesh, vertices, splitPos, end);
            
            int leftCount = splitPos - start;
            int rightCount = end - splitPos;
            
            float cost = leftBounds.surfaceArea() * leftCount + rightBounds.surfaceArea() * rightCount;
            
            if (cost < bestCost) {
                bestCost = cost;
                bestAxis = axis;
                bestSplit = splitPos;
            }
        }
    }
    std::sort(faceIndices.begin() + start, faceIndices.begin() + end,
        [&](int a, int b) {
            VectorFloatTriplet centroidA = computeFaceCentroid(mesh.faces[a], vertices);
            VectorFloatTriplet centroidB = computeFaceCentroid(mesh.faces[b], vertices);
            
            if (bestAxis == 0) return centroidA.x < centroidB.x;
            else if (bestAxis == 1) return centroidA.y < centroidB.y;
            else return centroidA.z < centroidB.z;
        }
    );
    
    bestSplit = start + faceCount / 2;
    
    if (bestSplit == start || bestSplit == end) {
        bestSplit = start + faceCount / 2;
    }

    node.leftChild = buildRecursive(mesh, vertices, start, bestSplit, depth + 1);
    node.rightChild = buildRecursive(mesh, vertices, bestSplit, end, depth + 1);
    
    return nodeIndex;
}

bool MeshBVH::traverse(
    Ray& ray,
    const Mesh& mesh,
    const std::vector<VectorFloatTriplet>& normals,
    const std::vector<VectorFloatTriplet>& vertices,
    const std::vector<float>& determinants,
    float& t_min,
    Intersection& intersection,
    float intersectionTestEpsilon,
    bool enableBackFaceCulling
) const {
    if (nodes.empty()) {
        return false;
    }
    
    bool hit = false;
    
    const int MAX_STACK_SIZE = 64;
    int stack[MAX_STACK_SIZE];
    int stackPtr = 0;
    
    stack[stackPtr++] = 0;
    
    while (stackPtr > 0) {
        int nodeIndex = stack[--stackPtr];
        const BVHNode& node = nodes[nodeIndex];
        
        if (!node.bounds.intersect(ray, 0.0f, std::numeric_limits<float>::max())) {
            continue;
        }
        
        if (node.isLeaf()) {
            for (int i = 0; i < node.faceCount; i++) {
                int faceIdx = faceIndices[node.firstFace + i];
                
                if (rayHitsTriangle(
                    ray,
                    mesh.faces[faceIdx],
                    normals[faceIdx],
                    vertices,
                    t_min,
                    intersection,
                    intersectionTestEpsilon,
                    determinants[faceIdx],
                    mesh.material,
                    enableBackFaceCulling
                )) {
                    hit = true;
                }
            }
        } else {
            if (stackPtr < MAX_STACK_SIZE - 2) {
                stack[stackPtr++] = node.rightChild;
                stack[stackPtr++] = node.leftChild;
            }
        }
    }
    
    return hit;
}

