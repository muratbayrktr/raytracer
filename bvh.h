#ifndef __BVH__
#define __BVH__

#include "scene.h"
#include <vector>
#include <algorithm>

namespace scene {
    struct AABB {
        VectorFloatTriplet min;
        VectorFloatTriplet max;
        
        AABB() : min{1e30f, 1e30f, 1e30f}, max{-1e30f, -1e30f, -1e30f} {}
        void expand(const VectorFloatTriplet& p) {
            min.x = std::min(min.x, p.x);
            min.y = std::min(min.y, p.y);
            min.z = std::min(min.z, p.z);
            max.x = std::max(max.x, p.x);
            max.y = std::max(max.y, p.y);
            max.z = std::max(max.z, p.z);
        }
        
        void expand(const AABB& other) {
            expand(other.min);
            expand(other.max);
        }
        
        VectorFloatTriplet center() const {
            return VectorFloatTriplet{
                (min.x + max.x) * 0.5f,
                (min.y + max.y) * 0.5f,
                (min.z + max.z) * 0.5f
            };
        }
        
        float surfaceArea() const {
            float dx = max.x - min.x;
            float dy = max.y - min.y;
            float dz = max.z - min.z;
            return 2.0f * (dx * dy + dy * dz + dz * dx);
        }
        
        bool intersect(const Ray& ray, float t_min, float t_max) const;
    };


    struct BVHNode {
        AABB bounds;
        int leftChild;
        int rightChild;
        int firstFace;
        int faceCount;
        
        BVHNode() : leftChild(-1), rightChild(-1), firstFace(0), faceCount(0) {}
        
        bool isLeaf() const { return leftChild == -1; }
    };

    class MeshBVH {
    public:
        std::vector<BVHNode> nodes;
        std::vector<int> faceIndices;
        
        MeshBVH() {}
        
        void build(const Mesh& mesh, const std::vector<VectorFloatTriplet>& vertices);
        
        bool traverse(
            Ray& ray,
            const Mesh& mesh,
            const std::vector<VectorFloatTriplet>& vertices,
            const std::vector<float>& determinants,
            float& t_min,
            Intersection& intersection,
            float intersectionTestEpsilon,
            bool enableBackFaceCulling,
            int meshIndex
        ) const;
        
    private:
        int buildRecursive(
            const Mesh& mesh,
            const std::vector<VectorFloatTriplet>& vertices,
            int start,
            int end,
            int depth
        );
        
        AABB computeBounds(
            const Mesh& mesh,
            const std::vector<VectorFloatTriplet>& vertices,
            int start,
            int end
        ) const;
        
        VectorFloatTriplet computeFaceCentroid(
            const VectorIntTriplet& face,
            const std::vector<VectorFloatTriplet>& vertices
        ) const;
    };
}
#endif

