#include <math.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <chrono>
#include "scene.h"
#include "utils.h"
#include "overloads.h"
#include "bvh.h"

using namespace std;
using namespace scene;

#define PROFILE_PERF 0

// Forward declare optimized inline transformation functions
static inline VectorFloatTriplet transformPointFast(const double* mat, double px, double py, double pz);
static inline VectorFloatTriplet transformDirectionFast(const double* mat, double dx, double dy, double dz);
static inline VectorFloatTriplet transformNormalFast(const double* mat, double nx, double ny, double nz);

#if PROFILE_PERF
#include <atomic>
std::atomic<long> g_intersectCalls(0);
std::atomic<long> g_bvhTraversals(0);
std::atomic<long> g_triangleTests(0);
std::atomic<long> g_transformedMeshCalls(0);
std::atomic<long> g_worldBoundsRejects(0);
std::atomic<long> g_worldBoundsAccepts(0);
std::atomic<long> g_rayTransforms(0);

// Timing atomics (in nanoseconds)
std::atomic<long long> g_timeIntersect(0);
std::atomic<long long> g_timeRayHitsMesh(0);
std::atomic<long long> g_timeBVHTraverse(0);
std::atomic<long long> g_timeRayTransform(0);
std::atomic<long long> g_timeAABBTest(0);
std::atomic<long long> g_timeTriangleTest(0);
std::atomic<long long> g_timeBackTransform(0);
std::atomic<long long> g_timeIntersectPlanes(0);
std::atomic<long long> g_timeIntersectTriangles(0);
std::atomic<long long> g_timeIntersectMeshes(0);
std::atomic<long long> g_timeIntersectSpheres(0);
std::atomic<long long> g_timeIntersectInstances(0);
std::atomic<long long> g_timeIntersectPostProcess(0);

void printPerfStats() {
    std::cout << "\n=== Performance Stats ===" << std::endl;
    std::cout << "Intersect calls: " << g_intersectCalls.load() << std::endl;
    std::cout << "BVH traversals: " << g_bvhTraversals.load() << std::endl;
    std::cout << "Triangle tests: " << g_triangleTests.load() << std::endl;
    std::cout << "Transformed mesh calls: " << g_transformedMeshCalls.load() << std::endl;
    std::cout << "World bounds rejects: " << g_worldBoundsRejects.load() << std::endl;
    std::cout << "World bounds accepts: " << g_worldBoundsAccepts.load() << std::endl;
    std::cout << "Ray transforms: " << g_rayTransforms.load() << std::endl;
    double rejectRate = 100.0 * g_worldBoundsRejects.load() / (g_worldBoundsRejects.load() + g_worldBoundsAccepts.load() + 1);
    std::cout << "World bounds reject rate: " << rejectRate << "%" << std::endl;
    
    std::cout << "\n=== Timing Breakdown (ms) ===" << std::endl;
    std::cout << "Time in intersect(): " << g_timeIntersect.load() / 1000000.0 << " ms" << std::endl;
    std::cout << "  - Planes: " << g_timeIntersectPlanes.load() / 1000000.0 << " ms" << std::endl;
    std::cout << "  - Triangles: " << g_timeIntersectTriangles.load() / 1000000.0 << " ms" << std::endl;
    std::cout << "  - Meshes: " << g_timeIntersectMeshes.load() / 1000000.0 << " ms" << std::endl;
    std::cout << "  - Spheres: " << g_timeIntersectSpheres.load() / 1000000.0 << " ms" << std::endl;
    std::cout << "  - Instances: " << g_timeIntersectInstances.load() / 1000000.0 << " ms" << std::endl;
    std::cout << "  - PostProcess: " << g_timeIntersectPostProcess.load() / 1000000.0 << " ms" << std::endl;
    std::cout << "Time in rayHitsMesh(): " << g_timeRayHitsMesh.load() / 1000000.0 << " ms" << std::endl;
    std::cout << "  - BVH traverse: " << g_timeBVHTraverse.load() / 1000000.0 << " ms" << std::endl;
    std::cout << "  - Ray transform: " << g_timeRayTransform.load() / 1000000.0 << " ms" << std::endl;
    std::cout << "  - AABB tests: " << g_timeAABBTest.load() / 1000000.0 << " ms" << std::endl;
    std::cout << "  - Triangle tests: " << g_timeTriangleTest.load() / 1000000.0 << " ms" << std::endl;
    std::cout << "  - Back transform: " << g_timeBackTransform.load() / 1000000.0 << " ms" << std::endl;
    std::cout << "=========================\n" << std::endl;
}

void printPerfStatsInline() {
    long long timeIntersect = g_timeIntersect.load();
    long long timeRayHitsMesh = g_timeRayHitsMesh.load();
    long long timeBVH = g_timeBVHTraverse.load();
    long long timeTransform = g_timeRayTransform.load();
    long long timeAABB = g_timeAABBTest.load();
    long long timeBackTransform = g_timeBackTransform.load();
    
    std::cout << " | Timing(ms): Intersect=" << timeIntersect/1e6 
              << " RayHitsMesh=" << timeRayHitsMesh/1e6
              << " (BVH=" << timeBVH/1e6 
              << " Transform=" << timeTransform/1e6
              << " AABB=" << timeAABB/1e6
              << " BackTransform=" << timeBackTransform/1e6 << ")";
}
#else
void printPerfStats() {}
void printPerfStatsInline() {}
#endif

// Optimized inline transformation functions - defined early for use throughout
static inline VectorFloatTriplet transformPointFast(const double* mat, double px, double py, double pz) {
    double x = mat[0]*px + mat[1]*py + mat[2]*pz + mat[3];
    double y = mat[4]*px + mat[5]*py + mat[6]*pz + mat[7];
    double z = mat[8]*px + mat[9]*py + mat[10]*pz + mat[11];
    double w = mat[12]*px + mat[13]*py + mat[14]*pz + mat[15];
    if (w > 1.00001f || w < 0.99999f) {
        double invW = 1.0 / w;
        return VectorFloatTriplet{x*invW, y*invW, z*invW};
    }
    return VectorFloatTriplet{x, y, z};
}

static inline VectorFloatTriplet transformDirectionFast(const double* mat, double dx, double dy, double dz) {
    return VectorFloatTriplet{
        mat[0]*dx + mat[1]*dy + mat[2]*dz,
        mat[4]*dx + mat[5]*dy + mat[6]*dz,
        mat[8]*dx + mat[9]*dy + mat[10]*dz
    };
}

static inline VectorFloatTriplet transformNormalFast(const double* mat, double nx, double ny, double nz) {
    // Use row-major access (NOT transposed) because normalMatrix is already transpose(inverse(M))
    return VectorFloatTriplet{
        mat[0]*nx + mat[1]*ny + mat[2]*nz,
        mat[4]*nx + mat[5]*ny + mat[6]*nz,
        mat[8]*nx + mat[9]*ny + mat[10]*nz
    };
}

void clamp(VectorFloatTriplet& color, int min, int max) {
    if (color.x < min) color.x = min;
    if (color.x > max) color.x = max;
    if (color.y < min) color.y = min;
    if (color.y > max) color.y = max;
    if (color.z < min) color.z = min;
    if (color.z > max) color.z = max;
}

Ray castRay(const Camera& camera, int x, int y, int width, int height) {
    VectorFloatTriplet w = -normalize(camera.gaze);
    VectorFloatTriplet v = normalize(camera.up);
    VectorFloatTriplet u = crossProduct(v, w);
    double l = camera.nearPlane.x;
    double r = camera.nearPlane.y;
    double b = camera.nearPlane.z;
    double t = camera.nearPlane.w;
    double s_u = (x+0.5)*(r - l) / width;
    double s_v = (y+0.5)*(t - b) / height;
    VectorFloatTriplet e = camera.position;
    VectorFloatTriplet m = e - w * camera.nearDistance;
    VectorFloatTriplet q = m + l*u + t*v;
    VectorFloatTriplet s = q + u*s_u - v*s_v;
    /* 
        From lecture slides:
        r(t) = e + t * (s - e) 
        ray_direction = s - e
        origin = e
    */
    VectorFloatTriplet ray_direction = s - e;
    VectorFloatTriplet origin = e;
    Ray ray = Ray{origin, normalize(ray_direction), 0, false, false, false};
    return ray;
}

bool rayHitsPlane(
    Ray& ray, 
    const Plane& plane, 
    const vector<VectorFloatTriplet>& vertices, 
    double& t_min, 
    Intersection& intersection,
    int planeIndex
) {
       Ray objectRay = ray;
       VectorFloatTriplet objectNormal = plane.normal;
       VectorFloatTriplet objectPoint = vertices[plane.point];
       
       if (plane.hasTransformation) {
           objectRay.origin = transformPoint(*plane.inverseTransformMatrix, ray.origin);
           objectRay.direction = normalize(transformDirection(*plane.inverseTransformMatrix, ray.direction));
       }
       
       VectorFloatTriplet n = objectNormal;
       VectorFloatTriplet a = objectPoint;
       VectorFloatTriplet d = objectRay.direction;
       VectorFloatTriplet o = objectRay.origin;
       double denom = dotProduct(d, n);
       if(fabs(denom) < 1e-9) return false;
       double t = dotProduct(a - o, n) / denom;
       
       if(t > 0 && t < t_min) {
            VectorFloatTriplet objPoint = o + d * t;
            
            if (plane.hasTransformation) {
                VectorFloatTriplet worldPoint = transformPoint(*plane.transformMatrix, objPoint);
                VectorFloatTriplet worldNormal = normalize(transformNormal(*plane.normalMatrix, objectNormal));
                double worldDistance = sqrt(dotProduct(worldPoint - ray.origin, worldPoint - ray.origin));
                
                if (worldDistance < t_min) {
                    t_min = worldDistance;
                    intersection.hit = true;
                    intersection.distance = worldDistance;
                    intersection.point = worldPoint;
                    intersection.geometricNormal = worldNormal;
                    intersection.shadingNormal = worldNormal;
                    intersection.material = plane.material;
                    intersection.kind = Intersection::Kind::Plane;
                    intersection.containerIndex = planeIndex;
                    return true;
                }
            } else {
                t_min = t;
                intersection.hit = true;
                intersection.distance = t;
                intersection.point = objPoint;
                intersection.geometricNormal = n;
                intersection.shadingNormal = n;
                intersection.material = plane.material;
                intersection.kind = Intersection::Kind::Plane;
                intersection.containerIndex = planeIndex;
                return true;
            }
       }
       return false;
}

bool rayHitsSphere(
    Ray& ray, 
    const Sphere& sphere, 
    const vector<VectorFloatTriplet>& vertices, 
    double& t_min, 
    Intersection& intersection,
    int sphereIndex
) {
    Ray objectRay = ray;
    if (sphere.hasTransformation) {
        const double* invPtr = sphere.inverseTransformMatrix->m;
        objectRay.origin = transformPointFast(invPtr, ray.origin.x, ray.origin.y, ray.origin.z);
        VectorFloatTriplet dir = transformDirectionFast(invPtr, ray.direction.x, ray.direction.y, ray.direction.z);
        double lenSq = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
        double invLen = 1.0 / sqrtf(lenSq);
        objectRay.direction.x = dir.x * invLen;
        objectRay.direction.y = dir.y * invLen;
        objectRay.direction.z = dir.z * invLen;
    }
    
    VectorFloatTriplet c = vertices[sphere.center];
    VectorFloatTriplet o = objectRay.origin;
    VectorFloatTriplet d = objectRay.direction;
    double r = sphere.radius;
    double D = dotProduct(d, o - c) * dotProduct(d, o - c) 
                - (dotProduct(d, d) * (dotProduct(o - c, o - c) - r * r));
    if(D < 0) return false;
    double t1 = (-dotProduct(d, o - c) + sqrt(D)) / dotProduct(d, d);
    double t2 = (-dotProduct(d, o - c) - sqrt(D)) / dotProduct(d, d);
    double t = min(t1, t2);
    VectorFloatTriplet objectPoint = o + d * t;
    VectorFloatTriplet objectNormal = normalize(objectPoint - c);
    
    if(t > 0 && t < t_min) {
        if (sphere.hasTransformation) {
            const double* transPtr = sphere.transformMatrix->m;
            const double* normPtr = sphere.normalMatrix->m;
            
            VectorFloatTriplet worldPoint = transformPointFast(transPtr, objectPoint.x, objectPoint.y, objectPoint.z);
            VectorFloatTriplet worldNormal = transformNormalFast(normPtr, objectNormal.x, objectNormal.y, objectNormal.z);
            double nlen = sqrtf(worldNormal.x*worldNormal.x + worldNormal.y*worldNormal.y + worldNormal.z*worldNormal.z);
            worldNormal.x /= nlen; worldNormal.y /= nlen; worldNormal.z /= nlen;
            
            double dx = worldPoint.x - ray.origin.x;
            double dy = worldPoint.y - ray.origin.y;
            double dz = worldPoint.z - ray.origin.z;
            double worldDistance = sqrtf(dx*dx + dy*dy + dz*dz);
            
            if (worldDistance < t_min) {
                t_min = worldDistance;
                intersection.hit = true;
                intersection.distance = worldDistance;
                intersection.point = worldPoint;
                intersection.geometricNormal = worldNormal;
                intersection.shadingNormal = worldNormal;
                intersection.material = sphere.material;
                intersection.kind = Intersection::Kind::Sphere;
                intersection.containerIndex = sphereIndex;
                return true;
            }
        } else {
            t_min = t;
            intersection.hit = true;
            intersection.distance = t;
            intersection.point = objectPoint;
            intersection.geometricNormal = objectNormal;
            intersection.shadingNormal = objectNormal;
            intersection.material = sphere.material;
            intersection.kind = Intersection::Kind::Sphere;
            intersection.containerIndex = sphereIndex;
            return true;
        }
    }
    return false;
}


bool rayHitsTriangle(
    Ray& ray, 
    const VectorIntTriplet& face, 
    const vector<VectorFloatTriplet>& vertices, 
    double& t_min, 
    Intersection& intersection, 
    double intersectionTestEpsilon, 
    double determinantT, 
    Material* material,
    bool enableBackFaceCulling,
    int containerIndex,
    int faceIndex
) {
#if PROFILE_PERF
    auto t_tri_start = std::chrono::high_resolution_clock::now();
    g_triangleTests++;
#endif
    const VectorFloatTriplet a = vertices[face.x];
    const VectorFloatTriplet b = vertices[face.y];
    const VectorFloatTriplet c = vertices[face.z];
    
    const VectorFloatTriplet e1 = b - a;
    const VectorFloatTriplet e2 = c - a;
    const VectorFloatTriplet geometricNormal = normalize(crossProduct(e1, e2));
    if (enableBackFaceCulling) {
        if (dotProduct(ray.direction, geometricNormal) > 0.0) {
            return false;
        }
    }
    const double ax=a.x, ay=a.y, az=a.z, bx=b.x, by=b.y, bz=b.z, cx=c.x, cy=c.y, cz=c.z;
    const double ox=ray.origin.x, oy=ray.origin.y, oz=ray.origin.z;
    const double dx=ray.direction.x, dy=ray.direction.y, dz=ray.direction.z;
    
    const double e1x = bx - ax, e1y = by - ay, e1z = bz - az;
    const double e2x = cx - ax, e2y = cy - ay, e2z = cz - az;
    const double rx  = ox - ax, ry  = oy - ay, rz  = oz - az;
 
    double determinant =
        -( e1x * (e2y * dz - e2z * dy)
        - e1y * (e2x * dz - e2z * dx)
        + e1z * (e2x * dy - e2y * dx) );

    // For micro-triangles, determinant can be very small
    // Use a much smaller epsilon than intersectionTestEpsilon
    if (std::fabs(determinant) < 1e-15) {
        return false;
    }
    const double invDet = 1.0 / determinant;

    double determinantBeta =
        -( rx * (e2y * dz - e2z * dy)
        - ry * (e2x * dz - e2z * dx)
        + rz * (e2x * dy - e2y * dx) );
    double beta = determinantBeta * invDet;
    if (beta < 0.0 || beta > 1.0) return false;
 
    double determinantGamma =
        -( e1x * (ry * dz - rz * dy)
        - e1y * (rx * dz - rz * dx)
        + e1z * (rx * dy - ry * dx) );
    double gamma = determinantGamma * invDet;
    if (gamma < 0.0 || gamma > 1.0 || beta + gamma > 1.0) return false;

    // if (std::fabs(determinantT) < 1e-9 || ray.shadowRay || ray.reflectionRay || ray.refractionRay) {
        determinantT =
            e1x * (e2y * rz - e2z * ry)
            - e1y * (e2x * rz - e2z * rx)
            + e1z * (e2x * ry - e2y * rx);
    // }
    double t = determinantT * invDet;
    if (t < intersectionTestEpsilon) return false;

    if(t < t_min) {
        t_min = t;
        intersection.hit = true;
        intersection.distance = t;
        intersection.point = a + beta * (b - a) + gamma * (c - a);
        intersection.geometricNormal = geometricNormal;      // Store geometric normal
        intersection.shadingNormal = geometricNormal;      // Default shading normal (will be updated for smooth shading)
        intersection.beta = beta;
        intersection.gamma = gamma;
        intersection.material = material;
        intersection.containerIndex = containerIndex;
        intersection.faceIndex = faceIndex;
#if PROFILE_PERF
        auto t_tri_end = std::chrono::high_resolution_clock::now();
        g_timeTriangleTest += std::chrono::duration_cast<std::chrono::nanoseconds>(t_tri_end - t_tri_start).count();
#endif
        return true;
    }
 
#if PROFILE_PERF
    auto t_tri_end = std::chrono::high_resolution_clock::now();
    g_timeTriangleTest += std::chrono::duration_cast<std::chrono::nanoseconds>(t_tri_end - t_tri_start).count();
#endif
    return false;
}

bool rayHitsMesh(
    Ray& ray, 
    const Mesh& mesh, 
    const vector<VectorFloatTriplet>& vertices, 
    const vector<double>& determinants, 
    double& t_min, 
    Intersection& intersection,
    double intersectionTestEpsilon,
    scene::MeshBVH* bvh,
    bool enableBackFaceCulling,
    int meshIndex,
    const Matrix4x4* transformMatrix,
    const Matrix4x4* inverseTransformMatrix,
    const Matrix4x4* normalMatrix,
    const Scene* scene,
    Material* materialOverride,
    const scene::AABB* worldSpaceBoundsOverride
) {
#if PROFILE_PERF
    auto t_func_start = std::chrono::high_resolution_clock::now();
#endif
    
    bool hasTransform = (mesh.hasTransformation || transformMatrix != nullptr);
    
#if PROFILE_PERF
    if (hasTransform) g_transformedMeshCalls++;
#endif
    
    // Use material override if provided (for instances), otherwise use mesh's material
    Material* materialToUse = materialOverride ? materialOverride : mesh.material;
    
    if (!hasTransform) {
        if (bvh != nullptr) {
#if PROFILE_PERF
            g_bvhTraversals++;
#endif
            return bvh->traverse(ray, mesh, vertices, determinants, t_min, intersection, intersectionTestEpsilon, enableBackFaceCulling, meshIndex, materialToUse);
        }
        bool hit = false;
        for(int i = 0; i < (int)mesh.faces.size(); i++) {
            if (rayHitsTriangle(ray, mesh.faces[i], vertices, t_min, intersection, intersectionTestEpsilon, determinants[i], materialToUse, enableBackFaceCulling, meshIndex, i)) {
                hit = true;
                intersection.kind = Intersection::Kind::Mesh;
            }
        }
        return hit;
    }
    
    // OPTIMIZATION: Early rejection using world-space bounding box
    // Use override bounds if provided (for instances), otherwise use mesh's bounds
    const scene::AABB* worldBounds = worldSpaceBoundsOverride ? worldSpaceBoundsOverride : mesh.worldSpaceBounds;
    if (worldBounds != nullptr) {
#if PROFILE_PERF
        auto t_aabb_start = std::chrono::high_resolution_clock::now();
#endif
        if (!worldBounds->intersect(ray, 0.0, t_min)) {
#if PROFILE_PERF
            auto t_aabb_end = std::chrono::high_resolution_clock::now();
            g_timeAABBTest += std::chrono::duration_cast<std::chrono::nanoseconds>(t_aabb_end - t_aabb_start).count();
            g_worldBoundsRejects++;
            auto t_func_end = std::chrono::high_resolution_clock::now();
            g_timeRayHitsMesh += std::chrono::duration_cast<std::chrono::nanoseconds>(t_func_end - t_func_start).count();
#endif
            return false;
        }
#if PROFILE_PERF
        auto t_aabb_end = std::chrono::high_resolution_clock::now();
        g_timeAABBTest += std::chrono::duration_cast<std::chrono::nanoseconds>(t_aabb_end - t_aabb_start).count();
        g_worldBoundsAccepts++;
#endif
    }
    
    double original_t_min = t_min;
    
#if PROFILE_PERF
    auto t_transform_start = std::chrono::high_resolution_clock::now();
    g_rayTransforms++;
#endif
    
    // Cache matrix pointers to avoid repeated conditionals
    const Matrix4x4* invTransMat = transformMatrix ? inverseTransformMatrix : mesh.inverseTransformMatrix;
    const double* invTransPtr = invTransMat->m;
    
    // Transform ray to object space using fast inline functions
    Ray objectRay;
    objectRay.origin = transformPointFast(invTransPtr, ray.origin.x, ray.origin.y, ray.origin.z);
    
    VectorFloatTriplet objDir = transformDirectionFast(invTransPtr, ray.direction.x, ray.direction.y, ray.direction.z);
    // Fast normalize
    double lenSq = objDir.x*objDir.x + objDir.y*objDir.y + objDir.z*objDir.z;
    double invLen = 1.0 / sqrtf(lenSq);
    objectRay.direction.x = objDir.x * invLen;
    objectRay.direction.y = objDir.y * invLen;
    objectRay.direction.z = objDir.z * invLen;
    
    objectRay.depth = ray.depth;
    objectRay.shadowRay = ray.shadowRay;
    objectRay.reflectionRay = ray.reflectionRay;
    objectRay.refractionRay = ray.refractionRay;
    
#if PROFILE_PERF
    auto t_transform_end = std::chrono::high_resolution_clock::now();
    g_timeRayTransform += std::chrono::duration_cast<std::chrono::nanoseconds>(t_transform_end - t_transform_start).count();
#endif
    
    bool hit = false;
    
    if (bvh != nullptr) {
#if PROFILE_PERF
        auto t_bvh_start = std::chrono::high_resolution_clock::now();
        g_bvhTraversals++;
#endif
        // For transformed meshes, don't use precomputed determinants (ray origin has changed in local space)
        const vector<double>& dets = hasTransform ? vector<double>() : determinants;
        hit = bvh->traverse(objectRay, mesh, vertices, dets, t_min, intersection, intersectionTestEpsilon, enableBackFaceCulling, meshIndex, materialToUse);
#if PROFILE_PERF
        auto t_bvh_end = std::chrono::high_resolution_clock::now();
        g_timeBVHTraverse += std::chrono::duration_cast<std::chrono::nanoseconds>(t_bvh_end - t_bvh_start).count();
#endif
    } else {
#if PROFILE_PERF
        auto t_tri_start = std::chrono::high_resolution_clock::now();
#endif
        for(int i = 0; i < (int)mesh.faces.size(); i++) {
            // For transformed meshes, always pass 0.0 to force determinant recomputation
            double det = hasTransform ? 0.0 : (i < (int)determinants.size() ? determinants[i] : 0.0);
            if (rayHitsTriangle(objectRay, mesh.faces[i], vertices, t_min, intersection, intersectionTestEpsilon, det, materialToUse, enableBackFaceCulling, meshIndex, i)) {
                hit = true;
                intersection.kind = Intersection::Kind::Mesh;
            }
        }
#if PROFILE_PERF
        auto t_tri_end = std::chrono::high_resolution_clock::now();
        g_timeTriangleTest += std::chrono::duration_cast<std::chrono::nanoseconds>(t_tri_end - t_tri_start).count();
#endif
    }
    
    if (!hit) {
#if PROFILE_PERF
        auto t_func_end = std::chrono::high_resolution_clock::now();
        g_timeRayHitsMesh += std::chrono::duration_cast<std::chrono::nanoseconds>(t_func_end - t_func_start).count();
#endif
        return false;
    }
    if (!hasTransform) {
#if PROFILE_PERF
        auto t_func_end = std::chrono::high_resolution_clock::now();
        g_timeRayHitsMesh += std::chrono::duration_cast<std::chrono::nanoseconds>(t_func_end - t_func_start).count();
#endif
        return true;
    }
    
#if PROFILE_PERF
    auto t_back_start = std::chrono::high_resolution_clock::now();
#endif
    
    // Cache matrix pointers for fast access
    const Matrix4x4* transMat = transformMatrix ? transformMatrix : mesh.transformMatrix;
    const Matrix4x4* normMat = normalMatrix ? normalMatrix : mesh.normalMatrix;
    const double* transPtr = transMat->m;
    const double* normPtr = normMat->m;
    
    // Fast transform intersection point to world space
    VectorFloatTriplet worldPoint = transformPointFast(transPtr, intersection.point.x, intersection.point.y, intersection.point.z);
    
    // Fast distance calculation
    double dx = worldPoint.x - ray.origin.x;
    double dy = worldPoint.y - ray.origin.y;
    double dz = worldPoint.z - ray.origin.z;
    double worldDistance = sqrtf(dx*dx + dy*dy + dz*dz);
    
    if (worldDistance >= original_t_min) {
        t_min = original_t_min;
        return false;
    }

    // Determine if we need to flip normals (for negative scale/reflection)
    bool shouldFlipNormals = transformMatrix ? hasNegativeScale(*transformMatrix) : mesh.hasNegativeScale;
    
    // Fast transform and normalize geometric normal
    VectorFloatTriplet worldGeomNormal = transformNormalFast(normPtr, intersection.geometricNormal.x, 
                                                               intersection.geometricNormal.y, 
                                                               intersection.geometricNormal.z);
    // DON'T flip at all - let the normal matrix handle it!
    // The problem is that ANY manual flipping breaks one case or the other
    double normLenSq = worldGeomNormal.x*worldGeomNormal.x + worldGeomNormal.y*worldGeomNormal.y + worldGeomNormal.z*worldGeomNormal.z;
    double invNormLen = 1.0 / sqrtf(normLenSq);
    worldGeomNormal.x *= invNormLen;
    worldGeomNormal.y *= invNormLen;
    worldGeomNormal.z *= invNormLen;
    
    VectorFloatTriplet worldShadingNormal;
    
    if (mesh.shadingMode == 's' && scene && meshIndex >= 0 && 
        meshIndex < (int)scene->meshVertexNormals.size() &&
        intersection.faceIndex >= 0 && intersection.faceIndex < (int)mesh.faces.size()) {
        
        const auto& face = mesh.faces[intersection.faceIndex];
        if (face.x < (int)scene->meshVertexNormals[meshIndex].size() &&
            face.y < (int)scene->meshVertexNormals[meshIndex].size() &&
            face.z < (int)scene->meshVertexNormals[meshIndex].size()) {
            
            double u = intersection.beta;
            double v = intersection.gamma;
            double w = 1.0 - u - v;
            
            const auto& n0 = scene->meshVertexNormals[meshIndex][face.x];
            const auto& n1 = scene->meshVertexNormals[meshIndex][face.y];
            const auto& n2 = scene->meshVertexNormals[meshIndex][face.z];
            
            // Interpolate and transform in one go
            double nx = w * n0.x + u * n1.x + v * n2.x;
            double ny = w * n0.y + u * n1.y + v * n2.y;
            double nz = w * n0.z + u * n1.z + v * n2.z;
            double nlen = sqrtf(nx*nx + ny*ny + nz*nz);
            double ninv = 1.0 / nlen;
            nx *= ninv; ny *= ninv; nz *= ninv;
            
            worldShadingNormal = transformNormalFast(normPtr, nx, ny, nz);
            double slen = sqrtf(worldShadingNormal.x*worldShadingNormal.x + 
                              worldShadingNormal.y*worldShadingNormal.y + 
                              worldShadingNormal.z*worldShadingNormal.z);
            double sinv = 1.0 / slen;
            worldShadingNormal.x *= sinv;
            worldShadingNormal.y *= sinv;
            worldShadingNormal.z *= sinv;
        } else {
            worldShadingNormal = transformNormalFast(normPtr, intersection.shadingNormal.x, 
                                                      intersection.shadingNormal.y, 
                                                      intersection.shadingNormal.z);
            double slen = sqrtf(worldShadingNormal.x*worldShadingNormal.x + 
                              worldShadingNormal.y*worldShadingNormal.y + 
                              worldShadingNormal.z*worldShadingNormal.z);
            double sinv = 1.0 / slen;
            worldShadingNormal.x *= sinv;
            worldShadingNormal.y *= sinv;
            worldShadingNormal.z *= sinv;
        }
    } else {
        worldShadingNormal = transformNormalFast(normPtr, intersection.shadingNormal.x, 
                                                  intersection.shadingNormal.y, 
                                                  intersection.shadingNormal.z);
        double slen = sqrtf(worldShadingNormal.x*worldShadingNormal.x + 
                          worldShadingNormal.y*worldShadingNormal.y + 
                          worldShadingNormal.z*worldShadingNormal.z);
        double sinv = 1.0 / slen;
        worldShadingNormal.x *= sinv;
        worldShadingNormal.y *= sinv;
        worldShadingNormal.z *= sinv;
    }
    
    // Ensure normals face the camera (critical for mirror reflections and lighting)
    if (dotProduct(worldGeomNormal, ray.direction) > 0.0) {
        worldGeomNormal = -worldGeomNormal;
        worldShadingNormal = -worldShadingNormal;
    }
    
    intersection.point = worldPoint;
    intersection.geometricNormal = worldGeomNormal;
    intersection.shadingNormal = worldShadingNormal;
    intersection.distance = worldDistance;
    t_min = worldDistance;
    intersection.kind = Intersection::Kind::Mesh;
    intersection.containerIndex = -2;
    
#if PROFILE_PERF
    auto t_back_end = std::chrono::high_resolution_clock::now();
    g_timeBackTransform += std::chrono::duration_cast<std::chrono::nanoseconds>(t_back_end - t_back_start).count();
    auto t_func_end = std::chrono::high_resolution_clock::now();
    g_timeRayHitsMesh += std::chrono::duration_cast<std::chrono::nanoseconds>(t_func_end - t_func_start).count();
#endif
    
    return true;
}

Intersection intersect(const Scene& scene, Ray& ray) {
#if PROFILE_PERF
    auto t_start = std::chrono::high_resolution_clock::now();
    g_intersectCalls++;
#endif
    double t_min = numeric_limits<double>::max();
    Intersection intersection;
    bool hit = false;
    
#if PROFILE_PERF
    auto t_planes_start = std::chrono::high_resolution_clock::now();
#endif
    for(int i = 0; i < (int)scene.planes.size(); i++) {
        hit = rayHitsPlane(ray, scene.planes[i], scene.vertices, t_min, intersection, i) || hit;
    }
#if PROFILE_PERF
    auto t_planes_end = std::chrono::high_resolution_clock::now();
    g_timeIntersectPlanes += std::chrono::duration_cast<std::chrono::nanoseconds>(t_planes_end - t_planes_start).count();
    
    auto t_tris_start = std::chrono::high_resolution_clock::now();
#endif
    for(int i = 0; i < (int)scene.triangles.size(); i++) {
        const Triangle& tri = scene.triangles[i];
        
        if (tri.hasTransformation) {
            Ray objectRay;
            objectRay.origin = transformPoint(*tri.inverseTransformMatrix, ray.origin);
            objectRay.direction = normalize(transformDirection(*tri.inverseTransformMatrix, ray.direction));
            objectRay.depth = ray.depth;
            objectRay.shadowRay = ray.shadowRay;
            objectRay.reflectionRay = ray.reflectionRay;
            objectRay.refractionRay = ray.refractionRay;
            
            double original_t_min = t_min;
            if (rayHitsTriangle(objectRay, tri.indices, scene.vertices, t_min, intersection, scene.intersectionTestEpsilon, 0.0, tri.material, scene.enableBackFaceCulling, -1, i)) {
                VectorFloatTriplet worldPoint = transformPoint(*tri.transformMatrix, intersection.point);
                VectorFloatTriplet worldGeomNormal = normalize(transformNormal(*tri.normalMatrix, intersection.geometricNormal));
                VectorFloatTriplet worldShadingNormal = normalize(transformNormal(*tri.normalMatrix, intersection.shadingNormal));
                double worldDistance = sqrt(dotProduct(worldPoint - ray.origin, worldPoint - ray.origin));
                
                if (worldDistance < original_t_min) {
                    intersection.point = worldPoint;
                    intersection.geometricNormal = worldGeomNormal;
                    intersection.shadingNormal = worldShadingNormal;
                    intersection.distance = worldDistance;
                    t_min = worldDistance;
                    hit = true;
                    intersection.kind = Intersection::Kind::Triangle;
                } else {
                    t_min = original_t_min;
                }
            }
        } else {
            if (rayHitsTriangle(ray, tri.indices, scene.vertices, t_min, intersection, scene.intersectionTestEpsilon, scene.cameraTriangleDeterminant[scene.currentCameraIndex][i], tri.material, scene.enableBackFaceCulling, -1, i)) {
                hit = true;
                intersection.kind = Intersection::Kind::Triangle;
            }
        }
    }
#if PROFILE_PERF
    auto t_tris_end = std::chrono::high_resolution_clock::now();
    g_timeIntersectTriangles += std::chrono::duration_cast<std::chrono::nanoseconds>(t_tris_end - t_tris_start).count();
    
    auto t_meshes_start = std::chrono::high_resolution_clock::now();
#endif
    for(int i = 0; i < (int)scene.meshes.size(); i++) {
        MeshBVH* bvh = (i < scene.meshBVHs.size()) ? scene.meshBVHs[i] : nullptr;
        hit = rayHitsMesh(ray, scene.meshes[i], scene.vertices, scene.cameraMeshDeterminant[scene.currentCameraIndex][i], t_min, intersection, scene.intersectionTestEpsilon, bvh, scene.enableBackFaceCulling, i, nullptr, nullptr, nullptr, &scene) || hit;
    }
#if PROFILE_PERF
    auto t_meshes_end = std::chrono::high_resolution_clock::now();
    g_timeIntersectMeshes += std::chrono::duration_cast<std::chrono::nanoseconds>(t_meshes_end - t_meshes_start).count();
    
    auto t_spheres_start = std::chrono::high_resolution_clock::now();
#endif
    for(int i = 0; i < (int)scene.spheres.size(); i++) {
        hit = rayHitsSphere(ray, scene.spheres[i], scene.vertices, t_min, intersection, i) || hit;
    }
#if PROFILE_PERF
    auto t_spheres_end = std::chrono::high_resolution_clock::now();
    g_timeIntersectSpheres += std::chrono::duration_cast<std::chrono::nanoseconds>(t_spheres_end - t_spheres_start).count();
    
    auto t_instances_start = std::chrono::high_resolution_clock::now();
#endif
    
    for(int i = 0; i < (int)scene.meshInstances.size(); i++) {
        const MeshInstance& instance = scene.meshInstances[i];
        if (instance.baseMesh && instance.baseMeshIndex >= 0) {
            // OPTIMIZATION: Early rejection using world-space bounding box
            if (instance.worldSpaceBounds != nullptr) {
                if (!instance.worldSpaceBounds->intersect(ray, 0.0, t_min)) {
#if PROFILE_PERF
                    g_worldBoundsRejects++;
#endif
                    continue;
                }
#if PROFILE_PERF
                g_worldBoundsAccepts++;
#endif
            }
            
            int baseMeshIndex = instance.baseMeshIndex;
            
            // CRITICAL FIX: Don't copy the mesh! Just use the base mesh directly with different matrices
            // Pass instance material AND world-space bounds as overrides (thread-safe)
            MeshBVH* bvh = (baseMeshIndex < scene.meshBVHs.size()) ? scene.meshBVHs[baseMeshIndex] : nullptr;
            
            bool instanceHit = rayHitsMesh(ray, *instance.baseMesh, scene.vertices, 
                                         scene.cameraMeshDeterminant[scene.currentCameraIndex][baseMeshIndex], 
                                         t_min, intersection, scene.intersectionTestEpsilon, 
                                         bvh, scene.enableBackFaceCulling, baseMeshIndex,
                                         instance.transformMatrix, instance.inverseTransformMatrix, instance.normalMatrix, &scene,
                                         instance.material,  // Pass instance material as override
                                         instance.worldSpaceBounds);  // Pass instance's world-space bounds
            
            if (instanceHit) {
                hit = true;
            }
        }
    }
#if PROFILE_PERF
    auto t_instances_end = std::chrono::high_resolution_clock::now();
    g_timeIntersectInstances += std::chrono::duration_cast<std::chrono::nanoseconds>(t_instances_end - t_instances_start).count();
#endif
    
    if (!hit) {
#if PROFILE_PERF
        auto t_end = std::chrono::high_resolution_clock::now();
        g_timeIntersect += std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count();
#endif
        return intersection;
    }
    
#if PROFILE_PERF
    auto t_postprocess_start = std::chrono::high_resolution_clock::now();
#endif
    const double u = intersection.beta;
    const double v = intersection.gamma;
    const double w = 1.0 - u - v;

    intersection.geometricNormal = normalize(intersection.geometricNormal);
    intersection.shadingNormal = normalize(intersection.shadingNormal);
    
    if (intersection.kind == Intersection::Kind::Mesh && intersection.containerIndex != -2) {
        const int meshIdx = intersection.containerIndex;
        
        if (meshIdx >= 0 && meshIdx < (int)scene.meshes.size()) {
            const Mesh& mesh = scene.meshes[meshIdx];
            
            if (mesh.shadingMode == 's' && 
                meshIdx < (int)scene.meshVertexNormals.size() &&
                intersection.faceIndex >= 0 && 
                intersection.faceIndex < (int)mesh.faces.size()) {
                
                const auto& face = mesh.faces[intersection.faceIndex];
                if (face.x < (int)scene.meshVertexNormals[meshIdx].size() &&
                    face.y < (int)scene.meshVertexNormals[meshIdx].size() &&
                    face.z < (int)scene.meshVertexNormals[meshIdx].size()) {
                    
                const auto& n0 = scene.meshVertexNormals[meshIdx][face.x];
                const auto& n1 = scene.meshVertexNormals[meshIdx][face.y];
                const auto& n2 = scene.meshVertexNormals[meshIdx][face.z];
                VectorFloatTriplet smoothNormal = normalize(w * n0 + u * n1 + v * n2);
                
                if (mesh.hasTransformation && mesh.normalMatrix) {
                    smoothNormal = normalize(transformNormal(*mesh.normalMatrix, smoothNormal));
                }
                
                intersection.shadingNormal = smoothNormal;
                
                // Ensure both normals face the camera
                // This is critical for mesh instances with negative scale transformations
                if (dotProduct(intersection.geometricNormal, ray.direction) > 0.0) {
                    intersection.geometricNormal = -intersection.geometricNormal;
                    intersection.shadingNormal = -intersection.shadingNormal;
                } else if (dotProduct(intersection.shadingNormal, intersection.geometricNormal) < 0.0) {
                    // Align shading normal with geometric normal
                    intersection.shadingNormal = -intersection.shadingNormal;
                }
                }
            }
        }
    }
#if PROFILE_PERF
    auto t_postprocess_end = std::chrono::high_resolution_clock::now();
    g_timeIntersectPostProcess += std::chrono::duration_cast<std::chrono::nanoseconds>(t_postprocess_end - t_postprocess_start).count();
#endif
    
#if PROFILE_PERF
    auto t_end = std::chrono::high_resolution_clock::now();
    g_timeIntersect += std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count();
#endif
    
    return intersection;
}


VectorFloatTriplet computeShading(const Scene& scene, Ray& ray, const Intersection& intersection) {
    Material* material = intersection.material;
    VectorFloatTriplet color{0.0, 0.0, 0.0};

    if (!material) {
        return color;
    }

    // Use shading normal for lighting
    VectorFloatTriplet normal = normalize(intersection.shadingNormal);

    // --- Ambient ---
    color += material->ambientReflectance * scene.ambientLight.intensity;

    // ==============================
    //     CONDUCTOR (METAL)
    // ==============================
    if (material->type == "conductor") {
        // Direct specular from point lights (Blinn–Phong, tinted by metal color)
        for (const PointLight& light : scene.pointLights) {
            if (isInShadow(scene, ray, light, intersection)) {
                continue;
            }

            VectorFloatTriplet lightVec = light.position - intersection.point;
            double distanceSq = dotProduct(lightVec, lightVec);
            double distance = std::sqrt(distanceSq);
            VectorFloatTriplet lightDir = lightVec * (1.0 / distance);

            double attenuation = 1.0 / distanceSq;

            VectorFloatTriplet viewDir = normalize(-ray.direction);
            VectorFloatTriplet halfVector = normalize(lightDir + viewDir);

            double specularFactor = std::pow(
                std::max(0.0, dotProduct(normal, halfVector)),
                material->phongExponent
            );

            // mirrorReflectance is the metal color (RGB)
            VectorFloatTriplet specular =
                material->mirrorReflectance * light.intensity *
                specularFactor * attenuation;

            color += specular;
        }

        // Environment / scene reflection with Fresnel
        double cosTheta = dotProduct(-ray.direction, normal);
        cosTheta = std::max(0.0, std::min(1.0, cosTheta));

        double F = fresnelConductor(cosTheta,
                                    material->refractionIndex,
                                    material->absorptionIndex);

        Ray reflectedRay = reflect(ray, normal, intersection.point, scene.shadowRayEpsilon);
        Intersection reflectionIntersection = intersect(scene, reflectedRay);
        VectorFloatTriplet reflectedColor =
            computePixelColor(scene, reflectedRay, reflectionIntersection);

        color += F * material->mirrorReflectance * reflectedColor;
        return color;
    }

    // ==============================
    //     DIELECTRIC (GLASS, etc.)
    // ==============================
    if (material->type == "dielectric") {
        double n1 = 1.0;
        double n2 = material->refractionIndex;
        bool entering = dotProduct(ray.direction, normal) < 0.0;

        // Flip normal if we’re exiting the medium
        if (!entering) {
            normal = -normal;
            std::swap(n1, n2);
        }

        double cosTheta = dotProduct(-ray.direction, normal);
        cosTheta = std::max(-1.0, std::min(1.0, cosTheta));

        double F = fresnelDielectric(cosTheta, n1, n2);

        // Reflection
        Ray reflectedRay = reflect(ray, normal, intersection.point, scene.shadowRayEpsilon);
        Intersection reflectionIntersection = intersect(scene, reflectedRay);
        VectorFloatTriplet reflectedColor =
            computePixelColor(scene, reflectedRay, reflectionIntersection);

        // Refraction
        bool totalInternalReflection = false;
        Ray refractedRay = refract(ray, normal, n1, n2,
                                   intersection.point,
                                   scene.shadowRayEpsilon,
                                   totalInternalReflection);

        if (totalInternalReflection) {
            // Pure reflection if we can’t refract
            color += reflectedColor;
            return color;
        }

        Intersection refractionIntersection = intersect(scene, refractedRay);
        VectorFloatTriplet refractedColor =
            computePixelColor(scene, refractedRay, refractionIntersection);

        // Beer–Lambert absorption when travelling inside the medium
        if (entering && refractionIntersection.hit) {
            double distance = refractionIntersection.distance;
            VectorFloatTriplet absorbance = material->absorptionCoefficient * distance;
            VectorFloatTriplet transmittance{
                std::exp(-absorbance.x),
                std::exp(-absorbance.y),
                std::exp(-absorbance.z)
            };
            refractedColor = refractedColor * transmittance;
        }

        color += F * reflectedColor + (1.0 - F) * refractedColor;
        return color;
    }

    // ==============================
    //   OPAQUE / PHONG MATERIALS
    // ==============================

    // Optional perfect mirror term on top
    if (material->isMirror) {
        // Use shading normal for reflection to ensure consistency with lighting
        Ray reflectedRay = reflect(ray,
                                   intersection.shadingNormal,
                                   intersection.point,
                                   scene.shadowRayEpsilon);
        Intersection reflectionIntersection = intersect(scene, reflectedRay);
        VectorFloatTriplet mirrorColor =
            computePixelColor(scene, reflectedRay, reflectionIntersection);
        color += material->mirrorReflectance * mirrorColor;
    }

    // Direct illumination from point lights (diffuse + specular)
    for (const PointLight& light : scene.pointLights) {
        if (isInShadow(scene, ray, light, intersection)) {
            continue;
        }

        VectorFloatTriplet lightVec = light.position - intersection.point;
        double distanceSq = dotProduct(lightVec, lightVec);
        double distance = std::sqrt(distanceSq);
        VectorFloatTriplet lightDir = lightVec * (1.0 / distance);

        double attenuation = 1.0 / distanceSq;

        // Diffuse
        double diffuseFactor = std::max(0.0, dotProduct(normal, lightDir));
        VectorFloatTriplet diffuse =
            material->diffuseReflectance * light.intensity *
            diffuseFactor * attenuation;

        // Blinn–Phong specular
        VectorFloatTriplet viewDir = normalize(-ray.direction);
        VectorFloatTriplet halfVector = normalize(lightDir + viewDir);
        double specularFactor = std::pow(
            std::max(0.0, dotProduct(normal, halfVector)),
            material->phongExponent
        );
        VectorFloatTriplet specular =
            material->specularReflectance * light.intensity *
            specularFactor * attenuation;

        color += diffuse + specular;
    }

    return color;
}

// =======================================================
//    RECURSION / PIXEL COLOR
// =======================================================
VectorFloatTriplet computePixelColor(const Scene& scene, Ray& ray, const Intersection& intersection) {
    if (ray.depth > scene.maxRecursionDepth) {
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }

    if (intersection.hit) {
        return computeShading(scene, ray, intersection);
    }

    // Only primary ray sees background color
    if (ray.depth == 0) {
        return scene.backgroundColor;
    }

    return VectorFloatTriplet{0.0, 0.0, 0.0};
}

// =======================================================
//    SHADOW TEST
// =======================================================
bool isInShadow(const Scene& scene, Ray& /*ray*/, const PointLight& light, const Intersection& intersection) {
    Ray shadowRay{
        intersection.point + scene.shadowRayEpsilon * intersection.geometricNormal,
        normalize(light.position - intersection.point),
        0,
        true,   // shadow ray
        false,  // reflection ray
        false   // refraction ray
    };

    Intersection shadowIntersection = intersect(scene, shadowRay);
    double distToLight = std::sqrt(
        dotProduct(light.position - intersection.point,
                   light.position - intersection.point)
    );

    // Anything between the point and the light blocks it (simple hard shadows)
    return shadowIntersection.hit && shadowIntersection.distance < distToLight;
}

// =======================================================
//    FRESNEL FUNCTIONS
// =======================================================
double fresnelConductor(double cosTheta, double n, double k) {
    cosTheta = std::max(0.0, std::min(1.0, cosTheta));
    double cosThetaSq = cosTheta * cosTheta;

    double n2 = n * n;
    double k2 = k * k;
    double n2PlusK2 = n2 + k2;

    double twoNCos = 2.0 * n * cosTheta;

    double RsNum = (n2PlusK2 - twoNCos + cosThetaSq);
    double RsDen = (n2PlusK2 + twoNCos + cosThetaSq);
    double Rs = RsNum / RsDen;

    double RpNum = (n2PlusK2 * cosThetaSq - twoNCos + 1.0);
    double RpDen = (n2PlusK2 * cosThetaSq + twoNCos + 1.0);
    double Rp = RpNum / RpDen;

    double R = 0.5 * (Rs + Rp);
    if (R < 0.0) R = 0.0;
    if (R > 1.0) R = 1.0;
    return R;
}

double fresnelDielectric(double cosTheta, double n1, double n2) {
    cosTheta = std::max(-1.0, std::min(1.0, cosTheta));
    double eta = n1 / n2;
    double sinThetaTSq = eta * eta * (1.0 - cosTheta * cosTheta);

    if (sinThetaTSq >= 1.0) {
        // Total internal reflection
        return 1.0;
    }

    double cosPhi = std::sqrt(std::max(0.0, 1.0 - sinThetaTSq));

    double rPerp = (n1 * cosTheta - n2 * cosPhi) / (n1 * cosTheta + n2 * cosPhi);
    double rPara = (n2 * cosTheta - n1 * cosPhi) / (n2 * cosTheta + n1 * cosPhi);

    return 0.5 * (rPerp * rPerp + rPara * rPara);
}

// =======================================================
//    REFLECT / REFRACT RAYS
// =======================================================
Ray refract(Ray& ray,
            const VectorFloatTriplet normal,
            double n1,
            double n2,
            VectorFloatTriplet point,
            const double shadowRayEpsilon,
            bool& totalInternalReflection) {
    double eta = n1 / n2;
    double cosTheta = -dotProduct(ray.direction, normal);
    cosTheta = std::max(-1.0, std::min(1.0, cosTheta));

    double sinThetaTSq = eta * eta * (1.0 - cosTheta * cosTheta);
    totalInternalReflection = (sinThetaTSq > 1.0);

    if (totalInternalReflection) {
        // Return a dummy ray; caller will ignore refracted path
        return Ray{point, ray.direction, ray.depth, false, false, false};
    }

    double cosPhi = std::sqrt(std::max(0.0, 1.0 - sinThetaTSq));
    VectorFloatTriplet wt =
        (ray.direction + normal * cosTheta) * eta - normal * cosPhi;

    return Ray{
        point - shadowRayEpsilon * normalize(normal),
        normalize(wt),
        ray.depth + 1,
        false,
        false,
        true
    };
}

Ray reflect(Ray& ray,
            const VectorFloatTriplet normal,
            VectorFloatTriplet point,
            const double shadowRayEpsilon) {
    VectorFloatTriplet n = normalize(normal);
    VectorFloatTriplet dir = ray.direction - 2.0 * dotProduct(ray.direction, n) * n;

    return Ray{
        point + shadowRayEpsilon * n,
        normalize(dir),
        ray.depth + 1,
        false,
        true,
        false
    };
}

Matrix4x4 identityMatrix() {
    Matrix4x4 m;
    for (int i = 0; i < 16; i++) m.m[i] = 0.0;
    m.m[0] = m.m[5] = m.m[10] = m.m[15] = 1.0;
    return m;
}

Matrix4x4 buildTranslationMatrix(const Translation& t) {
    Matrix4x4 m = identityMatrix();
    m.m[3] = t.data.x;
    m.m[7] = t.data.y;
    m.m[11] = t.data.z;
    return m;
}

Matrix4x4 buildScalingMatrix(const Scaling& s) {
    Matrix4x4 m = identityMatrix();
    m.m[0] = s.data.x;
    m.m[5] = s.data.y;
    m.m[10] = s.data.z;
    return m;
}

Matrix4x4 buildRotationMatrix(const Rotation& r) {
    double angle = r.angle * M_PI / 180.0;
    double c = cos(angle);
    double s = sin(angle);
    VectorFloatTriplet u = normalize(r.axis);
    
    Matrix4x4 m = identityMatrix();
    m.m[0] = u.x*u.x*(1-c) + c;
    m.m[1] = u.x*u.y*(1-c) - u.z*s;
    m.m[2] = u.x*u.z*(1-c) + u.y*s;
    
    m.m[4] = u.y*u.x*(1-c) + u.z*s;
    m.m[5] = u.y*u.y*(1-c) + c;
    m.m[6] = u.y*u.z*(1-c) - u.x*s;
    
    m.m[8] = u.z*u.x*(1-c) - u.y*s;
    m.m[9] = u.z*u.y*(1-c) + u.x*s;
    m.m[10] = u.z*u.z*(1-c) + c;
    
    return m;
}

Matrix4x4 buildCompositeMatrix(const Composite& c) {
    Matrix4x4 m;
    for (int i = 0; i < 16; i++) {
        m.m[i] = c.data[i];
    }
    return m;
}

Matrix4x4 multiplyMatrices(const Matrix4x4& a, const Matrix4x4& b) {
    Matrix4x4 result;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result.m[i*4 + j] = 0.0;
            for (int k = 0; k < 4; k++) {
                result.m[i*4 + j] += a.m[i*4 + k] * b.m[k*4 + j];
            }
        }
    }
    return result;
}

Matrix4x4 transposeMatrix(const Matrix4x4& m) {
    Matrix4x4 result;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result.m[i*4 + j] = m.m[j*4 + i];
        }
    }
    return result;
}

Matrix4x4 invertMatrix(const Matrix4x4& mat) {
    const double* m = mat.m;
    Matrix4x4 inv;
    double* invOut = inv.m;
    
    invOut[0] = m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15] + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
    invOut[4] = -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15] - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
    invOut[8] = m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15] + m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
    invOut[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14] - m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
    invOut[1] = -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15] - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
    invOut[5] = m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15] + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
    invOut[9] = -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15] - m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
    invOut[13] = m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14] + m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
    invOut[2] = m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15] + m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
    invOut[6] = -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15] - m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
    invOut[10] = m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15] + m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
    invOut[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14] - m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
    invOut[3] = -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11] - m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
    invOut[7] = m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11] + m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
    invOut[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11] - m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
    invOut[15] = m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10] + m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];
    
    double det = m[0]*invOut[0] + m[1]*invOut[4] + m[2]*invOut[8] + m[3]*invOut[12];
    
    if (fabs(det) < 1e-9) {
        return identityMatrix();
    }
    
    det = 1.0 / det;
    for (int i = 0; i < 16; i++) {
        invOut[i] *= det;
    }
    
    return inv;
}

VectorFloatTriplet transformPoint(const Matrix4x4& m, const VectorFloatTriplet& p) {
    return transformPointFast(m.m, p.x, p.y, p.z);
}

VectorFloatTriplet transformDirection(const Matrix4x4& m, const VectorFloatTriplet& d) {
    return transformDirectionFast(m.m, d.x, d.y, d.z);
}

VectorFloatTriplet transformNormal(const Matrix4x4& invTranspose, const VectorFloatTriplet& n) {
    return transformNormalFast(invTranspose.m, n.x, n.y, n.z);
}

Matrix4x4 buildObjectTransformMatrix(const Scene& scene, const std::vector<TransformationRef>& transforms) {
    Matrix4x4 result = identityMatrix();
    
    for (int i = transforms.size() - 1; i >= 0; i--) {
        const TransformationRef& ref = transforms[i];
        Matrix4x4 t;
        
        if (ref.type == 't') {
            auto it = scene.translationIdToIndex.find(ref.id);
            if (it != scene.translationIdToIndex.end()) {
                t = buildTranslationMatrix(scene.translations[it->second]);
            } else {
                continue;
            }
        } else if (ref.type == 's') {
            auto it = scene.scalingIdToIndex.find(ref.id);
            if (it != scene.scalingIdToIndex.end()) {
                t = buildScalingMatrix(scene.scalings[it->second]);
            } else {
                continue;
            }
        } else if (ref.type == 'r') {
            auto it = scene.rotationIdToIndex.find(ref.id);
            if (it != scene.rotationIdToIndex.end()) {
                t = buildRotationMatrix(scene.rotations[it->second]);
            } else {
                continue;
            }
        } else if (ref.type == 'c') {
            auto it = scene.compositeIdToIndex.find(ref.id);
            if (it != scene.compositeIdToIndex.end()) {
                t = buildCompositeMatrix(scene.composites[it->second]);
            } else {
                continue;
            }
        } else {
            continue;
        }
        
        result = multiplyMatrices(result, t);
    }
    
    return result;
}

bool hasNegativeScale(const Matrix4x4& m) {
    // Compute determinant of the 3x3 upper-left submatrix
    // If determinant is negative, the transformation includes a reflection
    const double* mat = m.m;
    double det = mat[0] * (mat[5] * mat[10] - mat[6] * mat[9])
              - mat[1] * (mat[4] * mat[10] - mat[6] * mat[8])
              + mat[2] * (mat[4] * mat[9] - mat[5] * mat[8]);
    return det < 0.0;
}