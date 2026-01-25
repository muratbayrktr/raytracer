#include <math.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <chrono>
#include <sstream>
#include <map>
#include <fstream>
#include "scene.h"
#include "utils.h"
#include "overloads.h"
#include "bvh.h"
#include "precompute.h"
#include "texture.h"
#include "brdf.h"

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

static inline Ray applyMotionBlurToRay(const Ray& ray, const VectorFloatTriplet& motionBlur) {
    Ray offsetRay = ray;
    offsetRay.origin.x = ray.origin.x - ray.time * motionBlur.x;
    offsetRay.origin.y = ray.origin.y - ray.time * motionBlur.y;
    offsetRay.origin.z = ray.origin.z - ray.time * motionBlur.z;
    return offsetRay;
}

static inline void correctHitPointForMotionBlur(Intersection& intersection, const VectorFloatTriplet& motionBlur, double time) {
    intersection.point.x += time * motionBlur.x;
    intersection.point.y += time * motionBlur.y;
    intersection.point.z += time * motionBlur.z;
}

static inline double recomputeDistanceFromOrigin(const VectorFloatTriplet& point, const VectorFloatTriplet& origin) {
    double dx = point.x - origin.x;
    double dy = point.y - origin.y;
    double dz = point.z - origin.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

static inline void finalizeMotionBlurHit(const Ray& originalRay,
                                         const VectorFloatTriplet& motionBlur,
                                         double time,
                                         Intersection& intersection,
                                         double& t_min) {
    correctHitPointForMotionBlur(intersection, motionBlur, time);
    double newDistance = recomputeDistanceFromOrigin(intersection.point, originalRay.origin);
    intersection.distance = newDistance;
    t_min = newDistance;
}

void clamp(VectorFloatTriplet& color, int min, int max) {
    if (color.x < min) color.x = min;
    if (color.x > max) color.x = max;
    if (color.y < min) color.y = min;
    if (color.y > max) color.y = max;
    if (color.z < min) color.z = min;
    if (color.z > max) color.z = max;
}

Ray castRay(const Camera& camera,
            double x,
            double y,
            int width,
            int height,
            double time,
            double random1,
            double random2) {
    // Safety check: ensure gaze is not zero
    VectorFloatTriplet gaze = camera.gaze;
    double gazeLen = std::sqrt(dotProduct(gaze, gaze));
    if (gazeLen < 1e-10) {
        gaze = VectorFloatTriplet{0, 0, -1};  // Default gaze
    } else {
        gaze = gaze * (1.0 / gazeLen);
    }
    
    VectorFloatTriplet w = -gaze;
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
    
    VectorFloatTriplet ray_direction = s - e;
    VectorFloatTriplet origin = e;
    
    // Depth-of-field: if aperture > 0, sample lens and aim at focal point
    if (camera.apertureSize > 0.0) {
        // Focal point: scale ray to reach focal plane at focusDistance
        double scale = camera.focusDistance / camera.nearDistance;
        VectorFloatTriplet focalPoint = e + ray_direction * scale;
        
        double halfAperture = camera.apertureSize / 2.0;
        double lensU = (random1 * 2.0 - 1.0) * halfAperture;
        double lensV = (random2 * 2.0 - 1.0) * halfAperture;
        origin = e + u * lensU + v * lensV;
        
        ray_direction = focalPoint - origin;
    }
    
    Ray ray = Ray(origin,
                  normalize(ray_direction),
                  /*depth*/ 0,
                  /*shadow*/ false,
                  /*reflection*/ false,
                  /*refraction*/ false,
                  time,
                  random1,
                  random2);
    return ray;
}

bool rayHitsPlane(
    Ray& ray, 
    const Plane& plane, 
    const vector<VectorFloatTriplet>& vertices, 
    double& t_min, 
    Intersection& intersection,
    int planeIndex,
    double minDistance
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
                
                if (worldDistance < minDistance) {
                    return false;
                }
                
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
                if (t < minDistance) {
                    return false;
                }
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
    int sphereIndex,
    double minDistance
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
    
    // Pick the closest positive t value
    // t2 <= t1 always (since t2 uses -sqrt)
    // If t2 > 0, use it (ray hits front surface from outside)
    // If t2 <= 0 but t1 > 0, use t1 (ray inside sphere, hits exit point)
    double t;
    if (t2 > 0) {
        t = t2;  // Front surface hit (entering sphere from outside)
    } else if (t1 > 0) {
        t = t1;  // Back surface hit (exiting sphere from inside)
    } else {
        return false;  // Both negative, sphere is behind the ray
    }
    
    VectorFloatTriplet objectPoint = o + d * t;
    VectorFloatTriplet objectNormal = normalize(objectPoint - c);
    
    if(t < t_min) {
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
            
            // Near-plane clipping for primary camera rays
            if (worldDistance < minDistance) {
                return false;
            }
            
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
            // For untransformed spheres, t is distance along (normalized) ray direction
            if (t < minDistance) {
                return false;
            }
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
    int faceIndex,
    double minDistance
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
    if (t < intersectionTestEpsilon) {
#if PROFILE_PERF
        auto t_tri_end = std::chrono::high_resolution_clock::now();
        g_timeTriangleTest += std::chrono::duration_cast<std::chrono::nanoseconds>(t_tri_end - t_tri_start).count();
#endif
        return false;
    }

    // Near-plane clipping for primary camera rays (in the current ray's metric)
    if (t < minDistance) {
#if PROFILE_PERF
        auto t_tri_end = std::chrono::high_resolution_clock::now();
        g_timeTriangleTest += std::chrono::duration_cast<std::chrono::nanoseconds>(t_tri_end - t_tri_start).count();
#endif
        return false;
    }

    if (t < t_min) {
        t_min = t;
        intersection.hit = true;
        intersection.distance = t;
        intersection.point = a + beta * (b - a) + gamma * (c - a);
        intersection.geometricNormal = geometricNormal;   // geometric normal
        intersection.shadingNormal = geometricNormal;     // default shading normal
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
    const scene::AABB* worldSpaceBoundsOverride,
    double minDistance
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

    // ---------------------------
    // UNTRANSFORMED MESH PATH
    // ---------------------------
    // For untransformed meshes, local t IS world distance (ray direction is normalized).
    // We use t_min directly for comparisons - the caller handles motion blur correction.
    if (!hasTransform) {
        bool hit = false;

        // BVH-accelerated path
        if (bvh != nullptr) {
#if PROFILE_PERF
            g_bvhTraversals++;
#endif
            hit = bvh->traverse(
                ray,
                mesh,
                vertices,
                determinants,
                t_min,
                intersection,
                intersectionTestEpsilon,
                enableBackFaceCulling,
                meshIndex,
                materialToUse,
                minDistance
            );

            if (hit) {
                intersection.kind = Intersection::Kind::Mesh;
            }

#if PROFILE_PERF
            auto t_func_end = std::chrono::high_resolution_clock::now();
            g_timeRayHitsMesh += std::chrono::duration_cast<std::chrono::nanoseconds>(t_func_end - t_func_start).count();
#endif
            return hit;
        }

        // Brute-force triangle path
        for (int i = 0; i < (int)mesh.faces.size(); i++) {
            if (rayHitsTriangle(ray,
                                mesh.faces[i],
                                vertices,
                                t_min,
                                intersection,
                                intersectionTestEpsilon,
                                determinants[i],
                                materialToUse,
                                enableBackFaceCulling,
                                meshIndex,
                                i,
                                minDistance)) {
                hit = true;
                intersection.kind = Intersection::Kind::Mesh;
            }
        }

#if PROFILE_PERF
        auto t_func_end = std::chrono::high_resolution_clock::now();
        g_timeRayHitsMesh += std::chrono::duration_cast<std::chrono::nanoseconds>(t_func_end - t_func_start).count();
#endif
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
    objectRay.time = ray.time;
    objectRay.random1 = ray.random1;
    objectRay.random2 = ray.random2;
    
#if PROFILE_PERF
    auto t_transform_end = std::chrono::high_resolution_clock::now();
    g_timeRayTransform += std::chrono::duration_cast<std::chrono::nanoseconds>(t_transform_end - t_transform_start).count();
#endif
    
    bool hit = false;

    // Local best hit for this mesh in OBJECT SPACE
    double local_t_min = std::numeric_limits<double>::max();
    Intersection localIntersection;

    if (bvh != nullptr) {
#if PROFILE_PERF
        auto t_bvh_start = std::chrono::high_resolution_clock::now();
        g_bvhTraversals++;
#endif
        // For transformed meshes, don't use precomputed determinants (ray origin has changed in local space)
        const vector<double>& dets = hasTransform ? vector<double>() : determinants;
        // For transformed meshes, near-plane clipping is applied in world space after back-transform,
        // so we pass minDistance = 0.0 here.
        hit = bvh->traverse(objectRay,
                            mesh,
                            vertices,
                            dets,
                            local_t_min,
                            localIntersection,
                            intersectionTestEpsilon,
                            enableBackFaceCulling,
                            meshIndex,
                            materialToUse,
                            0.0);
#if PROFILE_PERF
        auto t_bvh_end = std::chrono::high_resolution_clock::now();
        g_timeBVHTraverse += std::chrono::duration_cast<std::chrono::nanoseconds>(t_bvh_end - t_bvh_start).count();
#endif
    } else {
#if PROFILE_PERF
        auto t_tri_start = std::chrono::high_resolution_clock::now();
#endif
        for (int i = 0; i < (int)mesh.faces.size(); i++) {
            // For transformed meshes, always pass 0.0 to force determinant recomputation
            double det = hasTransform ? 0.0 : (i < (int)determinants.size() ? determinants[i] : 0.0);
            if (rayHitsTriangle(objectRay,
                                mesh.faces[i],
                                vertices,
                                local_t_min,
                                localIntersection,
                                intersectionTestEpsilon,
                                det,
                                materialToUse,
                                enableBackFaceCulling,
                                meshIndex,
                                i)) {
                hit = true;
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
    
#if PROFILE_PERF
    auto t_back_start = std::chrono::high_resolution_clock::now();
#endif
    
    // Cache matrix pointers for fast access
    const Matrix4x4* transMat = transformMatrix ? transformMatrix : mesh.transformMatrix;
    const Matrix4x4* normMat = normalMatrix ? normalMatrix : mesh.normalMatrix;
    const double* transPtr = transMat->m;
    const double* normPtr = normMat->m;
    
    // Fast transform intersection point to world space (from localIntersection)
    VectorFloatTriplet worldPoint = transformPointFast(
        transPtr,
        localIntersection.point.x,
        localIntersection.point.y,
        localIntersection.point.z);
    
    // Compute world-space distance from the ORIGINAL ray origin
    double worldDistance = recomputeDistanceFromOrigin(worldPoint, ray.origin);
    
    // Near-plane clipping for primary camera rays in world space
    if (worldDistance < minDistance) {
#if PROFILE_PERF
        auto t_back_end = std::chrono::high_resolution_clock::now();
        g_timeBackTransform += std::chrono::duration_cast<std::chrono::nanoseconds>(t_back_end - t_back_start).count();
        auto t_func_end = std::chrono::high_resolution_clock::now();
        g_timeRayHitsMesh += std::chrono::duration_cast<std::chrono::nanoseconds>(t_func_end - t_func_start).count();
#endif
        return false;
    }

    // Depth ordering: only commit if closer than the current global best (in world metric)
    if (worldDistance >= t_min) {
#if PROFILE_PERF
        auto t_back_end = std::chrono::high_resolution_clock::now();
        g_timeBackTransform += std::chrono::duration_cast<std::chrono::nanoseconds>(t_back_end - t_back_start).count();
        auto t_func_end = std::chrono::high_resolution_clock::now();
        g_timeRayHitsMesh += std::chrono::duration_cast<std::chrono::nanoseconds>(t_func_end - t_func_start).count();
#endif
        return false;
    }

    // Determine if we need to flip normals (for negative scale/reflection)
    bool shouldFlipNormals = transformMatrix ? hasNegativeScale(*transformMatrix) : mesh.hasNegativeScale;
    
    // Fast transform and normalize geometric normal
    VectorFloatTriplet worldGeomNormal = transformNormalFast(
        normPtr,
        localIntersection.geometricNormal.x,
        localIntersection.geometricNormal.y,
        localIntersection.geometricNormal.z);
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
        localIntersection.faceIndex >= 0 && localIntersection.faceIndex < (int)mesh.faces.size()) {
        
        const auto& face = mesh.faces[localIntersection.faceIndex];
        if (face.x < (int)scene->meshVertexNormals[meshIndex].size() &&
            face.y < (int)scene->meshVertexNormals[meshIndex].size() &&
            face.z < (int)scene->meshVertexNormals[meshIndex].size()) {
            
            double u = localIntersection.beta;
            double v = localIntersection.gamma;
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
            worldShadingNormal = transformNormalFast(
                normPtr,
                localIntersection.shadingNormal.x,
                localIntersection.shadingNormal.y,
                localIntersection.shadingNormal.z);
            double slen = sqrtf(worldShadingNormal.x*worldShadingNormal.x + 
                              worldShadingNormal.y*worldShadingNormal.y + 
                              worldShadingNormal.z*worldShadingNormal.z);
            double sinv = 1.0 / slen;
            worldShadingNormal.x *= sinv;
            worldShadingNormal.y *= sinv;
            worldShadingNormal.z *= sinv;
        }
    } else {
        worldShadingNormal = transformNormalFast(
            normPtr,
            localIntersection.shadingNormal.x,
            localIntersection.shadingNormal.y,
            localIntersection.shadingNormal.z);
        double slen = sqrtf(worldShadingNormal.x*worldShadingNormal.x + 
                          worldShadingNormal.y*worldShadingNormal.y + 
                          worldShadingNormal.z*worldShadingNormal.z);
        double sinv = 1.0 / slen;
        worldShadingNormal.x *= sinv;
        worldShadingNormal.y *= sinv;
        worldShadingNormal.z *= sinv;
    }
    
    // NOTE: Do NOT flip normals here to face the camera!
    // Dielectric materials need the original geometric normal direction 
    // to determine if the ray is entering or exiting the medium.
    // The shading code will handle normal orientation per material type.
    
    // Commit this mesh hit as the new global best
    // Commit this mesh hit as the new global best
    intersection.hit = true;  // CRITICAL: Mark as hit!
    intersection.point = worldPoint;
    intersection.geometricNormal = worldGeomNormal;
    intersection.shadingNormal = worldShadingNormal;
    intersection.distance = worldDistance;
    intersection.kind = Intersection::Kind::Mesh;
    intersection.containerIndex = -2;
    intersection.faceIndex = localIntersection.faceIndex;
    intersection.beta = localIntersection.beta;
    intersection.gamma = localIntersection.gamma;
    intersection.material = localIntersection.material;
    t_min = worldDistance;
    
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
    
    // NOTE: Near-plane clipping disabled. The nearDistance field defines the image plane
    // for ray generation, not a clipping plane. Objects closer than nearDistance should
    // still be visible (they're just between the camera and the image plane).
    // If you need near-plane clipping, enable this code and ensure scenes are designed accordingly.
    double minDistance = 0.0;
    bool isPrimaryRay = (ray.depth == 0 && !ray.shadowRay && !ray.reflectionRay && !ray.refractionRay);
    
#if PROFILE_PERF
    auto t_planes_start = std::chrono::high_resolution_clock::now();
#endif
    for(int i = 0; i < (int)scene.planes.size(); i++) {
        const Plane& plane = scene.planes[i];
        Ray testRay = plane.hasMotionBlur ? applyMotionBlurToRay(ray, plane.motionBlur) : ray;
        bool thisHit = rayHitsPlane(testRay, plane, scene.vertices, t_min, intersection, i, minDistance);
        if (thisHit && plane.hasMotionBlur) {
            finalizeMotionBlurHit(ray, plane.motionBlur, ray.time, intersection, t_min);
        }
        hit = thisHit || hit;
    }
#if PROFILE_PERF
    auto t_planes_end = std::chrono::high_resolution_clock::now();
    g_timeIntersectPlanes += std::chrono::duration_cast<std::chrono::nanoseconds>(t_planes_end - t_planes_start).count();
    
    auto t_tris_start = std::chrono::high_resolution_clock::now();
#endif
    for(int i = 0; i < (int)scene.triangles.size(); i++) {
        const Triangle& tri = scene.triangles[i];
        Ray baseRay = tri.hasMotionBlur ? applyMotionBlurToRay(ray, tri.motionBlur) : ray;
        
        if (tri.hasTransformation) {
            // Intersect in object space first, using a LOCAL t_min and intersection
            Ray objectRay;
            objectRay.origin = transformPoint(*tri.inverseTransformMatrix, baseRay.origin);
            objectRay.direction = normalize(transformDirection(*tri.inverseTransformMatrix, baseRay.direction));
            objectRay.depth = baseRay.depth;
            objectRay.shadowRay = baseRay.shadowRay;
            objectRay.reflectionRay = baseRay.reflectionRay;
            objectRay.refractionRay = baseRay.refractionRay;
            objectRay.time = baseRay.time;
            objectRay.random1 = baseRay.random1;
            objectRay.random2 = baseRay.random2;

            double local_t_min = std::numeric_limits<double>::max();
            Intersection localIntersection;

            if (rayHitsTriangle(objectRay,
                                tri.indices,
                                scene.vertices,
                                local_t_min,
                                localIntersection,
                                scene.intersectionTestEpsilon,
                                0.0,
                                tri.material,
                                scene.enableBackFaceCulling,
                                -1,
                                i,
                                0.0)) {

                // Back-transform hit to world space
                VectorFloatTriplet worldPoint = transformPoint(*tri.transformMatrix, localIntersection.point);
                VectorFloatTriplet worldGeomNormal = normalize(transformNormal(*tri.normalMatrix, localIntersection.geometricNormal));
                VectorFloatTriplet worldShadingNormal = normalize(transformNormal(*tri.normalMatrix, localIntersection.shadingNormal));

                double worldDistance = recomputeDistanceFromOrigin(worldPoint, baseRay.origin);

                // Near-plane clipping in world space for primary rays
                if (isPrimaryRay && worldDistance < minDistance) {
                    continue;
                }

                // Depth ordering in world space (global metric)
                if (worldDistance < t_min) {
                    intersection = localIntersection;
                    intersection.point = worldPoint;
                    intersection.geometricNormal = worldGeomNormal;
                    intersection.shadingNormal = worldShadingNormal;
                    intersection.distance = worldDistance;
                    intersection.kind = Intersection::Kind::Triangle;
                    t_min = worldDistance;
                    hit = true;

                    if (tri.hasMotionBlur) {
                        finalizeMotionBlurHit(ray, tri.motionBlur, ray.time, intersection, t_min);
                    }
                }
            }
        } else {
            bool thisHit = rayHitsTriangle(baseRay, tri.indices, scene.vertices, t_min, intersection, scene.intersectionTestEpsilon, scene.cameraTriangleDeterminant[scene.currentCameraIndex][i], tri.material, scene.enableBackFaceCulling, -1, i, minDistance);
            if (thisHit) {
                hit = true;
                intersection.kind = Intersection::Kind::Triangle;
                if (tri.hasMotionBlur) {
                    finalizeMotionBlurHit(ray, tri.motionBlur, ray.time, intersection, t_min);
                }
            }
        }
    }
#if PROFILE_PERF
    auto t_tris_end = std::chrono::high_resolution_clock::now();
    g_timeIntersectTriangles += std::chrono::duration_cast<std::chrono::nanoseconds>(t_tris_end - t_tris_start).count();
    
    auto t_meshes_start = std::chrono::high_resolution_clock::now();
#endif
    for(int i = 0; i < (int)scene.meshes.size(); i++) {
        const Mesh& mesh = scene.meshes[i];
        Ray testRay = mesh.hasMotionBlur ? applyMotionBlurToRay(ray, mesh.motionBlur) : ray;
        MeshBVH* bvh = (i < scene.meshBVHs.size()) ? scene.meshBVHs[i] : nullptr;
        bool thisHit = rayHitsMesh(testRay, mesh, scene.vertices, scene.cameraMeshDeterminant[scene.currentCameraIndex][i], t_min, intersection, scene.intersectionTestEpsilon, bvh, scene.enableBackFaceCulling, i, nullptr, nullptr, nullptr, &scene, nullptr, nullptr, minDistance);
        if (thisHit && mesh.hasMotionBlur) {
            finalizeMotionBlurHit(ray, mesh.motionBlur, ray.time, intersection, t_min);
        }
        hit = thisHit || hit;
    }
#if PROFILE_PERF
    auto t_meshes_end = std::chrono::high_resolution_clock::now();
    g_timeIntersectMeshes += std::chrono::duration_cast<std::chrono::nanoseconds>(t_meshes_end - t_meshes_start).count();
    
    auto t_spheres_start = std::chrono::high_resolution_clock::now();
#endif
    for(int i = 0; i < (int)scene.spheres.size(); i++) {
        const Sphere& sphere = scene.spheres[i];
        Ray testRay = sphere.hasMotionBlur ? applyMotionBlurToRay(ray, sphere.motionBlur) : ray;
        bool thisHit = rayHitsSphere(testRay, sphere, scene.vertices, t_min, intersection, i, minDistance);
        if (thisHit && sphere.hasMotionBlur) {
            finalizeMotionBlurHit(ray, sphere.motionBlur, ray.time, intersection, t_min);
        }
        hit = thisHit || hit;
    }
#if PROFILE_PERF
    auto t_spheres_end = std::chrono::high_resolution_clock::now();
    g_timeIntersectSpheres += std::chrono::duration_cast<std::chrono::nanoseconds>(t_spheres_end - t_spheres_start).count();
    
    auto t_instances_start = std::chrono::high_resolution_clock::now();
#endif
    
    for(int i = 0; i < (int)scene.meshInstances.size(); i++) {
        const MeshInstance& instance = scene.meshInstances[i];
        if (instance.baseMesh && instance.baseMeshIndex >= 0) {
            Ray testRay = instance.hasMotionBlur ? applyMotionBlurToRay(ray, instance.motionBlur) : ray;
            
            // OPTIMIZATION: Early rejection using world-space bounding box
            // Note: for motion blur, AABB check uses offset ray but this is conservative
            if (instance.worldSpaceBounds != nullptr) {
                if (!instance.worldSpaceBounds->intersect(testRay, 0.0, t_min)) {
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
            
            bool instanceHit = rayHitsMesh(testRay, *instance.baseMesh, scene.vertices, 
                                         scene.cameraMeshDeterminant[scene.currentCameraIndex][baseMeshIndex], 
                                         t_min, intersection, scene.intersectionTestEpsilon, 
                                         bvh, scene.enableBackFaceCulling, baseMeshIndex,
                                         instance.transformMatrix, instance.inverseTransformMatrix, instance.normalMatrix, &scene,
                                         instance.material,
                                         instance.worldSpaceBounds,
                                         minDistance);
            
            if (instanceHit) {
                hit = true;
                if (instance.hasMotionBlur) {
                    finalizeMotionBlurHit(ray, instance.motionBlur, ray.time, intersection, t_min);
                }
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
                
                // NOTE: Do NOT flip normals here to face the camera!
                // Dielectric materials need the original geometric normal direction 
                // to determine if the ray is entering or exiting the medium.
                // Just ensure shading normal is aligned with geometric normal.
                if (dotProduct(intersection.shadingNormal, intersection.geometricNormal) < 0.0) {
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
    
    // Test LightSpheres (emissive spheres)
    for(int i = 0; i < (int)scene.lightSpheres.size(); i++) {
        const LightSphere& lightSphere = scene.lightSpheres[i];
        // Create a temporary Sphere for intersection testing
        Sphere tempSphere;
        tempSphere.center = lightSphere.center;
        tempSphere.radius = lightSphere.radius;
        tempSphere.hasMotionBlur = lightSphere.hasMotionBlur;
        tempSphere.motionBlur = lightSphere.motionBlur;
        tempSphere.hasTransformation = lightSphere.hasTransformation;
        tempSphere.transformMatrix = lightSphere.transformMatrix;
        tempSphere.inverseTransformMatrix = lightSphere.inverseTransformMatrix;
        tempSphere.normalMatrix = lightSphere.normalMatrix;
        tempSphere.material = lightSphere.material;
        
        Ray testRay = lightSphere.hasMotionBlur ? applyMotionBlurToRay(ray, lightSphere.motionBlur) : ray;
        Intersection tempIntersection;
        double temp_t_min = t_min;
        bool thisHit = rayHitsSphere(testRay, tempSphere, scene.vertices, temp_t_min, tempIntersection, -10000 - i, minDistance);  // Use negative index to identify as light
        if (thisHit && temp_t_min < t_min) {
            intersection = tempIntersection;
            intersection.kind = Intersection::Kind::LightSphere;
            intersection.containerIndex = i;  // Store lightSphere index
            t_min = temp_t_min;
            hit = true;
        }
        if (thisHit && lightSphere.hasMotionBlur) {
            finalizeMotionBlurHit(ray, lightSphere.motionBlur, ray.time, intersection, t_min);
        }
    }
    
    // Test LightMeshes (emissive meshes)
    for(int i = 0; i < (int)scene.lightMeshes.size(); i++) {
        const LightMesh& lightMesh = scene.lightMeshes[i];
        // Create a temporary Mesh for intersection testing
        Mesh tempMesh;
        tempMesh.faces = lightMesh.faces;
        tempMesh.shadingMode = lightMesh.shadingMode;
        tempMesh.hasMotionBlur = lightMesh.hasMotionBlur;
        tempMesh.motionBlur = lightMesh.motionBlur;
        tempMesh.hasTransformation = lightMesh.hasTransformation;
        tempMesh.transformMatrix = lightMesh.transformMatrix;
        tempMesh.inverseTransformMatrix = lightMesh.inverseTransformMatrix;
        tempMesh.normalMatrix = lightMesh.normalMatrix;
        tempMesh.material = lightMesh.material;
        
        Ray testRay = lightMesh.hasMotionBlur ? applyMotionBlurToRay(ray, lightMesh.motionBlur) : ray;
        // Use a simple BVH or direct testing - for now, test without BVH
        MeshBVH* bvh = nullptr;  // LightMeshes don't have BVH yet
        Intersection tempIntersection;
        double temp_t_min = t_min;
        // Get determinant for light mesh (use 0.0 as fallback)
        // cameraMeshDeterminant is [camera][mesh][determinant], but for light meshes we don't have precomputed determinants
        // So we'll use 0.0 and let rayHitsMesh compute it
        std::vector<double> emptyDeterminants;  // Empty vector - rayHitsMesh will compute
        
        bool thisHit = rayHitsMesh(testRay, tempMesh, scene.vertices, 
                                    emptyDeterminants,
                                    temp_t_min, tempIntersection, scene.intersectionTestEpsilon, bvh, 
                                    false, -20000 - i,  // Disable back-face culling for emissive objects
                                    nullptr, nullptr, nullptr, &scene, nullptr, nullptr, minDistance);
        if (thisHit && temp_t_min < t_min) {
            intersection = tempIntersection;
            intersection.kind = Intersection::Kind::LightMesh;
            intersection.containerIndex = i;  // Store lightMesh index
            t_min = temp_t_min;
            hit = true;
        }
        if (thisHit && lightMesh.hasMotionBlur) {
            finalizeMotionBlurHit(ray, lightMesh.motionBlur, ray.time, intersection, t_min);
        }
    }
    
#if PROFILE_PERF
    auto t_end = std::chrono::high_resolution_clock::now();
    g_timeIntersect += std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count();
#endif
    
    return intersection;
}

void orthonormalBasis(const VectorFloatTriplet& n, VectorFloatTriplet& u, VectorFloatTriplet& v) {
    if (std::abs(n.x) > std::abs(n.y)) {
        double invLen = 1.0 / std::sqrt(n.x * n.x + n.z * n.z);
        u = VectorFloatTriplet{-n.z * invLen, 0.0, n.x * invLen};
    } else {
        double invLen = 1.0 / std::sqrt(n.y * n.y + n.z * n.z);
        u = VectorFloatTriplet{0.0, -n.z * invLen, n.y * invLen};
    }
    v = crossProduct(n, u);
    u = normalize(u);
    v = normalize(v);
}

VectorFloatTriplet sampleHemisphereUniform(const VectorFloatTriplet& N, double xi1, double xi2) {
    // Uniform hemisphere sampling
    // phi = 2*pi*xi1, theta = arccos(xi2)
    // PDF = 1/(2*pi)
    double phi = 2.0 * M_PI * xi1;
    double cosTheta = xi2;  // Uniform in [0,1]
    double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    
    // Sample in local coordinate system (upper hemisphere)
    double x = sinTheta * std::cos(phi);
    double y = cosTheta;  // y is up
    double z = sinTheta * std::sin(phi);
    
    // Transform to world space using surface normal
    VectorFloatTriplet u, v;
    orthonormalBasis(N, u, v);
    VectorFloatTriplet direction = u * x + N * y + v * z;
    return normalize(direction);
}

VectorFloatTriplet sampleHemisphereCosine(const VectorFloatTriplet& N, double xi1, double xi2) {
    // Cosine-weighted hemisphere sampling
    // phi = 2*pi*xi1, theta = arcsin(sqrt(xi2))
    // PDF = cos(theta)/pi
    double phi = 2.0 * M_PI * xi1;
    double cosTheta = std::sqrt(xi2);  // Cosine-weighted
    double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    
    // Sample in local coordinate system (upper hemisphere)
    double x = sinTheta * std::cos(phi);
    double y = cosTheta;  // y is up (cosine-weighted)
    double z = sinTheta * std::sin(phi);
    
    // Transform to world space using surface normal
    VectorFloatTriplet u, v;
    orthonormalBasis(N, u, v);
    VectorFloatTriplet direction = u * x + N * y + v * z;
    return normalize(direction);
}

void precomputeLightMeshSampling(LightMesh& mesh, const vector<VectorFloatTriplet>& vertices) {
    // Compute triangle areas and build CDF for importance sampling
    mesh.totalArea = 0.0;
    mesh.cdfTriangleAreas.clear();
    
    for (const auto& face : mesh.faces) {
        const VectorFloatTriplet& v0 = vertices[face.x];
        const VectorFloatTriplet& v1 = vertices[face.y];
        const VectorFloatTriplet& v2 = vertices[face.z];
        
        VectorFloatTriplet e1 = v1 - v0;
        VectorFloatTriplet e2 = v2 - v0;
        VectorFloatTriplet cross = crossProduct(e1, e2);
        double area = 0.5 * std::sqrt(dotProduct(cross, cross));
        
        mesh.totalArea += area;
        mesh.cdfTriangleAreas.push_back(mesh.totalArea);
    }
    
    // Normalize CDF to [0,1]
    if (mesh.totalArea > 1e-10) {
        for (auto& cdf : mesh.cdfTriangleAreas) {
            cdf /= mesh.totalArea;
        }
    }
}

VectorFloatTriplet sampleLightMesh(const LightMesh& mesh,
                                    const vector<VectorFloatTriplet>& vertices,
                                    double xi1, double xi2, double xi3,
                                    VectorFloatTriplet& lightNormal,
                                    double& pdf) {
    // Step 1: Select triangle based on area (using CDF)
    int triangleIndex = 0;
    if (!mesh.cdfTriangleAreas.empty()) {
        for (size_t i = 0; i < mesh.cdfTriangleAreas.size(); i++) {
            if (xi1 <= mesh.cdfTriangleAreas[i]) {
                triangleIndex = i;
                break;
            }
        }
    }
    
    if (triangleIndex >= (int)mesh.faces.size()) {
        triangleIndex = mesh.faces.size() - 1;
    }
    
    const VectorIntTriplet& face = mesh.faces[triangleIndex];
    const VectorFloatTriplet& v0 = vertices[face.x];
    const VectorFloatTriplet& v1 = vertices[face.y];
    const VectorFloatTriplet& v2 = vertices[face.z];
    
    // Step 2: Sample point uniformly on triangle using barycentric coordinates
    double sqrt_xi2 = std::sqrt(xi2);
    double u = 1.0 - sqrt_xi2;
    double v = xi3 * sqrt_xi2;
    double w = 1.0 - u - v;
    
    VectorFloatTriplet samplePoint = v0 * u + v1 * v + v2 * w;
    
    // Compute triangle normal
    VectorFloatTriplet e1 = v1 - v0;
    VectorFloatTriplet e2 = v2 - v0;
    VectorFloatTriplet normal = crossProduct(e1, e2);
    double normalLen = std::sqrt(dotProduct(normal, normal));
    if (normalLen > 1e-10) {
        lightNormal = normal * (1.0 / normalLen);
    } else {
        lightNormal = VectorFloatTriplet{0, 1, 0};  // Fallback
    }
    
    // Compute area PDF: probability of selecting this triangle * probability of point on triangle
    double triangleArea = 0.5 * normalLen;
    if (mesh.totalArea > 1e-10 && triangleArea > 1e-10) {
        double triangleProb = triangleArea / mesh.totalArea;
        double pointProb = 1.0 / triangleArea;  // Uniform on triangle
        pdf = triangleProb * pointProb;  // = 1 / totalArea
    } else {
        pdf = 0.0;
    }
    
    return samplePoint;
}

VectorFloatTriplet sampleLightSphere(const LightSphere& sphere,
                                      const vector<VectorFloatTriplet>& vertices,
                                      const VectorFloatTriplet& shadingPoint,
                                      double xi1, double xi2,
                                      VectorFloatTriplet& lightNormal,
                                      double& pdf) {
    // Get sphere center and radius
    VectorFloatTriplet center = vertices[sphere.center];
    double radius = sphere.radius;
    
    // Sample uniformly on sphere surface using spherical coordinates
    // Uniform sampling: phi = 2*pi*xi1, theta = arccos(1 - 2*xi2)
    double phi = 2.0 * M_PI * xi1;
    double cosTheta = 1.0 - 2.0 * xi2;  // Uniform in [-1, 1]
    double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    
    // Point on unit sphere in local coordinates
    VectorFloatTriplet localPoint;
    localPoint.x = sinTheta * std::cos(phi);
    localPoint.y = cosTheta;
    localPoint.z = sinTheta * std::sin(phi);
    
    // Transform to world space: scale by radius and translate by center
    VectorFloatTriplet samplePoint = center + localPoint * radius;
    
    // Normal points from center to surface point
    lightNormal = normalize(localPoint);
    
    // Area PDF: 1 / (4 * pi * r^2) for uniform sphere sampling
    double sphereArea = 4.0 * M_PI * radius * radius;
    if (sphereArea > 1e-10) {
        pdf = 1.0 / sphereArea;
    } else {
        pdf = 0.0;
    }
    
    return samplePoint;
}

double areaPDFToSolidAnglePDF(double areaPDF, double distance, double cosAtLight) {
    // Convert area PDF to solid angle PDF
    // p(w) = p(x) * r^2 / |cos(theta_light)|
    if (cosAtLight < 1e-10 || distance < 1e-10) {
        return 0.0;
    }
    double distanceSq = distance * distance;
    return areaPDF * distanceSq / std::abs(cosAtLight);
}

double misWeight(double pdf1, double pdf2, const std::string& heuristic) {
    if (pdf1 <= 0.0 && pdf2 <= 0.0) {
        return 0.0;
    }
    if (pdf1 <= 0.0) {
        return 0.0;
    }
    if (pdf2 <= 0.0) {
        return 1.0;
    }
    
    if (heuristic == "balance") {
        return pdf1 / (pdf1 + pdf2);
    } else if (heuristic == "power") {
        double pdf1Sq = pdf1 * pdf1;
        double pdf2Sq = pdf2 * pdf2;
        return pdf1Sq / (pdf1Sq + pdf2Sq);
    } else if (heuristic == "01") {
        return (pdf1 > pdf2) ? 1.0 : 0.0;
    }
    // Default: balance heuristic
    return pdf1 / (pdf1 + pdf2);
}

VectorFloatTriplet sampleDirectLight(const Scene& scene,
                                      const VectorFloatTriplet& shadingPoint,
                                      const VectorFloatTriplet& shadingNormal,
                                      double xi1, double xi2, double xi3, double xi4,
                                      VectorFloatTriplet& lightDir,
                                      double& pdfLight) {
    // Count total number of light sources (point lights, area lights, object lights)
    int numPointLights = scene.pointLights.size();
    int numAreaLights = scene.areaLights.size();
    int numLightSpheres = scene.lightSpheres.size();
    int numLightMeshes = scene.lightMeshes.size();
    int totalLights = numPointLights + numAreaLights + numLightSpheres + numLightMeshes;
    
    if (totalLights == 0) {
        pdfLight = 0.0;
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }
    
    // Select a light uniformly
    double lightSelect = xi1 * totalLights;
    int lightIndex = (int)std::floor(lightSelect);
    if (lightIndex >= totalLights) {
        lightIndex = totalLights - 1;
    }
    
    VectorFloatTriplet lightPoint;
    VectorFloatTriplet lightNormal;
    VectorFloatTriplet radiance;
    double areaPDF = 0.0;
    
    // Sample from selected light
    if (lightIndex < numPointLights) {
        // Point light - sample the point light position
        const PointLight& light = scene.pointLights[lightIndex];
        lightPoint = light.position;
        lightNormal = VectorFloatTriplet{0, 0, 0};  // Point lights don't have normals
        radiance = light.intensity;  // Point lights use intensity, not radiance
        areaPDF = 1.0;  // Dirac delta - will be handled specially
        // For point lights, PDF conversion is different
    } else if (lightIndex < numPointLights + numAreaLights) {
        // Area light
        const AreaLight& light = scene.areaLights[lightIndex - numPointLights];
        VectorFloatTriplet u, v;
        orthonormalBasis(normalize(light.normal), u, v);
        double halfSize = light.size / 2.0;
        double offsetU = (xi2 * 2.0 - 1.0) * halfSize;
        double offsetV = (xi3 * 2.0 - 1.0) * halfSize;
        lightPoint = light.position + u * offsetU + v * offsetV;
        lightNormal = normalize(light.normal);
        radiance = light.radiance;
        double area = light.size * light.size;
        areaPDF = 1.0 / area;
    } else if (lightIndex < numPointLights + numAreaLights + numLightSpheres) {
        // LightSphere
        const LightSphere& light = scene.lightSpheres[lightIndex - numPointLights - numAreaLights];
        lightPoint = sampleLightSphere(light, scene.vertices, shadingPoint, xi2, xi3, lightNormal, areaPDF);
        radiance = light.radiance;
    } else {
        // LightMesh
        const LightMesh& light = scene.lightMeshes[lightIndex - numPointLights - numAreaLights - numLightSpheres];
        lightPoint = sampleLightMesh(light, scene.vertices, xi2, xi3, xi4, lightNormal, areaPDF);
        radiance = light.radiance;
    }
    
    // Compute direction from shading point to light point
    VectorFloatTriplet toLight = lightPoint - shadingPoint;
    double distanceSq = dotProduct(toLight, toLight);
    double distance = std::sqrt(distanceSq);
    
    if (distance < 1e-10) {
        pdfLight = 0.0;
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }
    
    lightDir = toLight * (1.0 / distance);
    
    // Check visibility
    VectorFloatTriplet offsetNormal = shadingNormal;
    if (dotProduct(offsetNormal, lightDir) < 0.0) {
        offsetNormal = -offsetNormal;
    }
    Ray shadowRay(shadingPoint + scene.shadowRayEpsilon * offsetNormal, lightDir, 0, true, false, false, 0.0);
    Intersection shadowHit = intersect(scene, shadowRay);
    if (shadowHit.hit && shadowHit.distance < distance - scene.shadowRayEpsilon) {
        // Occluded
        pdfLight = 0.0;
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }
    
    // Compute cosine term at shading point
    double cosTheta = std::max(0.0, dotProduct(shadingNormal, lightDir));
    if (cosTheta <= 0.0) {
        pdfLight = 0.0;
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }
    
    if (numPointLights > 0 && lightIndex < numPointLights) {
        // Point light: treat as delta distribution with uniform light selection
        pdfLight = 1.0 / totalLights;
        // Apply distance falloff here; cosine handled by caller
        return radiance * (1.0 / distanceSq);
    }
    
    // Convert area PDF to solid angle PDF for area/object lights
    double cosAtLight = std::abs(dotProduct(lightNormal, -lightDir));
    if (areaPDF > 0.0 && cosAtLight > 1e-10) {
        pdfLight = areaPDFToSolidAnglePDF(areaPDF, distance, cosAtLight) / totalLights;
    } else {
        pdfLight = 0.0;
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }
    
    // Return radiance (cosine handled by caller)
    return radiance;
}

VectorFloatTriplet perturbDirection(const VectorFloatTriplet& idealDir,
                                    double roughness,
                                    double random1,
                                    double random2) {
    if (roughness <= 0.0) return idealDir;
    
    VectorFloatTriplet u, v;
    orthonormalBasis(idealDir, u, v);
    
    double phi = 2.0 * M_PI * random1;
    double cosTheta = 1.0 - random2 * (1.0 - std::cos(roughness));
    double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    
    VectorFloatTriplet perturbed = u * (sinTheta * std::cos(phi)) 
                                 + v * (sinTheta * std::sin(phi)) 
                                 + idealDir * cosTheta;
    return normalize(perturbed);
}

static std::vector<unsigned int> getTextureIdsFromIntersection(const Scene& scene, const Intersection& intersection) {
    std::vector<unsigned int> textureIds;
    
    switch (intersection.kind) {
        case Intersection::Kind::Triangle:
            if (intersection.containerIndex >= 0 && intersection.containerIndex < (int)scene.triangles.size()) {
                textureIds = scene.triangles[intersection.containerIndex].textureIds;
            }
            break;
        case Intersection::Kind::Mesh: {
            if (intersection.containerIndex >= 0 && intersection.containerIndex < (int)scene.meshes.size()) {
                textureIds = scene.meshes[intersection.containerIndex].textureIds;
            } else if (intersection.containerIndex == -2) {
                for (size_t i = 0; i < scene.meshes.size(); i++) {
                    if (intersection.faceIndex >= 0 && 
                        intersection.faceIndex < (int)scene.meshes[i].faces.size() &&
                        scene.meshes[i].material == intersection.material) {
                        textureIds = scene.meshes[i].textureIds;
                        break;
                    }
                }
                if (textureIds.empty()) {
                    for (size_t i = 0; i < scene.meshes.size(); i++) {
                        if (intersection.faceIndex >= 0 && 
                            intersection.faceIndex < (int)scene.meshes[i].faces.size()) {
                            textureIds = scene.meshes[i].textureIds;
                            break;
                        }
                    }
                }
            }
            break;
        }
        case Intersection::Kind::Sphere:
            if (intersection.containerIndex >= 0 && intersection.containerIndex < (int)scene.spheres.size()) {
                textureIds = scene.spheres[intersection.containerIndex].textureIds;
            }
            break;
        case Intersection::Kind::Plane:
            if (intersection.containerIndex >= 0 && intersection.containerIndex < (int)scene.planes.size()) {
                textureIds = scene.planes[intersection.containerIndex].textureIds;
            }
            break;
        default:
            break;
    }
    
    return textureIds;
}

VectorFloatTriplet computeShading(const Scene& scene, Ray& ray, const Intersection& intersection) {
    Material* material = intersection.material;
    VectorFloatTriplet color{0.0, 0.0, 0.0};

    if (!material) {
        return color;
    }

    std::vector<unsigned int> textureIds = getTextureIdsFromIntersection(scene, intersection);
    
    VectorFloatTriplet objectSpacePoint = intersection.point;
    VectorFloatTriplet objectSpaceNormal = intersection.geometricNormal;
    
    // For transformed objects, compute object-space coordinates for proper texture mapping
    if (intersection.kind == Intersection::Kind::Sphere && 
        intersection.containerIndex >= 0 && 
        intersection.containerIndex < (int)scene.spheres.size()) {
        const Sphere& sphere = scene.spheres[intersection.containerIndex];
        if (sphere.hasTransformation && sphere.inverseTransformMatrix) {
            // Transform world point to object space
            objectSpacePoint = transformPoint(*sphere.inverseTransformMatrix, intersection.point);
            // Compute object-space normal from object-space point and sphere center
            VectorFloatTriplet sphereCenter = scene.vertices[sphere.center];
            objectSpaceNormal = normalize(objectSpacePoint - sphereCenter);
        }
    } else if (intersection.kind == Intersection::Kind::Mesh) {
        int meshIdx = intersection.containerIndex;
        if (meshIdx == -2) {
            // Find mesh by material
            for (size_t i = 0; i < scene.meshes.size(); i++) {
                if (intersection.faceIndex >= 0 && 
                    intersection.faceIndex < (int)scene.meshes[i].faces.size() &&
                    scene.meshes[i].material == intersection.material) {
                    meshIdx = i;
                    break;
                }
            }
        }
        if (meshIdx >= 0 && meshIdx < (int)scene.meshes.size()) {
            const Mesh& mesh = scene.meshes[meshIdx];
            if (mesh.hasTransformation && mesh.inverseTransformMatrix) {
                objectSpacePoint = transformPoint(*mesh.inverseTransformMatrix, intersection.point);
            }
        }
    } else if (intersection.kind == Intersection::Kind::Plane &&
               intersection.containerIndex >= 0 &&
               intersection.containerIndex < (int)scene.planes.size()) {
        const Plane& plane = scene.planes[intersection.containerIndex];
        if (plane.hasTransformation && plane.inverseTransformMatrix) {
            objectSpacePoint = transformPoint(*plane.inverseTransformMatrix, intersection.point);
        }
    } else if (intersection.kind == Intersection::Kind::Triangle &&
               intersection.containerIndex >= 0 &&
               intersection.containerIndex < (int)scene.triangles.size()) {
        const Triangle& tri = scene.triangles[intersection.containerIndex];
        if (tri.hasTransformation && tri.inverseTransformMatrix) {
            objectSpacePoint = transformPoint(*tri.inverseTransformMatrix, intersection.point);
        }
    }
    
    // Compute UV coordinates using object-space normal for spheres
    VectorFloatPair uv = computeUVCoordinates(intersection, scene);
    
    // Override UV for spheres with object-space normal
    if (intersection.kind == Intersection::Kind::Sphere) {
        // Spherical UV mapping: 
        // u wraps around Y axis, v goes from north to south pole
        // Negate atan2 to match common texture orientation (Americas on left side)
        double u = -atan2(objectSpaceNormal.z, objectSpaceNormal.x) / (2.0 * M_PI) + 0.5;
        double v = acos(std::max(-1.0, std::min(1.0, objectSpaceNormal.y))) / M_PI;
        // Ensure UV is in [0,1] range (handle wrapping)
        u = u - floor(u);
        uv.x = u;
        uv.y = v;
    }
    
    // Create a mutable copy of material properties for texture application
    VectorFloatTriplet ambientReflectance = material->ambientReflectance;
    VectorFloatTriplet diffuseReflectance = material->diffuseReflectance;
    VectorFloatTriplet specularReflectance = material->specularReflectance;
    VectorFloatTriplet normal = normalize(intersection.shadingNormal);
    
    for (unsigned int textureId : textureIds) {
        const TextureMap* textureMap = scene.getTextureMapById(textureId);
        if (!textureMap) continue;
        
        VectorFloatTriplet texturePosition = (textureMap->type == "perlin" || textureMap->type == "checkerboard") 
                                              ? objectSpacePoint : intersection.point;
        VectorFloatTriplet textureValue = sampleTexture(textureMap, uv, texturePosition, &scene, false);
        
        if (textureMap->decalMode == DecalMode::ReplaceNormal) {
            VectorFloatTriplet tangent, bitangent;
            VectorFloatTriplet up{0.0, 1.0, 0.0};
            VectorFloatTriplet right{1.0, 0.0, 0.0};
            if (std::abs(dotProduct(intersection.geometricNormal, up)) < 0.999) {
                tangent = normalize(crossProduct(intersection.geometricNormal, up));
            } else {
                tangent = normalize(crossProduct(intersection.geometricNormal, right));
            }
            bitangent = normalize(crossProduct(intersection.geometricNormal, tangent));
            normal = transformNormalFromTangentSpace(textureValue, tangent, bitangent, intersection.geometricNormal);
        } else if (textureMap->decalMode == DecalMode::BumpNormal) {
            VectorFloatTriplet bumpPos = (textureMap->type == "perlin" || textureMap->type == "checkerboard")
                ? objectSpacePoint
                : intersection.point;
            
            VectorFloatTriplet bumpN = intersection.geometricNormal;
            if (intersection.kind == Intersection::Kind::Sphere &&
                (textureMap->type == "perlin" || textureMap->type == "checkerboard")) {
                bumpN = objectSpaceNormal;
            }
            
            normal = applyBumpMapping(textureMap, uv, bumpPos, bumpN, &scene, textureMap->bumpFactor);
        }
    }
    
    if (material->type != "dielectric") {
        if (dotProduct(normal, ray.direction) > 0.0) {
            normal = -normal;
        }
    }
    
    bool replaceAllMode = false;
    VectorFloatTriplet replaceAllColor{0.0, 0.0, 0.0};
    
    for (unsigned int textureId : textureIds) {
        const TextureMap* textureMap = scene.getTextureMapById(textureId);
        if (!textureMap) continue;
        
        VectorFloatTriplet texturePosition = (textureMap->type == "perlin" || textureMap->type == "checkerboard") 
                                              ? objectSpacePoint : intersection.point;
        VectorFloatTriplet textureValue = sampleTexture(textureMap, uv, texturePosition, &scene, true);
        
        switch (textureMap->decalMode) {
            case DecalMode::ReplaceKd:
                diffuseReflectance = textureValue;
                break;
            case DecalMode::BlendKd:
                diffuseReflectance = (diffuseReflectance + textureValue) * 0.5;
                break;
            case DecalMode::ReplaceKs:
                specularReflectance = textureValue;
                break;
            case DecalMode::ReplaceAll:
                replaceAllMode = true;
                replaceAllColor = textureValue;
                break;
            default:
                break;
        }
    }
    
    if (replaceAllMode) {
        return replaceAllColor;
    }
    
    color += ambientReflectance * scene.ambientLight.intensity;
    if (material->type == "conductor") {
        VectorFloatTriplet viewDir = normalize(-ray.direction);
        double cosTheta = dotProduct(normal, viewDir);
        double F = fresnelConductor(
            cosTheta,
            material->refractionIndex,
            material->absorptionIndex
        );
        for (const PointLight& light : scene.pointLights) {
            if (isInShadow(scene, ray, light, intersection)) {
                continue;
            }

            VectorFloatTriplet lightVec = light.position - intersection.point;
            double distanceSq = dotProduct(lightVec, lightVec);
            double distance   = std::sqrt(distanceSq);
            VectorFloatTriplet lightDir = lightVec * (1.0 / distance);

            double attenuation = 1.0 / distanceSq;

            double NdotL = std::max(0.0, dotProduct(normal, lightDir));
            if (NdotL <= 0.0) continue;

            VectorFloatTriplet viewDir   = normalize(-ray.direction);
            VectorFloatTriplet halfVector = normalize(lightDir + viewDir);
            double NdotH = std::max(0.0, dotProduct(normal, halfVector));
            double specularFactor = std::pow(NdotH, material->phongExponent);

            VectorFloatTriplet specular =
                specularReflectance * light.intensity *
                (specularFactor * NdotL * attenuation);

            color += specular;
        }

        for (const AreaLight& light : scene.areaLights) {
            VectorFloatTriplet lightNormal = normalize(light.normal);
            VectorFloatTriplet u, v;
            orthonormalBasis(lightNormal, u, v);

            double halfSize = light.size * 0.5;

            double ksi_1 = (ray.random1 * 2.0 - 1.0) * halfSize;
            double ksi_2 = (ray.random2 * 2.0 - 1.0) * halfSize;
            VectorFloatTriplet samplePoint = light.position + u * ksi_1 + v * ksi_2;

            VectorFloatTriplet toLight = samplePoint - intersection.point;
            double distanceSq = dotProduct(toLight, toLight);
            double distance   = std::sqrt(distanceSq);
            VectorFloatTriplet lightDir = toLight * (1.0 / distance);

            double cosThetaSurf  = dotProduct(normal, lightDir);
            double cosThetaLight = dotProduct(lightNormal, -lightDir);
            cosThetaLight = std::abs(cosThetaLight);
            if (cosThetaSurf <= 0.0) continue;

            VectorFloatTriplet offsetNormal = intersection.geometricNormal;
            if (dotProduct(offsetNormal, lightDir) < 0.0) {
                offsetNormal = -offsetNormal;
            }
            Ray shadowRay{
                intersection.point + scene.shadowRayEpsilon * offsetNormal,
                lightDir, ray.depth + 1, true, false, false, ray.time
            };
            Intersection shadowHit = intersect(scene, shadowRay);
            if (shadowHit.hit && shadowHit.distance < distance - scene.shadowRayEpsilon) continue;

            double area = light.size * light.size;
            VectorFloatTriplet irradiance =
                light.radiance * (cosThetaSurf * cosThetaLight * area / distanceSq);

            VectorFloatTriplet halfVector = normalize(lightDir + viewDir);
            double NdotH = std::max(0.0, dotProduct(normal, halfVector));
            double specFactor = std::pow(NdotH, material->phongExponent);

            color += material->mirrorReflectance * irradiance * specFactor;
        }

        Ray reflectedRay = reflect(ray, normal,
                                    intersection.point,
                                    scene.shadowRayEpsilon);

        reflectedRay.direction = normalize(reflectedRay.direction);
        reflectedRay.shadowRay  = false;
        reflectedRay.reflectionRay = true;
        reflectedRay.refractionRay = false;

        if (material->roughness > 0.0) {
            reflectedRay.direction =
                perturbDirection(reflectedRay.direction,
                                 material->roughness,
                                 ray.random1,
                                 ray.random2);
        }

        Intersection reflectionIntersection = intersect(scene, reflectedRay);
        VectorFloatTriplet reflectedColor =
            computePixelColor(scene, reflectedRay, reflectionIntersection);

        // F = std::max(0.0, std::min(1.0, F));

        color += F * material->mirrorReflectance * reflectedColor;

        return color;
    }

    else if (material->type == "dielectric") {
        VectorFloatTriplet N = normalize(normal);
    
        double n1 = 1.0;                          // index of refraction of current medium (air)
        double n2 = material->refractionIndex;    // index of refraction of the material
    
        bool entering = dotProduct(ray.direction, N) < 0.0;
    
        if (!entering) {
            N = -N;
            std::swap(n1, n2);
        }
    
        double cosThetaI = dotProduct(-ray.direction, N);
        cosThetaI = std::max(-1.0, std::min(1.0, cosThetaI));
    
        // Fresnel reflectance for dielectrics
        double F = fresnelDielectric(cosThetaI, n1, n2);
    
        Ray reflectedRay = reflect(ray, N, intersection.point, scene.shadowRayEpsilon);
    
        if (material->roughness > 0.0) {
            reflectedRay.direction = perturbDirection(reflectedRay.direction,
                                                      material->roughness,
                                                      ray.random1,
                                                      ray.random2);
        }
    
        Intersection reflectionIntersection = intersect(scene, reflectedRay);
        VectorFloatTriplet reflectedColor =
            computePixelColor(scene, reflectedRay, reflectionIntersection);
    
        const double TIR_EPS = 1e-4;
        if (F >= 1.0 - TIR_EPS) {
            color += reflectedColor;    // pure reflection
            return color;
        }
    
        bool totalInternalReflection = false;
        Ray refractedRay = refract(
            ray,
            N,
            n1,
            n2,
            intersection.point,
            scene.shadowRayEpsilon,
            totalInternalReflection
        );
    
        if (totalInternalReflection) {
            // Numerical safety net — should be rare if we handled F ~ 1 above
            color += reflectedColor;
            return color;
        }
    
        if (material->roughness > 0.0) {
            refractedRay.direction = perturbDirection(refractedRay.direction,
                                                      material->roughness,
                                                      ray.random1,
                                                      ray.random2);
        }
    
        Intersection refractionIntersection = intersect(scene, refractedRay);
        VectorFloatTriplet refractedColor =
            computePixelColor(scene, refractedRay, refractionIntersection);
    
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
    
        // Fresnel mix between reflection and refraction
        color += F * reflectedColor + (1.0 - F) * refractedColor;
        return color;
    }

    if (material->isMirror) {
        Ray reflectedRay = reflect(ray,
                                   normal,
                                   intersection.point,
                                   scene.shadowRayEpsilon);
        
        if (material->roughness > 0.0) {
            reflectedRay.direction = perturbDirection(reflectedRay.direction,
                                                      material->roughness,
                                                      ray.random1,
                                                      ray.random2);
        }
        
        Intersection reflectionIntersection = intersect(scene, reflectedRay);
        VectorFloatTriplet mirrorColor =
            computePixelColor(scene, reflectedRay, reflectionIntersection);
        color += material->mirrorReflectance * mirrorColor;
    }

    for (const PointLight& light : scene.pointLights) {
        if (isInShadow(scene, ray, light, intersection)) {
            continue;
        }

        VectorFloatTriplet lightVec = light.position - intersection.point;
        double distanceSq = dotProduct(lightVec, lightVec);
        double distance = std::sqrt(distanceSq);
        VectorFloatTriplet lightDir = lightVec * (1.0 / distance);

        double attenuation = 1.0 / distanceSq;

        double NdotL = std::max(0.0, dotProduct(normal, lightDir));
        if (NdotL <= 0.0) continue;

        VectorFloatTriplet diffuse =
            diffuseReflectance * light.intensity *
            (NdotL * attenuation);

        VectorFloatTriplet viewDir   = normalize(-ray.direction);
        VectorFloatTriplet halfVector = normalize(lightDir + viewDir);
        double NdotH = std::max(0.0, dotProduct(normal, halfVector));
        double specularFactor = std::pow(NdotH, material->phongExponent);

        VectorFloatTriplet specular =
            specularReflectance * light.intensity *
            (specularFactor * NdotL * attenuation);

        color += diffuse + specular;
    }

    for (const AreaLight& light : scene.areaLights) {
        VectorFloatTriplet lightNormal = normalize(light.normal);
        VectorFloatTriplet u, v;
        orthonormalBasis(lightNormal, u, v);

        double halfSize = light.size / 2.0;
        double ksi_1 = (ray.random1 * 2.0 - 1.0) * halfSize;
        double ksi_2 = (ray.random2 * 2.0 - 1.0) * halfSize;
        VectorFloatTriplet samplePoint = light.position + u * ksi_1 + v * ksi_2;

        VectorFloatTriplet toLight = samplePoint - intersection.point;
        double distanceSq = dotProduct(toLight, toLight);
        double distance = std::sqrt(distanceSq);
        VectorFloatTriplet lightDir = toLight * (1.0 / distance);

        double cosTheta = dotProduct(normal, lightDir);
        
        double cosThetaLight = std::abs(dotProduct(lightNormal, -lightDir));
        
        if (cosTheta <= 0.0) continue;

        VectorFloatTriplet offsetNormal = intersection.geometricNormal;
        if (dotProduct(offsetNormal, lightDir) < 0.0) {
            offsetNormal = -offsetNormal;
        }
        Ray shadowRay{
            intersection.point + scene.shadowRayEpsilon * offsetNormal,
            lightDir,
            0,
            true,
            false,
            false,
            ray.time,
            ray.random1,
            ray.random2
        };
        Intersection shadowHit = intersect(scene, shadowRay);
        if (shadowHit.hit && shadowHit.distance < distance - scene.shadowRayEpsilon) continue;

        double area = light.size * light.size;
        VectorFloatTriplet irradiance = light.radiance * (cosTheta * cosThetaLight * area / distanceSq);

        VectorFloatTriplet diffuse = diffuseReflectance * irradiance;

        VectorFloatTriplet viewDir = normalize(-ray.direction);
        VectorFloatTriplet halfVector = normalize(lightDir + viewDir);
        double specFactor = std::pow(std::max(0.0, dotProduct(normal, halfVector)), material->phongExponent);
        VectorFloatTriplet specular = specularReflectance * irradiance * specFactor;

        color += diffuse + specular;
    }

    // *** FIX: Handle LightMeshes (emissive mesh geometry) for default shading ***
    for (int i = 0; i < (int)scene.lightMeshes.size(); i++) {
        const LightMesh& lightMesh = scene.lightMeshes[i];
        
        // Sample a point on the light mesh
        VectorFloatTriplet lightNormal;
        double pdf;
        VectorFloatTriplet samplePoint = sampleLightMesh(lightMesh, scene.vertices, 
                                                          ray.random1, ray.random2, 0.5, 
                                                          lightNormal, pdf);
        
        VectorFloatTriplet toLight = samplePoint - intersection.point;
        double distanceSq = dotProduct(toLight, toLight);
        double distance = std::sqrt(distanceSq);
        if (distance < 1e-10) continue;
        
        VectorFloatTriplet lightDir = toLight * (1.0 / distance);
        
        double cosTheta = dotProduct(normal, lightDir);
        if (cosTheta <= 0.0) continue;
        
        double cosThetaLight = std::abs(dotProduct(lightNormal, -lightDir));
        if (cosThetaLight < 1e-10) continue;
        
        // Shadow test
        VectorFloatTriplet offsetNormal = intersection.geometricNormal;
        if (dotProduct(offsetNormal, lightDir) < 0.0) {
            offsetNormal = -offsetNormal;
        }
        Ray shadowRay{
            intersection.point + scene.shadowRayEpsilon * offsetNormal,
            lightDir,
            0,
            true,
            false,
            false,
            ray.time,
            ray.random1,
            ray.random2
        };
        Intersection shadowHit = intersect(scene, shadowRay);
        if (shadowHit.hit && shadowHit.distance < distance - scene.shadowRayEpsilon) continue;
        
        // Compute irradiance: L * cos_theta_surface * cos_theta_light * area / r^2
        double totalArea = lightMesh.totalArea > 0.0 ? lightMesh.totalArea : 1.0;
        VectorFloatTriplet irradiance = lightMesh.radiance * (cosTheta * cosThetaLight * totalArea / distanceSq);
        
        VectorFloatTriplet diffuse = diffuseReflectance * irradiance;
        
        VectorFloatTriplet viewDir = normalize(-ray.direction);
        VectorFloatTriplet halfVector = normalize(lightDir + viewDir);
        double specFactor = std::pow(std::max(0.0, dotProduct(normal, halfVector)), material->phongExponent);
        VectorFloatTriplet specular = specularReflectance * irradiance * specFactor;
        
        color += diffuse + specular;
    }

    // *** FIX: Handle LightSpheres (emissive sphere geometry) for default shading ***
    for (int i = 0; i < (int)scene.lightSpheres.size(); i++) {
        const LightSphere& lightSphere = scene.lightSpheres[i];
        
        // Sample a point on the light sphere
        VectorFloatTriplet lightNormal;
        double pdf;
        VectorFloatTriplet samplePoint = sampleLightSphere(lightSphere, scene.vertices, 
                                                            intersection.point,
                                                            ray.random1, ray.random2, 
                                                            lightNormal, pdf);
        
        VectorFloatTriplet toLight = samplePoint - intersection.point;
        double distanceSq = dotProduct(toLight, toLight);
        double distance = std::sqrt(distanceSq);
        if (distance < 1e-10) continue;
        
        VectorFloatTriplet lightDir = toLight * (1.0 / distance);
        
        double cosTheta = dotProduct(normal, lightDir);
        if (cosTheta <= 0.0) continue;
        
        double cosThetaLight = std::abs(dotProduct(lightNormal, -lightDir));
        if (cosThetaLight < 1e-10) continue;
        
        // Shadow test
        VectorFloatTriplet offsetNormal = intersection.geometricNormal;
        if (dotProduct(offsetNormal, lightDir) < 0.0) {
            offsetNormal = -offsetNormal;
        }
        Ray shadowRay{
            intersection.point + scene.shadowRayEpsilon * offsetNormal,
            lightDir,
            0,
            true,
            false,
            false,
            ray.time,
            ray.random1,
            ray.random2
        };
        Intersection shadowHit = intersect(scene, shadowRay);
        if (shadowHit.hit && shadowHit.distance < distance - scene.shadowRayEpsilon) continue;
        
        // Compute irradiance
        double sphereArea = 4.0 * M_PI * lightSphere.radius * lightSphere.radius;
        VectorFloatTriplet irradiance = lightSphere.radiance * (cosTheta * cosThetaLight * sphereArea / distanceSq);
        
        VectorFloatTriplet diffuse = diffuseReflectance * irradiance;
        
        VectorFloatTriplet viewDir = normalize(-ray.direction);
        VectorFloatTriplet halfVector = normalize(lightDir + viewDir);
        double specFactor = std::pow(std::max(0.0, dotProduct(normal, halfVector)), material->phongExponent);
        VectorFloatTriplet specular = specularReflectance * irradiance * specFactor;
        
        color += diffuse + specular;
    }

    // Directional lights
    for (const DirectionalLight& light : scene.directionalLights) {
        VectorFloatTriplet lightDir = normalize(light.direction);
        
        // Check if surface is facing the light
        double NdotL = std::max(0.0, dotProduct(normal, lightDir));
        if (NdotL <= 0.0) continue;
        
        // Shadow test
        VectorFloatTriplet offsetNormal = intersection.geometricNormal;
        if (dotProduct(offsetNormal, lightDir) < 0.0) {
            offsetNormal = -offsetNormal;
        }
        Ray shadowRay{
            intersection.point + scene.shadowRayEpsilon * offsetNormal,
            lightDir,
            ray.depth + 1,
            true,
            false,
            false,
            ray.time
        };
        Intersection shadowHit = intersect(scene, shadowRay);
        if (shadowHit.hit) continue;  // In shadow
        
        // No distance attenuation for directional lights
        VectorFloatTriplet diffuse = diffuseReflectance * light.radiance * NdotL;
        
        VectorFloatTriplet viewDir = normalize(-ray.direction);
        VectorFloatTriplet halfVector = normalize(lightDir + viewDir);
        double NdotH = std::max(0.0, dotProduct(normal, halfVector));
        double specularFactor = std::pow(NdotH, material->phongExponent);
        VectorFloatTriplet specular = specularReflectance * light.radiance * specularFactor;
        
        color += diffuse + specular;
    }

    // Spot lights
    for (const SpotLight& light : scene.spotLights) {
        VectorFloatTriplet toLight = light.position - intersection.point;
        double distanceSq = dotProduct(toLight, toLight);
        double distance = std::sqrt(distanceSq);
        VectorFloatTriplet lightDir = toLight * (1.0 / distance);
        VectorFloatTriplet lightDirection = normalize(light.direction);
        
        // Calculate angle between light direction and vector to point
        double cosAlpha = dotProduct(-lightDir, lightDirection);
        double alpha = std::acos(std::max(-1.0, std::min(1.0, cosAlpha)));
        
        // Convert angles to radians
        double coverageRad = light.coverageAngle * M_PI / 180.0;
        double falloffRad = light.falloffAngle * M_PI / 180.0;
        
        // Check if point is within coverage angle
        if (alpha > coverageRad / 2.0) continue;  // Outside coverage cone
        
        // Calculate spot attenuation
        double spotAttenuation = 1.0;
        if (alpha > falloffRad / 2.0) {
            // Between falloff and coverage - apply attenuation
            double cosFalloff = std::cos(falloffRad / 2.0);
            double cosCoverage = std::cos(coverageRad / 2.0);
            double cosAlphaVal = std::cos(alpha);
            spotAttenuation = std::pow((cosAlphaVal - cosCoverage) / (cosFalloff - cosCoverage), 4.0);
        }
        
        // Check if surface is facing the light
        double NdotL = std::max(0.0, dotProduct(normal, lightDir));
        if (NdotL <= 0.0) continue;
        
        // Shadow test
        VectorFloatTriplet offsetNormal = intersection.geometricNormal;
        if (dotProduct(offsetNormal, lightDir) < 0.0) {
            offsetNormal = -offsetNormal;
        }
        Ray shadowRay{
            intersection.point + scene.shadowRayEpsilon * offsetNormal,
            lightDir,
            ray.depth + 1,
            true,
            false,
            false,
            ray.time
        };
        Intersection shadowHit = intersect(scene, shadowRay);
        if (shadowHit.hit && shadowHit.distance < distance - scene.shadowRayEpsilon) continue;
        
        // Distance-based attenuation
        double attenuation = 1.0 / distanceSq;
        
        VectorFloatTriplet diffuse = diffuseReflectance * light.intensity * (NdotL * attenuation * spotAttenuation);
        
        VectorFloatTriplet viewDir = normalize(-ray.direction);
        VectorFloatTriplet halfVector = normalize(lightDir + viewDir);
        double NdotH = std::max(0.0, dotProduct(normal, halfVector));
        double specularFactor = std::pow(NdotH, material->phongExponent);
        VectorFloatTriplet specular = specularReflectance * light.intensity * (specularFactor * NdotL * attenuation * spotAttenuation);
        
        color += diffuse + specular;
    }

    // Environment lights (Spherical Directional Lights)
    // Only apply to direct illumination (depth 0) to avoid double-counting and reduce noise
    // Environment lights provide indirect lighting, so they should only be sampled for primary rays
    if (ray.depth == 0) {
        for (const SphericalDirectionalLight& light : scene.sphericalDirectionalLights) {
        const Image* envImage = scene.getImageById(light.imageId);
        if (!envImage || !envImage->isHDR || !envImage->hdrData) {
            continue;
        }
        
        // Sample direction from hemisphere
        VectorFloatTriplet direction;
        double pdf;
        bool isCosineWeighted = false;
        
        if (light.sampler == "cosine") {
            // Cosine-weighted hemisphere sampling
            isCosineWeighted = true;
            double r1 = ray.random1;
            double r2 = ray.random2;
            double cosTheta = std::sqrt(r1);
            double sinTheta = std::sqrt(1.0 - r1);
            double phi = 2.0 * M_PI * r2;
            
            // Sample in local coordinate system (upper hemisphere)
            double x = sinTheta * std::cos(phi);
            double y = cosTheta;  // y is up (cosine-weighted)
            double z = sinTheta * std::sin(phi);
            
            // Transform to world space using surface normal
            VectorFloatTriplet u, v;
            orthonormalBasis(normal, u, v);
            direction = u * x + normal * y + v * z;
            direction = normalize(direction);
            
            pdf = cosTheta / M_PI;  // cos(theta) / pi
        } else {
            // Uniform hemisphere sampling
            isCosineWeighted = false;
            double r1 = ray.random1;
            double r2 = ray.random2;
            double cosTheta = r1;
            double sinTheta = std::sqrt(1.0 - r1 * r1);
            double phi = 2.0 * M_PI * r2;
            
            double x = sinTheta * std::cos(phi);
            double y = cosTheta;
            double z = sinTheta * std::sin(phi);
            
            VectorFloatTriplet u, v;
            orthonormalBasis(normal, u, v);
            direction = u * x + normal * y + v * z;
            direction = normalize(direction);
            
            pdf = 1.0 / (2.0 * M_PI);
        }
        
        // Convert direction to UV coordinates
        double u, v;
        if (light.type == "latlong") {
            // Latitude-longitude (equirectangular) mapping
            double clampedY = std::max(-1.0, std::min(1.0, direction.y));
            // u wraps around Y-axis  
            u = 0.5 + std::atan2(direction.x, -direction.z) / (2.0 * M_PI);
            v = std::acos(clampedY) / M_PI;
        } else {
            // Spherical (probe) mapping
            double denom = std::sqrt(direction.x * direction.x + direction.y * direction.y);
            double r = (1.0 / M_PI) * std::acos(-direction.z);
            if (denom > 1e-10) {
                r = r / denom;
            } else {
                r = 0.0;  // Handle division by zero
            }
            u = (r * direction.x + 1.0) / 2.0;
            v = (-r * direction.y + 1.0) / 2.0;
        }
        
        // Clamp UV to [0,1]
        u = std::max(0.0, std::min(1.0, u));
        v = std::max(0.0, std::min(1.0, v));
        
        // Sample from HDR environment map using bilinear interpolation
        double x = u * (envImage->width - 1);
        double y = v * (envImage->height - 1);
        int x0 = (int)floor(x);
        int y0 = (int)floor(y);
        int x1 = std::min(x0 + 1, envImage->width - 1);
        int y1 = std::min(y0 + 1, envImage->height - 1);
        
        double fx = x - x0;
        double fy = y - y0;
        
        // Bounds checking - ensure indices are valid and account for channel access
        int maxIdx = envImage->width * envImage->height * envImage->channels;
        int idx00 = (y0 * envImage->width + x0) * envImage->channels;
        int idx10 = (y0 * envImage->width + x1) * envImage->channels;
        int idx01 = (y1 * envImage->width + x0) * envImage->channels;
        int idx11 = (y1 * envImage->width + x1) * envImage->channels;
        
        // Safety check - skip if indices are out of bounds (account for channel access)
        int maxChannelOffset = (envImage->channels >= 3) ? 2 : 0;
        if (idx00 < 0 || (idx00 + maxChannelOffset) >= maxIdx ||
            idx10 < 0 || (idx10 + maxChannelOffset) >= maxIdx ||
            idx01 < 0 || (idx01 + maxChannelOffset) >= maxIdx ||
            idx11 < 0 || (idx11 + maxChannelOffset) >= maxIdx) {
            continue;  // Skip this environment light sample
        }
        
        VectorFloatTriplet c00, c10, c01, c11;
        if (envImage->channels >= 3) {
            c00 = VectorFloatTriplet{envImage->hdrData[idx00], envImage->hdrData[idx00 + 1], envImage->hdrData[idx00 + 2]};
            c10 = VectorFloatTriplet{envImage->hdrData[idx10], envImage->hdrData[idx10 + 1], envImage->hdrData[idx10 + 2]};
            c01 = VectorFloatTriplet{envImage->hdrData[idx01], envImage->hdrData[idx01 + 1], envImage->hdrData[idx01 + 2]};
            c11 = VectorFloatTriplet{envImage->hdrData[idx11], envImage->hdrData[idx11 + 1], envImage->hdrData[idx11 + 2]};
        } else {
            double g00 = envImage->hdrData[idx00];
            double g10 = envImage->hdrData[idx10];
            double g01 = envImage->hdrData[idx01];
            double g11 = envImage->hdrData[idx11];
            c00 = VectorFloatTriplet{g00, g00, g00};
            c10 = VectorFloatTriplet{g10, g10, g10};
            c01 = VectorFloatTriplet{g01, g01, g01};
            c11 = VectorFloatTriplet{g11, g11, g11};
        }
        
        VectorFloatTriplet c0 = c00 * (1.0 - fx) + c10 * fx;
        VectorFloatTriplet c1 = c01 * (1.0 - fx) + c11 * fx;
        VectorFloatTriplet radiance = c0 * (1.0 - fy) + c1 * fy;
        
        // Check for NaN or Inf values
        if (std::isnan(radiance.x) || std::isnan(radiance.y) || std::isnan(radiance.z) ||
            std::isinf(radiance.x) || std::isinf(radiance.y) || std::isinf(radiance.z)) {
            continue;  // Skip invalid radiance values
        }
        
        // Compute cosine term (NdotL) for material interaction
        double NdotL = std::max(0.0, dotProduct(normal, direction));
        if (NdotL <= 0.0) continue;  // Skip if direction is below surface
        
        // Multiply by material properties
        // Monte Carlo estimator: (L_i(ω) * f_r(ω) * cos(θ)) / PDF(ω)
        // Based on other lights (diffuseReflectance * light.radiance * NdotL), BRDF is f_r = ρ_d (not ρ_d/π)
        // For cosine-weighted sampling: PDF = cos(θ)/π
        // Estimator = (L_i * ρ_d * cos(θ)) / (cos(θ)/π) = L_i * ρ_d * π
        // For uniform sampling: PDF = 1/(2π)
        // Estimator = (L_i * ρ_d * cos(θ)) / (1/(2π)) = L_i * ρ_d * cos(θ) * 2π
        VectorFloatTriplet envContribution;
        if (isCosineWeighted) {
            // Cosine-weighted sampling with physically correct Lambertian BRDF (ρ/π):
            // PDF = cos(θ)/π, BRDF = ρ/π, Integrand = L * (ρ/π) * cos(θ)
            // Estimator = Integrand/PDF = L * (ρ/π) * cos(θ) / (cos(θ)/π) = L * ρ
            // The π terms cancel out, so no π multiplier needed
            envContribution = radiance * diffuseReflectance;
        } else {
            // Uniform hemisphere sampling with physically correct Lambertian BRDF:
            // PDF = 1/(2π), BRDF = ρ/π, Integrand = L * (ρ/π) * cos(θ)
            // Estimator = Integrand/PDF = L * (ρ/π) * cos(θ) * 2π = L * ρ * 2 * cos(θ)
            envContribution = radiance * diffuseReflectance * NdotL * 2.0;
        }
        
        color += envContribution;
        }  // End of for loop over environment lights
    }  // End of if (ray.depth == 0)

    return color;
}

// Helper function to check if an intersection is an emissive object (LightSphere or LightMesh)
static bool isEmissiveObject(const Scene& scene, const Intersection& intersection) {
    return (intersection.kind == Intersection::Kind::LightSphere || 
            intersection.kind == Intersection::Kind::LightMesh);
}

// Helper function to get emission radiance from an emissive object
static VectorFloatTriplet getEmission(const Scene& scene, const Intersection& intersection) {
    if (intersection.kind == Intersection::Kind::LightSphere && 
        intersection.containerIndex >= 0 && 
        intersection.containerIndex < (int)scene.lightSpheres.size()) {
        return scene.lightSpheres[intersection.containerIndex].radiance;
    } else if (intersection.kind == Intersection::Kind::LightMesh && 
               intersection.containerIndex >= 0 && 
               intersection.containerIndex < (int)scene.lightMeshes.size()) {
        return scene.lightMeshes[intersection.containerIndex].radiance;
    }
    return VectorFloatTriplet{0.0, 0.0, 0.0};
}

VectorFloatTriplet computePathTracing(const Scene& scene, const Camera& camera, Ray& ray, const Intersection& intersection) {
    VectorFloatTriplet L = {0.0, 0.0, 0.0};
    VectorFloatTriplet throughput = {1.0, 1.0, 1.0};
    
    Ray currentRay = ray;
    Intersection currentHit = intersection;
    
    // Use camera's maxRecursionDepth if set, otherwise use scene's
    int maxDepth = (camera.maxRecursionDepth > 0) ? camera.maxRecursionDepth : scene.maxRecursionDepth;
    
    // Splitting factor: for primary rays, send multiple indirect rays
    int splittingFactor = (ray.depth == 0) ? camera.splittingFactor : 1;
    
    for (int depth = 0; depth < maxDepth; depth++) {
        // If we have a hit from previous iteration, use it; otherwise intersect
        if (depth > 0 || !currentHit.hit) {
            currentHit = intersect(scene, currentRay);
        }
        
        if (!currentHit.hit) {
            // Ray hit background - add background color weighted by throughput
            L = L + throughput * scene.backgroundColor;
            break;
        }
        
        // Check if we hit an emissive object (light source)
        if (isEmissiveObject(scene, currentHit)) {
            VectorFloatTriplet emission = getEmission(scene, currentHit);
            L = L + throughput * emission;
            // For pure path tracing, we can terminate here or continue
            // Typically we continue to allow indirect lighting
        }
        
        // Get material and BRDF
        Material* material = currentHit.material;
        if (!material) {
            break;  // No material, terminate path
        }
        
        const BRDF* brdf = nullptr;
        if (material->brdfId > 0) {
            brdf = scene.getBRDFById(material->brdfId);
        }
        
        // Sample next direction based on importance sampling setting
        VectorFloatTriplet N = normalize(currentHit.shadingNormal);
        VectorFloatTriplet wo = normalize(-currentRay.direction);  // Outgoing direction (to camera)
        // Ensure normal faces outgoing direction to keep BRDF terms valid
        if (dotProduct(N, wo) < 0.0) {
            N = -N;
        }
        
        // Generate random numbers for this bounce
        // Use a simple hash-based approach to generate new random numbers from depth and original random values
        double xi1, xi2;
        if (depth == 0) {
            // First bounce: use precomputed random values
            xi1 = currentRay.random1;
            xi2 = currentRay.random2;
        } else {
            // Subsequent bounces: generate new random numbers using a simple hash
            // This is a simple LCG-like approach
            double seed1 = currentRay.random1 * 1000.0 + depth * 17.0;
            double seed2 = currentRay.random2 * 1000.0 + depth * 23.0;
            xi1 = std::fmod(seed1 * 1103515245.0 + 12345.0, 2147483648.0) / 2147483648.0;
            xi2 = std::fmod(seed2 * 1103515245.0 + 12345.0, 2147483648.0) / 2147483648.0;
            // Ensure in [0,1)
            xi1 = std::max(0.0, std::min(0.999999, xi1));
            xi2 = std::max(0.0, std::min(0.999999, xi2));
        }
        
        // Next Event Estimation: sample light directly (only once, not per split)
        VectorFloatTriplet directContribution = {0.0, 0.0, 0.0};
        if (camera.nextEventEstimation) {
            double pdfLight = 0.0;
            VectorFloatTriplet lightDir;
            // Generate additional random numbers for light sampling
            double xi3, xi4;
            if (depth == 0) {
                // Reuse random values - but we need 4 random numbers, so generate more
                xi3 = std::fmod((currentRay.random1 + 0.5) * 1000.0, 1.0);
                xi4 = std::fmod((currentRay.random2 + 0.5) * 1000.0, 1.0);
            } else {
                double seed3 = currentRay.random1 * 1000.0 + depth * 31.0;
                double seed4 = currentRay.random2 * 1000.0 + depth * 37.0;
                xi3 = std::fmod(seed3 * 1103515245.0 + 12345.0, 2147483648.0) / 2147483648.0;
                xi4 = std::fmod(seed4 * 1103515245.0 + 12345.0, 2147483648.0) / 2147483648.0;
                xi3 = std::max(0.0, std::min(0.999999, xi3));
                xi4 = std::max(0.0, std::min(0.999999, xi4));
            }
            
            VectorFloatTriplet lightRadiance = sampleDirectLight(scene, currentHit.point, N, 
                                                                   xi1, xi2, xi3, xi4, lightDir, pdfLight);
            
            if (pdfLight > 1e-10) {
                // Evaluate BRDF for light direction
                VectorFloatTriplet diffuseLight, specularLight;
                VectorFloatTriplet brdfLight = brdf::evaluateBRDF(*material, brdf, N, lightDir, wo, diffuseLight, specularLight);
                double cosThetaLight = std::max(0.0, dotProduct(N, lightDir));
                
                // Compute BRDF PDF for this light direction
                double pdfBRDFLight = brdf::pdfBRDF(brdf, N, lightDir, wo);
                
                // MIS weight for light sample
                double wLight = 1.0;
                if (!camera.misHeuristic.empty()) {
                    wLight = misWeight(pdfLight, pdfBRDFLight, camera.misHeuristic);
                }
                
                // Direct light contribution: w_light * (L * BRDF * cos) / pdf_light
                if (cosThetaLight > 1e-10) {
                    VectorFloatTriplet brdfCosLight = brdfLight * cosThetaLight;
                    double invPdfLight = 1.0 / pdfLight;
                    directContribution = lightRadiance * brdfCosLight * invPdfLight * wLight;
                }
            }
        }
        
        // Splitting: for first bounce from primary ray, accumulate multiple indirect contributions
        VectorFloatTriplet indirectAccum = {0.0, 0.0, 0.0};
        int numSplits = (depth == 0 && ray.depth == 0) ? splittingFactor : 1;
        VectorFloatTriplet offsetPoint = currentHit.point + scene.shadowRayEpsilon * N;
        
        for (int splitIdx = 0; splitIdx < numSplits; splitIdx++) {
            // Generate random numbers for this split (vary for splitting)
            double xi1_split = xi1, xi2_split = xi2;
            if (splitIdx > 0) {
                // For additional splits, generate different random numbers
                double seed1 = currentRay.random1 * 1000.0 + splitIdx * 17.0;
                double seed2 = currentRay.random2 * 1000.0 + splitIdx * 23.0;
                xi1_split = std::fmod(seed1 * 1103515245.0 + 12345.0, 2147483648.0) / 2147483648.0;
                xi2_split = std::fmod(seed2 * 1103515245.0 + 12345.0, 2147483648.0) / 2147483648.0;
                xi1_split = std::max(0.0, std::min(0.999999, xi1_split));
                xi2_split = std::max(0.0, std::min(0.999999, xi2_split));
            }
        
            // Sample BRDF direction (indirect illumination) for this split
            VectorFloatTriplet wi_split;  // Incoming direction (from light)
            double pdfBRDF_split;
            
            if (camera.importanceSampling) {
                // Cosine-weighted hemisphere sampling
                wi_split = sampleHemisphereCosine(N, xi1_split, xi2_split);
                double cosTheta_split = std::max(0.0, dotProduct(N, wi_split));
                pdfBRDF_split = cosTheta_split / M_PI;
            } else {
                // Uniform hemisphere sampling
                wi_split = sampleHemisphereUniform(N, xi1_split, xi2_split);
                pdfBRDF_split = 1.0 / (2.0 * M_PI);
            }
            
            // Evaluate BRDF for sampled direction
            VectorFloatTriplet diffuse_split, specular_split;
            VectorFloatTriplet brdfValue_split = brdf::evaluateBRDF(*material, brdf, N, wi_split, wo, diffuse_split, specular_split);
            
            // Compute cosine term
            double cosTheta_split = std::max(0.0, dotProduct(N, wi_split));
            
            // Check if BRDF-sampled direction hits a light
            VectorFloatTriplet indirectContribution_split = {0.0, 0.0, 0.0};
            Ray testRay(offsetPoint, wi_split, currentRay.depth + 1, false, false, false, currentRay.time, 0.0, 0.0);
            Intersection testHit = intersect(scene, testRay);
            
            if (testHit.hit && isEmissiveObject(scene, testHit)) {
                // BRDF direction hit a light - compute contribution with MIS
                VectorFloatTriplet emission = getEmission(scene, testHit);
                
                // Compute light PDF for this direction
                double pdfLightBRDF = 0.0;
                if (testHit.kind == Intersection::Kind::LightSphere) {
                    const LightSphere& light = scene.lightSpheres[testHit.containerIndex];
                    VectorFloatTriplet toLight = testHit.point - currentHit.point;
                    double distanceSq = dotProduct(toLight, toLight);
                    double distance = std::sqrt(distanceSq);
                    double sphereArea = 4.0 * M_PI * light.radius * light.radius;
                    int totalLights = scene.pointLights.size() + scene.areaLights.size() + 
                                      scene.lightSpheres.size() + scene.lightMeshes.size();
                    if (sphereArea > 1e-10 && totalLights > 0) {
                        double areaPDF = 1.0 / sphereArea / totalLights;
                        double cosAtLight = 1.0;  // Approximation
                        pdfLightBRDF = areaPDFToSolidAnglePDF(areaPDF, distance, cosAtLight);
                    }
                } else if (testHit.kind == Intersection::Kind::LightMesh) {
                    const LightMesh& light = scene.lightMeshes[testHit.containerIndex];
                    VectorFloatTriplet toLight = testHit.point - currentHit.point;
                    double distanceSq = dotProduct(toLight, toLight);
                    double distance = std::sqrt(distanceSq);
                    int totalLights = scene.pointLights.size() + scene.areaLights.size() + 
                                      scene.lightSpheres.size() + scene.lightMeshes.size();
                    if (light.totalArea > 1e-10 && totalLights > 0) {
                        double areaPDF = 1.0 / light.totalArea / totalLights;
                        double cosAtLight = 1.0;  // Approximation
                        pdfLightBRDF = areaPDFToSolidAnglePDF(areaPDF, distance, cosAtLight);
                    }
                }
                
                // MIS weight for BRDF sample
                double wBRDF = 1.0;
                if (!camera.misHeuristic.empty() && pdfLightBRDF > 1e-10) {
                    wBRDF = misWeight(pdfBRDF_split, pdfLightBRDF, camera.misHeuristic);
                }
                
                // Indirect light contribution: w_brdf * (L * BRDF * cos) / pdf_brdf
                if (cosTheta_split > 1e-10 && pdfBRDF_split > 1e-10) {
                    VectorFloatTriplet brdfCos = brdfValue_split * cosTheta_split;
                    double invPdfBRDF = 1.0 / pdfBRDF_split;
                    indirectContribution_split = emission * brdfCos * invPdfBRDF * wBRDF;
                }
            }
            
            // Accumulate indirect contribution for this split
            indirectAccum = indirectAccum + indirectContribution_split;
        }  // End of splitting loop
        
        // Average indirect contributions from splitting
        if (numSplits > 1) {
            double invSplits = 1.0 / numSplits;
            indirectAccum = indirectAccum * invSplits;
        }
        
        // Use first split's direction for path continuation
        VectorFloatTriplet wi;  // Will be set below
        double pdfBRDF;
        VectorFloatTriplet brdfValue;
        double cosTheta;
        
        // Sample direction for continuation (using first random numbers)
        if (camera.importanceSampling) {
            wi = sampleHemisphereCosine(N, xi1, xi2);
            cosTheta = std::max(0.0, dotProduct(N, wi));
            pdfBRDF = cosTheta / M_PI;
        } else {
            wi = sampleHemisphereUniform(N, xi1, xi2);
            cosTheta = std::max(0.0, dotProduct(N, wi));
            pdfBRDF = 1.0 / (2.0 * M_PI);
        }
        
        // Evaluate BRDF for continuation direction
        VectorFloatTriplet diffuse, specular;
        brdfValue = brdf::evaluateBRDF(*material, brdf, N, wi, wo, diffuse, specular);
        
        // Add contributions (direct + averaged indirect from splitting)
        L = L + throughput * (directContribution + indirectAccum);
        
        // Update throughput: throughput *= BRDF * cos(theta) / PDF
        if (pdfBRDF > 1e-10 && cosTheta > 1e-10) {
            VectorFloatTriplet brdfCos = brdfValue * cosTheta;
            double invPdf = 1.0 / pdfBRDF;
            throughput = throughput * brdfCos * invPdf;
        } else {
            break;  // Zero PDF, terminate path
        }
        
        // Check for NaN/Inf in throughput
        if (std::isnan(throughput.x) || std::isnan(throughput.y) || std::isnan(throughput.z) ||
            std::isinf(throughput.x) || std::isinf(throughput.y) || std::isinf(throughput.z)) {
            break;
        }
        
        // Russian Roulette: probabilistically terminate path after minRecursionDepth
        if (camera.russianRoulette && depth >= camera.minRecursionDepth) {
            // Compute survival probability based on throughput
            double throughputMax = std::max(throughput.x, std::max(throughput.y, throughput.z));
            double survivalProb = std::min(std::max(throughputMax, 0.01), 0.99);  // Clamp to [0.01, 0.99]
            
            // Generate random number for Russian Roulette
            double rrRandom;
            if (depth == 0) {
                rrRandom = std::fmod((currentRay.random1 + currentRay.random2) * 1000.0, 1.0);
            } else {
                double seedRR = currentRay.random1 * 1000.0 + depth * 41.0;
                rrRandom = std::fmod(seedRR * 1103515245.0 + 12345.0, 2147483648.0) / 2147483648.0;
                rrRandom = std::max(0.0, std::min(0.999999, rrRandom));
            }
            
            if (rrRandom > survivalProb) {
                // Terminate path
                break;
            }
            
            // Scale throughput by survival probability
            double invSurvivalProb = 1.0 / survivalProb;
            throughput = throughput * invSurvivalProb;
        }
        
        // Create new ray for next bounce
        currentRay = Ray(offsetPoint, wi, currentRay.depth + 1, false, false, false, 
                        currentRay.time, 0.0, 0.0);  // Will need new random values for next bounce
        
        // Reset hit for next iteration
        currentHit.hit = false;
    }
    
    return L;
}

VectorFloatTriplet computePixelColor(const Scene& scene, Ray& ray, const Intersection& intersection) {
    // Check if we should use path tracing
    if (!scene.cameras.empty() && scene.currentCameraIndex < scene.cameras.size()) {
        const Camera& camera = scene.cameras[scene.currentCameraIndex];
        if (camera.renderer == "PathTracing") {
            return computePathTracing(scene, camera, ray, intersection);
        }
    }
    
    // Default shading (original implementation)
    if (ray.depth > scene.maxRecursionDepth) {
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }

    if (intersection.hit) {
        // *** FIX: If we hit an emissive object (LightMesh/LightSphere), return emission ***
        if (isEmissiveObject(scene, intersection)) {
            return getEmission(scene, intersection);
        }
        return computeShading(scene, ray, intersection);
    }

    // For background rays (no hit), sample environment lights directly
    // Environment lights should contribute to background when rays don't hit anything
    if (!scene.sphericalDirectionalLights.empty()) {
        VectorFloatTriplet envColor{0.0, 0.0, 0.0};
        for (const SphericalDirectionalLight& light : scene.sphericalDirectionalLights) {
            const Image* envImage = scene.getImageById(light.imageId);
            if (!envImage || !envImage->isHDR || !envImage->hdrData) continue;
            
            // Convert ray direction to UV coordinates for environment map
            VectorFloatTriplet direction = normalize(ray.direction);
            double u, v;
            
            if (light.type == "latlong") {
                // Latitude-longitude (equirectangular) mapping
                double clampedY = std::max(-1.0, std::min(1.0, direction.y));
                // u wraps around Y-axis
                u = 0.5 + std::atan2(direction.x, -direction.z) / (2.0 * M_PI);
                v = std::acos(clampedY) / M_PI;
            } else {
                // Spherical (probe) mapping
                double denom = std::sqrt(direction.x * direction.x + direction.y * direction.y);
                double r = (1.0 / M_PI) * std::acos(-direction.z);
                if (denom > 1e-10) {
                    r = r / denom;
                } else {
                    r = 0.0;
                }
                u = (r * direction.x + 1.0) / 2.0;
                v = (-r * direction.y + 1.0) / 2.0;
            }
            
            // Clamp UV to [0,1]
            u = std::max(0.0, std::min(1.0, u));
            v = std::max(0.0, std::min(1.0, v));
            
            // Sample from HDR environment map using bilinear interpolation
            double x = u * (envImage->width - 1);
            double y = v * (envImage->height - 1);
            int x0 = (int)floor(x);
            int y0 = (int)floor(y);
            int x1 = std::min(x0 + 1, envImage->width - 1);
            int y1 = std::min(y0 + 1, envImage->height - 1);
            
            double fx = x - x0;
            double fy = y - y0;
            
            int maxIdx = envImage->width * envImage->height * envImage->channels;
            int idx00 = (y0 * envImage->width + x0) * envImage->channels;
            int idx10 = (y0 * envImage->width + x1) * envImage->channels;
            int idx01 = (y1 * envImage->width + x0) * envImage->channels;
            int idx11 = (y1 * envImage->width + x1) * envImage->channels;
            
            // Bounds check
            int maxChannelOffset = (envImage->channels >= 3) ? 2 : 0;
            if (idx00 < 0 || (idx00 + maxChannelOffset) >= maxIdx ||
                idx10 < 0 || (idx10 + maxChannelOffset) >= maxIdx ||
                idx01 < 0 || (idx01 + maxChannelOffset) >= maxIdx ||
                idx11 < 0 || (idx11 + maxChannelOffset) >= maxIdx) {
                continue;
            }
            
            VectorFloatTriplet c00, c10, c01, c11;
            if (envImage->channels >= 3) {
                c00 = VectorFloatTriplet{envImage->hdrData[idx00], envImage->hdrData[idx00 + 1], envImage->hdrData[idx00 + 2]};
                c10 = VectorFloatTriplet{envImage->hdrData[idx10], envImage->hdrData[idx10 + 1], envImage->hdrData[idx10 + 2]};
                c01 = VectorFloatTriplet{envImage->hdrData[idx01], envImage->hdrData[idx01 + 1], envImage->hdrData[idx01 + 2]};
                c11 = VectorFloatTriplet{envImage->hdrData[idx11], envImage->hdrData[idx11 + 1], envImage->hdrData[idx11 + 2]};
            } else {
                double g00 = envImage->hdrData[idx00];
                double g10 = envImage->hdrData[idx10];
                double g01 = envImage->hdrData[idx01];
                double g11 = envImage->hdrData[idx11];
                c00 = VectorFloatTriplet{g00, g00, g00};
                c10 = VectorFloatTriplet{g10, g10, g10};
                c01 = VectorFloatTriplet{g01, g01, g01};
                c11 = VectorFloatTriplet{g11, g11, g11};
            }
            
            VectorFloatTriplet c0 = c00 * (1.0 - fx) + c10 * fx;
            VectorFloatTriplet c1 = c01 * (1.0 - fx) + c11 * fx;
            VectorFloatTriplet radiance = c0 * (1.0 - fy) + c1 * fy;
            
            // Check for NaN or Inf values
            if (std::isnan(radiance.x) || std::isnan(radiance.y) || std::isnan(radiance.z) ||
                std::isinf(radiance.x) || std::isinf(radiance.y) || std::isinf(radiance.z)) {
                continue;
            }
            
            envColor += radiance;
        }
        
        // If we have environment light contribution, return it (environment maps are the background)
        if (envColor.x > 0.0 || envColor.y > 0.0 || envColor.z > 0.0) {
            return envColor;
        }
    }

    if (scene.backgroundTextureId != 0) {
        const TextureMap* bgTex = scene.getTextureMapById(scene.backgroundTextureId);
        if (bgTex && !scene.cameras.empty()) {
            const Camera& cam = scene.cameras[scene.currentCameraIndex];
            
            // Build camera coordinate system
            VectorFloatTriplet w = -normalize(cam.gaze);
            VectorFloatTriplet vCam = normalize(cam.up);
            VectorFloatTriplet uCam = crossProduct(vCam, w);
            
            VectorFloatTriplet d = normalize(ray.direction);
            VectorFloatTriplet gazeDir = normalize(cam.gaze);
            
            double denom = dotProduct(d, gazeDir);
            
            VectorFloatPair uv;
            if (std::fabs(denom) > 1e-9) {
                double numer = cam.nearDistance - dotProduct(ray.origin - cam.position, gazeDir);
                double t = numer / denom;
                
                VectorFloatTriplet hitPoint = ray.origin + d * t;
                
                VectorFloatTriplet relPoint = hitPoint - cam.position;
                
                double uCoord = dotProduct(relPoint, uCam);
                double vCoord = dotProduct(relPoint, vCam);
                
                double l = cam.nearPlane.x;
                double r = cam.nearPlane.y;
                double b = cam.nearPlane.z;
                double tTop = cam.nearPlane.w;
                
                // Convert to UV coordinates [0, 1]
                uv.x = (uCoord - l) / (r - l);
                uv.y = (tTop - vCoord) / (tTop - b);
                
                // Clamp to [0, 1] for rays that might be outside the image bounds
                uv.x = std::max(0.0, std::min(1.0, uv.x));
                uv.y = std::max(0.0, std::min(1.0, uv.y));
            } else {
                uv.x = 0.5;
                uv.y = 0.5;
            }
            
            // The renderer expects [0,255] here, so scale sampled color
            return sampleTexture(bgTex, uv, ray.origin, &scene, true) * 255.0;
        }
    }

    return scene.backgroundColor;
}

bool isInShadow(const Scene& scene, Ray& ray, const PointLight& light, const Intersection& intersection) {
    VectorFloatTriplet lightDir = normalize(light.position - intersection.point);
    // Use geometric normal and ensure it points in the same hemisphere as the light
    VectorFloatTriplet offsetNormal = intersection.geometricNormal;
    if (dotProduct(offsetNormal, lightDir) < 0.0) {
        offsetNormal = -offsetNormal;
    }
    Ray shadowRay{
        intersection.point + scene.shadowRayEpsilon * offsetNormal,
        lightDir,
        0,
        true,   // shadow ray
        false,  // reflection ray
        false,  // refraction ray
        ray.time,  // Propagate time for motion blur consistency
        ray.random1,
        ray.random2
    };

    Intersection shadowIntersection = intersect(scene, shadowRay);
    double distToLight = std::sqrt(
        dotProduct(light.position - intersection.point,
                   light.position - intersection.point)
    );

    return shadowIntersection.hit && shadowIntersection.distance < distToLight;
}

double fresnelConductor(double cosTheta, double n, double k) {
    cosTheta = std::max(0.0, std::min(1.0, cosTheta));
    double cosThetaSq = cosTheta * cosTheta;

    n = std::max(n, 0.0);
    k = std::max(k, 0.0);

    double n2 = n * n;
    double k2 = k * k;
    double n2PlusK2 = n2 + k2;

    double twoNCos = 2.0 * n * cosTheta;

    double RsNum = n2PlusK2 - twoNCos + cosThetaSq;
    double RsDen = n2PlusK2 + twoNCos + cosThetaSq;
    double Rs = RsNum / RsDen;

    double RpNum = n2PlusK2 * cosThetaSq - twoNCos + 1.0;
    double RpDen = n2PlusK2 * cosThetaSq + twoNCos + 1.0;
    double Rp = RpNum / RpDen;

    double R = 0.5 * (Rs + Rp);

    return R;
}

double fresnelDielectric(double cosThetaI, double n1, double n2) {
    cosThetaI = std::max(-1.0, std::min(1.0, cosThetaI));
    n1 = std::max(n1, 0.0);
    n2 = std::max(n2, 0.0);

    double absCosThetaI = std::fabs(cosThetaI);

    double eta = n1 / n2;
    double sinThetaTSq = eta * eta * (1.0 - absCosThetaI * absCosThetaI);

    if (sinThetaTSq >= 1.0) {
        return 1.0;
    }

    double cosThetaT = std::sqrt(std::max(0.0, 1.0 - sinThetaTSq));

    double rPerp = (n1 * absCosThetaI - n2 * cosThetaT) /
                   (n1 * absCosThetaI + n2 * cosThetaT);

    double rPara = (n2 * absCosThetaI - n1 * cosThetaT) /
                   (n2 * absCosThetaI + n1 * cosThetaT);

    double R = 0.5 * (rPerp * rPerp + rPara * rPara);
    return std::max(0.0, std::min(1.0, R));
}

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
        return Ray(point,
                   ray.direction,
                   ray.depth,
                   false,
                   false,
                   false,
                   ray.time,
                   ray.random1,
                   ray.random2);
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
        true,
        ray.time,          // Propagate time for motion blur consistency
        ray.random1,
        ray.random2
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
        false,
        ray.time,  // Propagate time for motion blur consistency
        ray.random1,
        ray.random2
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