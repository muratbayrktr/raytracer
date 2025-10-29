#include <math.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "scene.h"
#include "utils.h"
#include "overloads.h"
#include "bvh.h"

using namespace std;
using namespace scene;

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
    float l = camera.nearPlane.x;
    float r = camera.nearPlane.y;
    float b = camera.nearPlane.z;
    float t = camera.nearPlane.w;
    float s_u = (x+0.5)*(r - l) / width;
    float s_v = (y+0.5)*(t - b) / height;
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
    float& t_min, 
    Intersection& intersection,
    int planeIndex
) {
        /*
            From lecture slides:
            Formula is simple (origin + t * direction - a) * normal = 0
            Solving for t gives us t = (a - origin) * normal / (direction * normal)
        */
       VectorFloatTriplet n = plane.normal;
       VectorFloatTriplet a = vertices[plane.point];
       VectorFloatTriplet d = ray.direction;
       VectorFloatTriplet o = ray.origin;
       float denom = dotProduct(d, n);
       if(fabs(denom) < 1e-9) return false;
       float t = dotProduct(a - o, n) / denom;
       if(t > 0 && t < t_min) {
        t_min = t;
        intersection.hit = true;
        intersection.distance = t;
        intersection.point = o + d * t;
        intersection.geometricNormal = n;
        intersection.shadingNormal = n;
        intersection.material = plane.material;
        intersection.kind = Intersection::Kind::Plane;
        intersection.containerIndex = planeIndex;
        return true;
       }
       return false;
}

bool rayHitsSphere(
    Ray& ray, 
    const Sphere& sphere, 
    const vector<VectorFloatTriplet>& vertices, 
    float& t_min, 
    Intersection& intersection,
    int sphereIndex
) {
    /*
        From the lecture slides:
        After solving for t we get the following formula:
            t = [-d . (o - c)  +- sqrt((d . (o - c))^2 - (d . d) * ((o - c)^2 - r^2))] / (d . d)
    */
    VectorFloatTriplet c = vertices[sphere.center];
    VectorFloatTriplet o = ray.origin;
    VectorFloatTriplet d = ray.direction;
    float r = sphere.radius;
    float D = dotProduct(d, o - c) * dotProduct(d, o - c) 
                - (dotProduct(d, d) * (dotProduct(o - c, o - c) - r * r));
    if(D < 0) return false;
    float t1 = (-dotProduct(d, o - c) + sqrt(D)) / dotProduct(d, d);
    float t2 = (-dotProduct(d, o - c) - sqrt(D)) / dotProduct(d, d);
    float t = min(t1, t2);
    VectorFloatTriplet intersection_point = o + d * t;
    VectorFloatTriplet intersection_normal = normalize(intersection_point - c);
    if(t > 0 && t < t_min) {
        t_min = t;
        intersection.hit = true;
        intersection.distance = t;
        intersection.point = intersection_point;
        intersection.geometricNormal = intersection_normal;
        intersection.shadingNormal = intersection_normal;
        intersection.material = sphere.material;
        intersection.kind = Intersection::Kind::Sphere;
        intersection.containerIndex = sphereIndex;
        return true;
    }
    return false;
}


bool rayHitsTriangle(
    Ray& ray, 
    const VectorIntTriplet& face, 
    const vector<VectorFloatTriplet>& vertices, 
    float& t_min, 
    Intersection& intersection, 
    float intersectionTestEpsilon, 
    float determinantT, 
    Material* material,
    bool enableBackFaceCulling,
    int containerIndex,
    int faceIndex
) {
    const VectorFloatTriplet a = vertices[face.x];
    const VectorFloatTriplet b = vertices[face.y];
    const VectorFloatTriplet c = vertices[face.z];
    
    const VectorFloatTriplet e1 = b - a;
    const VectorFloatTriplet e2 = c - a;
    const VectorFloatTriplet geometricNormal = normalize(crossProduct(e1, e2));
    if (enableBackFaceCulling) {
        if (dotProduct(ray.direction, geometricNormal) > 0.0f) {
            return false;
        }
    }
    const float ax=a.x, ay=a.y, az=a.z, bx=b.x, by=b.y, bz=b.z, cx=c.x, cy=c.y, cz=c.z;
    const float ox=ray.origin.x, oy=ray.origin.y, oz=ray.origin.z;
    const float dx=ray.direction.x, dy=ray.direction.y, dz=ray.direction.z;
    
    const float e1x = bx - ax, e1y = by - ay, e1z = bz - az;
    const float e2x = cx - ax, e2y = cy - ay, e2z = cz - az;
    const float rx  = ox - ax, ry  = oy - ay, rz  = oz - az;
 
    float determinant =
        -( e1x * (e2y * dz - e2z * dy)
        - e1y * (e2x * dz - e2z * dx)
        + e1z * (e2x * dy - e2y * dx) );

    if (std::fabs(determinant) < intersectionTestEpsilon) {
        return false;
    }
    const float invDet = 1.0f / determinant;

    float determinantBeta =
        -( rx * (e2y * dz - e2z * dy)
        - ry * (e2x * dz - e2z * dx)
        + rz * (e2x * dy - e2y * dx) );
    float beta = determinantBeta * invDet;
    if (beta < 0.0f || beta > 1.0f) return false;
 
    float determinantGamma =
        -( e1x * (ry * dz - rz * dy)
        - e1y * (rx * dz - rz * dx)
        + e1z * (rx * dy - ry * dx) );
    float gamma = determinantGamma * invDet;
    if (gamma < 0.0f || gamma > 1.0f || beta + gamma > 1.0f) return false;

    if (ray.shadowRay || ray.reflectionRay || ray.refractionRay) {
        determinantT =
            e1x * (e2y * rz - e2z * ry)
            - e1y * (e2x * rz - e2z * rx)
            + e1z * (e2x * ry - e2y * rx);
    }
    float t = determinantT * invDet;
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
        return true;
    }
 
    return false;
}

bool rayHitsMesh(
    Ray& ray, 
    const Mesh& mesh, 
    const vector<VectorFloatTriplet>& vertices, 
    const vector<float>& determinants, 
    float& t_min, 
    Intersection& intersection,
    float intersectionTestEpsilon,
    scene::MeshBVH* bvh,
    bool enableBackFaceCulling,
    int meshIndex
) {
    if (bvh != nullptr) {
        return bvh->traverse(ray, mesh, vertices, determinants, t_min, intersection, intersectionTestEpsilon, enableBackFaceCulling, meshIndex);
    }
    bool hit = false;
    for(int i = 0; i < (int)mesh.faces.size(); i++) {
        if (rayHitsTriangle(ray, mesh.faces[i], vertices, t_min, intersection, intersectionTestEpsilon, determinants[i], mesh.material, enableBackFaceCulling, meshIndex, i)) {
            hit = true;
            intersection.kind = Intersection::Kind::Mesh;
        }
    }
    return hit;
}

Intersection intersect(const Scene& scene, Ray& ray) {
    float t_min = numeric_limits<float>::max();
    Intersection intersection;
    bool hit = false;
    
    // Test all geometry
    for(int i = 0; i < (int)scene.planes.size(); i++) {
        hit = rayHitsPlane(ray, scene.planes[i], scene.vertices, t_min, intersection, i) || hit;
    }
    for(int i = 0; i < (int)scene.spheres.size(); i++) {
        hit = rayHitsSphere(ray, scene.spheres[i], scene.vertices, t_min, intersection, i) || hit;
    }
    for(int i = 0; i < (int)scene.triangles.size(); i++) {
        if (rayHitsTriangle(ray, scene.triangles[i].indices, scene.vertices, t_min, intersection, scene.intersectionTestEpsilon, scene.cameraTriangleDeterminant[scene.currentCameraIndex][i], scene.triangles[i].material, scene.enableBackFaceCulling, -1, i)) {
            hit = true;
            intersection.kind = Intersection::Kind::Triangle;
        }
    }
    for(int i = 0; i < (int)scene.meshes.size(); i++) {
        MeshBVH* bvh = (i < scene.meshBVHs.size()) ? scene.meshBVHs[i] : nullptr;
        hit = rayHitsMesh(ray, scene.meshes[i], scene.vertices, scene.cameraMeshDeterminant[scene.currentCameraIndex][i], t_min, intersection, scene.intersectionTestEpsilon, bvh, scene.enableBackFaceCulling, i) || hit;
    }
    
    if (!hit) return intersection;
    const float u = intersection.beta;
    const float v = intersection.gamma;
    const float w = 1.0f - u - v;

    intersection.geometricNormal = normalize(intersection.geometricNormal);
    intersection.shadingNormal = intersection.geometricNormal;
    
    if (intersection.kind == Intersection::Kind::Mesh) {
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
                    intersection.shadingNormal = normalize(w * n0 + u * n1 + v * n2);
                    if (dotProduct(intersection.shadingNormal, intersection.geometricNormal) < 0.0f) {
                        intersection.shadingNormal = -intersection.shadingNormal;
                    }
                }
            }
        }
    }
    return intersection;
}

VectorFloatTriplet computeShading(const Scene& scene, Ray& ray, const Intersection& intersection) {
    Material* material = intersection.material;
    VectorFloatTriplet color = VectorFloatTriplet{0, 0, 0};
    
    // Ambient component
    VectorFloatTriplet ambient = material->ambientReflectance * scene.ambientLight.intensity;
    color += ambient;

    if (material->isMirror) {
        VectorFloatTriplet mirrorColor;
        Ray reflectedRay = reflect(ray, intersection.geometricNormal, intersection.point, scene.shadowRayEpsilon);
        Intersection reflectionIntersection = intersect(scene, reflectedRay);
        color +=  material->mirrorReflectance * (mirrorColor = computePixelColor(scene, reflectedRay, reflectionIntersection));
    }
    
    // Iterate through all point lights
    for (const PointLight& light : scene.pointLights) {
        if (!isInShadow(scene, ray, light, intersection)) {
            VectorFloatTriplet lightVec = light.position - intersection.point;
            float distance = sqrt(dotProduct(lightVec, lightVec));
            VectorFloatTriplet lightDir = normalize(lightVec);
            
            VectorFloatTriplet normal = intersection.shadingNormal;
            
            // Distance attenuation (1/r^2)
            float attenuation = 1.0f / (distance * distance);
            
            // Diffuse component
            float diffuseFactor = max(0.0f, dotProduct(normal, lightDir));
            VectorFloatTriplet diffuse = material->diffuseReflectance * light.intensity * diffuseFactor * attenuation;
            
            // Specular component (Blinn-Phong)
            VectorFloatTriplet viewDir = normalize(-ray.direction);
            VectorFloatTriplet halfVector = normalize(lightDir + viewDir);
            float specularFactor = pow(max(0.0f, dotProduct(normal, halfVector)), material->phongExponent);
            VectorFloatTriplet specular = material->specularReflectance * light.intensity * specularFactor * attenuation;
            
            color += diffuse + specular;
        }
    }
    
    return color;
}

VectorFloatTriplet computePixelColor(const Scene& scene, Ray& ray, const Intersection& intersection) {
    if (ray.depth > scene.maxRecursionDepth) {
        return VectorFloatTriplet({0, 0, 0});
    }
    if(intersection.hit) {
        return computeShading(scene, ray, intersection);
    } else if (ray.depth == 0) {
        return scene.backgroundColor;
    }
    return VectorFloatTriplet({0, 0, 0});
}

bool isInShadow(const Scene& scene, Ray& ray, const PointLight& light, const Intersection& intersection) {
    Ray shadowRay = Ray{
        intersection.point + scene.shadowRayEpsilon * intersection.geometricNormal,
        normalize(light.position - intersection.point),
        0,
        true, // shadow ray
        false, // reflection ray
        false // refraction ray
    };
    Intersection shadowIntersection = intersect(scene, shadowRay);
    float distToLight = sqrt(dotProduct(light.position - intersection.point, light.position - intersection.point));
    return shadowIntersection.hit && shadowIntersection.distance < distToLight;
}

Ray reflect(Ray& ray, const VectorFloatTriplet normal, VectorFloatTriplet point, const float shadowRayEpsilon) {
    Ray reflectedRay = Ray{
        point + shadowRayEpsilon * normalize(normal),
        // wr = -wo + 2ncosΘ = -wo + 2n(n.wo)
        ray.direction - 2 * dotProduct(ray.direction, normal) * normal,
        ray.depth + 1,
        false,
        true,
        false
    };
    return reflectedRay;
}