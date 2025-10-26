#include <math.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "scene.h"
#include "utils.h"
#include "overloads.h"

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
    Ray ray = Ray{origin, normalize(ray_direction), 0};
    return ray;
}

bool rayHitsPlane(const Ray& ray, const Plane& plane, const vector<VectorFloatTriplet>& vertices, float& t_min, Intersection& intersection) {
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
        intersection = Intersection{true, t, o + d * t, n, plane.material};
        return true;
       }
       return false;
}

bool rayHitsSphere(const Ray& ray, const Sphere& sphere, const vector<VectorFloatTriplet>& vertices, float& t_min, Intersection& intersection) {
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
        intersection = Intersection{true, t, intersection_point, intersection_normal, sphere.material};
        return true;
    }
    return false;
}


bool rayHitsTriangle(const Ray& ray, const VectorIntTriplet& face, const VectorFloatTriplet& triangleNormal, const vector<VectorFloatTriplet>& vertices, float& t_min, Intersection& intersection, float intersectionTestEpsilon, float determinantT, Material* material) {
    VectorFloatTriplet a = vertices[face.x];
    VectorFloatTriplet b = vertices[face.y];
    VectorFloatTriplet c = vertices[face.z];
    float area = dotProduct(crossProduct(b - a, c - a), triangleNormal);
    if (fabs(area) < 1e-8) return false; // degenerate triangle
    float ax=a.x, ay=a.y, az=a.z, bx=b.x, by=b.y, bz=b.z, cx=c.x, cy=c.y, cz=c.z;
    float ox=ray.origin.x, oy=ray.origin.y, oz=ray.origin.z, dx=ray.direction.x, dy=ray.direction.y, dz=ray.direction.z;
    /*
    For the sake of simplicity, I am naming some parts of the formula i.e. e1x = bx - ax, e1y = by - ay, e1z = bz - az, etc.

    @TODO: I need to discard the precomputed determinant whenever the ray is reflected or refracted. ORIGIN CHANGES.
    */
    const float e1x = bx - ax, e1y = by - ay, e1z = bz - az;
    const float e2x = cx - ax, e2y = cy - ay, e2z = cz - az;
    const float rx  = ox - ax, ry  = oy - ay, rz  = oz - az;
 
    float determinant =
        -( e1x * (e2y * dz - e2z * dy)
        - e1y * (e2x * dz - e2z * dx)
        + e1z * (e2x * dy - e2y * dx) );

    // Early quit
    if (std::fabs(determinant) < intersectionTestEpsilon) {
        return false;
    }
    const float invDet = 1.0f / determinant;

    float determinantBeta =
        -( rx * (e2y * dz - e2z * dy)
        - ry * (e2x * dz - e2z * dx)
        + rz * (e2x * dy - e2y * dx) );
    float beta = determinantBeta * invDet;

    // Early quit
    if (beta < intersectionTestEpsilon || beta > 1.0f) return false;
 
    float determinantGamma =
        -( e1x * (ry * dz - rz * dy)
        - e1y * (rx * dz - rz * dx)
        + e1z * (rx * dy - ry * dx) );
    float gamma = determinantGamma * invDet;

    // Early quits
    if (gamma < intersectionTestEpsilon || gamma > 1.0f) return false;
    if (beta + gamma > 1.0f) return false; 

    /*
    @TODO: I need to discard the precomputed determinant whenever the ray is reflected or refracted. ORIGIN CHANGES.
    determinantT =
          e1x * (e2y * rz - e2z * ry)
        - e1y * (e2x * rz - e2z * rx)
        + e1z * (e2x * ry - e2y * rx);
    */
    float t = determinantT * invDet;

     // Early quit
    if (t < intersectionTestEpsilon) return false;

    if(t < t_min) {
        t_min = t;
        intersection = Intersection{true, t, a + beta * (b - a) + gamma * (c - a), triangleNormal, material};
        return true;
    }
 
     return false;
}

bool rayHitsMesh(
    const Ray& ray, 
    const Mesh& mesh, 
    const vector<VectorFloatTriplet>& normals, 
    const vector<VectorFloatTriplet>& vertices, 
    const vector<float>& determinants, 
    float& t_min, 
    Intersection& intersection,
    float intersectionTestEpsilon
) {
    for(int i = 0; i < mesh.faces.size(); i++) {
        if(rayHitsTriangle(ray, mesh.faces[i], normals[i], vertices, t_min, intersection, intersectionTestEpsilon, determinants[i], mesh.material)) {
            return true;
        }
    }
    return false;
}

Intersection intersect(const Scene& scene, const Ray& ray) {
    float t_min = numeric_limits<float>::max();
    Intersection intersection = Intersection{false, 0, VectorFloatTriplet{0, 0, 0}, VectorFloatTriplet{0, 0, 0}, nullptr};
    bool hit = false;
    for(int i = 0; i < scene.planes.size(); i++) {
        hit = rayHitsPlane(ray, scene.planes[i], scene.vertices, t_min, intersection) || hit;
    }
    for(int i = 0; i < scene.spheres.size(); i++) {
        hit = rayHitsSphere(ray, scene.spheres[i], scene.vertices, t_min, intersection) || hit;
    }
    for(int i = 0; i < scene.triangles.size(); i++) {
        hit = rayHitsTriangle(ray, scene.triangles[i].indices, scene.triangleNormals[i], scene.vertices, t_min, intersection, scene.intersectionTestEpsilon, scene.cameraTriangleDeterminant[scene.currentCameraIndex][i], scene.triangles[i].material) || hit;
    }
    for(int i = 0; i < scene.meshes.size(); i++) {
        hit = rayHitsMesh(ray, scene.meshes[i], scene.meshNormals[i], scene.vertices, scene.cameraMeshDeterminant[scene.currentCameraIndex][i], t_min, intersection, scene.intersectionTestEpsilon) || hit;
    }
    return intersection;
}

VectorFloatTriplet computeShading(const Scene& scene, const Ray& ray, const Intersection& intersection) {
    Material* material = intersection.material;
    VectorFloatTriplet color = VectorFloatTriplet{0, 0, 0};
    
    // Ambient component
    VectorFloatTriplet ambient = material->ambientReflectance * scene.ambientLight.intensity;
    color += ambient;
    
    // Iterate through all point lights
    for (const PointLight& light : scene.pointLights) {
        VectorFloatTriplet lightVec = light.position - intersection.point;
        float distance = sqrt(dotProduct(lightVec, lightVec));
        VectorFloatTriplet lightDir = normalize(lightVec);
        VectorFloatTriplet normal = normalize(intersection.normal);
        
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
    
    return color;
}

VectorFloatTriplet computePixelColor(const Scene& scene, const Ray& ray, const Intersection& intersection) {
    if(intersection.hit) {
        return computeShading(scene, ray, intersection);
    }
    return scene.backgroundColor;
}
