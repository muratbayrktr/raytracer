#include <math.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "scene.h"
#include "utils.h"

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

void precomputeCameraTriangleDeterminant(const Scene& scene, vector<vector<float>>& cameraTriangleDeterminant) {
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

void precomputeCameraMeshDeterminant(const Scene& scene, vector<vector<vector<float>>>& cameraMeshDeterminant) {
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

void writePPM(const string& filename, unsigned char* image, int width, int height) {
    // Remove following lines to write png file, I did this for debugging
    size_t dotPos = filename.find_last_of('.');
    string ppmFilename;
    if (dotPos != string::npos) {
        ppmFilename = filename.substr(0, dotPos) + ".ppm";
    }
    FILE *outfile;
    if ((outfile = fopen(ppmFilename.c_str(), "w")) == NULL) {
        throw runtime_error("Error: The ppm file cannot be opened for writing: " + filename);
    }
    fprintf(outfile, "P3\n%d %d\n255\n", width, height);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < 3; c++) {
                fprintf(outfile, "%u ", image[(y * width + x) * 3 + c]);
            }
        }
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
    fclose(outfile);
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
    /*
        The idea is to go from the area of the triangle to the barycentric coordinates of the intersection point.
        This will allow me to early quit with some checks as well.

        We can sit down and write cramer's rule function or determinant function to solve the equation.
        HOWEVER, I want speed and I won't be calculating the bigger determinants or I won't need cramer for other things.
        Instead, I will just plug in the formula for alpha, beta and gamma in the derived form.
    */
    VectorFloatTriplet a = vertices[face.x];
    VectorFloatTriplet b = vertices[face.y];
    VectorFloatTriplet c = vertices[face.z];
    float area = dotProduct(crossProduct(b - a, c - a), triangleNormal);
    if (fabs(area) < 1e-8) return false; // degenerate triangle
    float ax=a.x, ay=a.y, az=a.z, bx=b.x, by=b.y, bz=b.z, cx=c.x, cy=c.y, cz=c.z;
    float ox=ray.origin.x, oy=ray.origin.y, oz=ray.origin.z, dx=ray.direction.x, dy=ray.direction.y, dz=ray.direction.z;
    /*

    Taken from lecture slides:
    beta = | ax-ox ax-cx dx |
           | ay-oy ay-cy dy |
           | az-oz az-cz dz |
           -------------------
                  |area|

    gamma = | ax-bx ax-ox dx |
           | ay-by ay-oy dy |
           | az-bz az-oz dz |
           -------------------
                  |area|

    alpha = 1 - beta - gamma

    determinantBeta = (ax - ox) * [(ay-cy)*dz - (az-cz)*dy] 
                    + (ay - oy) * [(az-cz)*dx - (ax-cx)*dz]
                    + (az - oz) * [(ax-cx)*dy - (ay-cy)*dx]
    

    determinantGamma = (ax - bx) * [(ay-oy)*dz - (az-oz)*dy] 
                    + (ay - by) * [(az-oz)*dx - (ax-ox)*dz]
                    + (az - bz) * [(ax-ox)*dy - (ay-oy)*dx]
                    
    beta = determinantBeta / |area|
    gamma = determinantGamma / |area|
    alpha = 1 - beta - gamma

    For the t calculation we have (replacing first column with o-a):

    determinantT = | ax-bx ax-cx ax-ox |
                   | ay-by ay-cy ay-oy |
                   | az-bz az-cz az-oz |

    determinantT = (ax-bx) * [(ay-cy)*(az-oz) - (az-cz)*(ay-oy)]
                 + (ay-by) * [(az-cz)*(ax-ox) - (ax-cx)*(az-oz)]
                 + (az-bz) * [(ax-cx)*(ay-oy) - (ay-cy)*(ax-ox)]

    t = determinantT / area

    For the sake of simplicity, I am naming some parts of the formula i.e. e1x = bx - ax, e1y = by - ay, e1z = bz - az, etc.
    This is done to avoid having to write the same formula multiple times.
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

/* Operator Overloads */
/* VectorFloatTriplet */
VectorFloatTriplet crossProduct(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; 
}

float dotProduct(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return a.x * b.x + a.y * b.y + a.z * b.z; 
}

VectorFloatTriplet normalize(const VectorFloatTriplet& a) { 
    float length = sqrt(a.x * a.x + a.y * a.y + a.z * a.z); 
    return VectorFloatTriplet{a.x / length, a.y / length, a.z / length}; 
}

VectorFloatTriplet operator-(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.x - b.x, a.y - b.y, a.z - b.z}; 
}

VectorFloatTriplet operator+(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.x + b.x, a.y + b.y, a.z + b.z}; 
}

VectorFloatTriplet operator*(const VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a.x * b.x, a.y * b.y, a.z * b.z}; 
}

VectorFloatTriplet& operator+=(VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    a.x += b.x; a.y += b.y; a.z += b.z;
    return a;
}

VectorFloatTriplet& operator-=(VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
    return a;
}

VectorFloatTriplet& operator*=(VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
    return a;
}

VectorFloatTriplet& operator/=(VectorFloatTriplet& a, const VectorFloatTriplet& b) { 
    a.x /= b.x; a.y /= b.y; a.z /= b.z;
    return a;
}

VectorFloatTriplet operator-(const VectorFloatTriplet& a) { 
    return VectorFloatTriplet{-a.x, -a.y, -a.z}; 
}

VectorFloatTriplet operator*(const VectorFloatTriplet& a, const float& b) { 
    return VectorFloatTriplet{a.x * b, a.y * b, a.z * b}; 
}

VectorFloatTriplet operator*(const float& a, const VectorFloatTriplet& b) { 
    return VectorFloatTriplet{a * b.x, a * b.y, a * b.z}; 
}
/* VectorIntTriplet */
VectorIntTriplet crossProduct(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; 
}

int dotProduct(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return a.x * b.x + a.y * b.y + a.z * b.z; 
}

VectorIntTriplet normalize(const VectorIntTriplet& a) { 
    int length = sqrt(a.x * a.x + a.y * a.y + a.z * a.z); 
    return VectorIntTriplet{a.x / length, a.y / length, a.z / length}; 
}

VectorIntTriplet operator-(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.x - b.x, a.y - b.y, a.z - b.z}; 
}

VectorIntTriplet operator+(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.x + b.x, a.y + b.y, a.z + b.z}; 
}

VectorIntTriplet operator*(const VectorIntTriplet& a, const VectorIntTriplet& b) { 
    return VectorIntTriplet{a.x * b.x, a.y * b.y, a.z * b.z}; 
}

VectorIntTriplet& operator+=(VectorIntTriplet& a, const VectorIntTriplet& b) { 
    a.x += b.x; a.y += b.y; a.z += b.z;
    return a;
}

VectorIntTriplet& operator-=(VectorIntTriplet& a, const VectorIntTriplet& b) { 
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
    return a;
}

VectorIntTriplet& operator*=(VectorIntTriplet& a, const VectorIntTriplet& b) { 
    a.x *= b.x; a.y *= b.y; a.z *= b.z;
    return a;
}

VectorIntTriplet& operator/=(VectorIntTriplet& a, const VectorIntTriplet& b) { 
    a.x /= b.x; a.y /= b.y; a.z /= b.z;
    return a;
}
