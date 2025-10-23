
#include <iostream>
#include "scene.h"
#include "utils.h"
#define OUTPUT_PATH "../my_outputs/"
using namespace std;
using namespace scene;


int main(int argc, char* argv[])
{
    Scene scene;
    scene.loadSceneFromFile(argv[1]);

    scene.getSummary();


    vector<VectorFloatTriplet> meshNormals;
    vector<VectorFloatTriplet> triangleNormals;

    precomputeMeshNormals(scene.meshes, meshNormals, scene.vertices);
    precomputeTriangleNormals(scene.triangles, triangleNormals, scene.vertices);

    for (auto camera : scene.cameras) {
        int width = camera.imageResolution.x;
        int height = camera.imageResolution.y;
        vector<unsigned char> image(width * height * 3);
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                /*
                    TODO: I have to do the following to compute pixelColor:
                    1. Cast a ray from the camera to the pixel
                    2. Check if the ray intersects with any object, I'll prolly need to write separate function for different objects
                    3. If it intersects, compute the color of the pixel
                    Somewhere around here I gotta apply shading
                    4. If it doesn't intersect, return the background color
                */
                Ray ray = castRay(camera, x, y, width, height);
                Intersection intersection = intersect(scene, ray, meshNormals, triangleNormals);
                VectorFloatTriplet pixelColor = computePixelColor(scene, ray, intersection);
                image[y * width + x * 3] = pixelColor.x;
                image[y * width + x * 3 + 1] = pixelColor.y;
                image[y * width + x * 3 + 2] = pixelColor.z;
            }
        }
        writePPM((OUTPUT_PATH + camera.imageName).c_str(), image, width, height);
    }
    return 0;
}
 