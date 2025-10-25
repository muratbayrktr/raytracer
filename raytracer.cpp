
#include <iostream>
#include "scene.h"
#include "utils.h"
#include <cmath>
#define OUTPUT_PATH "../my_outputs/"
using namespace std;
using namespace scene;



int main(int argc, char* argv[])
{
    Scene scene;
    scene.loadSceneFromFile(argv[1]);

    scene.getSummary();


    vector<vector<VectorFloatTriplet>> meshNormals;
    vector<VectorFloatTriplet> triangleNormals;

    precomputeMeshNormals(scene.meshes, meshNormals, scene.vertices);
    precomputeTriangleNormals(scene.triangles, triangleNormals, scene.vertices);

    for (auto camera : scene.cameras) {
        int width = camera.imageResolution.x;
        int height = camera.imageResolution.y;
        unsigned char* image = new unsigned char[width * height * 3];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                Ray ray = castRay(camera, x, y, width, height);
                Intersection intersection = intersect(scene, ray, meshNormals, triangleNormals);
                VectorFloatTriplet pixelColor = computePixelColor(scene, ray, intersection, meshNormals, triangleNormals);
                
                image[(y * width + x) * 3] = (unsigned char) round(pixelColor.x);
                image[(y * width + x) * 3 + 1] = (unsigned char) round(pixelColor.y);
                image[(y * width + x) * 3 + 2] = (unsigned char) round(pixelColor.z);
            }
        }
        string outputName = camera.imageName;
        writePPM((OUTPUT_PATH + outputName).c_str(), image, width, height);
        delete[] image;
    }
    return 0;
}
 