
#include <iostream>
#include "scene.h"
#include "utils.h"
#include <cmath>
#include "precompute.h"
#include <thread>
#include <vector>
#include <chrono>
#define OUTPUT_PATH "../my_outputs/"
using namespace std;
using namespace scene;


VectorFloatTriplet __compute(Scene& scene, Camera& camera, int x, int y, int width, int height) {
    Ray ray = castRay(camera, x, y, width, height);
    Intersection intersection = intersect(scene, ray);
    VectorFloatTriplet pixelColor = computePixelColor(scene, ray, intersection);
    clamp(pixelColor, 0, 255);
    return pixelColor;
}

void multiThreadedRayTracing(Scene& scene, Camera& camera, int width, int height, unsigned char* image) {
    auto start = std::chrono::high_resolution_clock::now();
    unsigned int nThreads = std::thread::hardware_concurrency();
    if (nThreads == 0) nThreads = 4;
    std::cout << "Number of threads: " << nThreads << std::endl;

    auto threadFunction = [&](int startY, int endY) {
        for (int y = startY; y < endY; ++y) {
            for (int x = 0; x < width; ++x) {
                VectorFloatTriplet pixelColor = __compute(scene, camera, x, y, width, height);
                image[(y * width + x) * 3] = (unsigned char) round(pixelColor.x);
                image[(y * width + x) * 3 + 1] = (unsigned char) round(pixelColor.y);
                image[(y * width + x) * 3 + 2] = (unsigned char) round(pixelColor.z);
            }
        }
    };

    std::vector<std::thread> threads;
    int rowsPerThread = height / nThreads;
    int extra = height % nThreads;
    int currentY = 0;

    for (unsigned int t = 0; t < nThreads; ++t) {
        int startY = currentY;
        int endY = startY + rowsPerThread + (t < extra ? 1 : 0);
        threads.emplace_back(threadFunction, startY, endY);
        currentY = endY;
    }
    for (auto& t : threads) {
        t.join();
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Multi-threaded ray tracing time: " << duration.count() << " milliseconds" << std::endl;
}

void singleThreadedRayTracing(Scene& scene, Camera& camera, int width, int height, unsigned char* image) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            VectorFloatTriplet pixelColor = __compute(scene, camera, x, y, width, height);
            image[(y * width + x) * 3] = (unsigned char) round(pixelColor.x);
            image[(y * width + x) * 3 + 1] = (unsigned char) round(pixelColor.y);
            image[(y * width + x) * 3 + 2] = (unsigned char) round(pixelColor.z);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Single-threaded ray tracing time: " << duration.count() << " milliseconds" << std::endl;
}


int main(int argc, char* argv[])
{
    Scene scene;
    scene.loadSceneFromFile(argv[1]);

    scene.getSummary();


    vector<vector<VectorFloatTriplet>> meshNormals;
    vector<VectorFloatTriplet> triangleNormals;
    vector<vector<float>> cameraTriangleDeterminant;
    vector<vector<vector<float>>> cameraMeshDeterminant;

    precomputeMeshNormals(scene.meshes, meshNormals, scene.vertices);
    precomputeTriangleNormals(scene.triangles, triangleNormals, scene.vertices);
    precomputeCameraTriangleDeterminant(scene, cameraTriangleDeterminant);
    precomputeCameraMeshDeterminant(scene, cameraMeshDeterminant);

    scene.cameraTriangleDeterminant = cameraTriangleDeterminant;
    scene.cameraMeshDeterminant = cameraMeshDeterminant;
    scene.meshNormals = meshNormals;
    scene.triangleNormals = triangleNormals;


    for (int i = 0; i < scene.cameras.size(); i++) {
        Camera camera = scene.cameras[i];
        scene.currentCameraIndex = i;
        int width = camera.imageResolution.x;
        int height = camera.imageResolution.y;
        unsigned char* image = new unsigned char[width * height * 3];
        if (argc > 2 && strcmp(argv[2], "multi") == 0) {
            multiThreadedRayTracing(scene, camera, width, height, image);
        } else {
            singleThreadedRayTracing(scene, camera, width, height, image);
        }
        string outputName = camera.imageName;
        scene.writePPM((OUTPUT_PATH + outputName).c_str(), image, width, height);
        delete[] image;
    }
    return 0;
}
 