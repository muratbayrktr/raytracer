
#include <iostream>
#include "scene.h"
#include "utils.h"
#include <cmath>
#include "precompute.h"
#include <pthread.h>
#include <vector>
#include <chrono>
#include <unistd.h>
#define OUTPUT_PATH "./"
using namespace std;
using namespace scene;


VectorFloatTriplet __compute(Scene& scene, Camera& camera, int x, int y, int width, int height) {
    Ray ray = castRay(camera, x, y, width, height);
    Intersection intersection = intersect(scene, ray);
    VectorFloatTriplet pixelColor = computePixelColor(scene, ray, intersection);
    clamp(pixelColor, 0, 255);
    return pixelColor;
}

// Thread argument structure for pthread
struct ThreadArgs {
    Scene* scene;
    Camera* camera;
    int startY;
    int endY;
    int width;
    int height;
    unsigned char* image;
};

// Thread function for pthread
void* threadFunction(void* arg) {
    ThreadArgs* args = (ThreadArgs*)arg;
    Scene* scene = args->scene;
    Camera* camera = args->camera;
    
    for (int y = args->startY; y < args->endY; ++y) {
        for (int x = 0; x < args->width; ++x) {
            VectorFloatTriplet pixelColor = __compute(*scene, *camera, x, y, args->width, args->height);
            args->image[(y * args->width + x) * 3] = (unsigned char) round(pixelColor.x);
            args->image[(y * args->width + x) * 3 + 1] = (unsigned char) round(pixelColor.y);
            args->image[(y * args->width + x) * 3 + 2] = (unsigned char) round(pixelColor.z);
        }
    }
    return NULL;
}

void multiThreadedRayTracing(Scene& scene, Camera& camera, int width, int height, unsigned char* image) {
    auto start = std::chrono::high_resolution_clock::now();

    long nThreads = sysconf(_SC_NPROCESSORS_ONLN);
    if (nThreads <= 0) nThreads = 4;
    std::cout << "Number of threads: " << nThreads << std::endl;

    std::vector<pthread_t> threads(nThreads);
    std::vector<ThreadArgs> threadArgs(nThreads);
    
    int rowsPerThread = height / nThreads;
    int extra = height % nThreads;
    int currentY = 0;

    for (long t = 0; t < nThreads; ++t) {
        int startY = currentY;
        int endY = startY + rowsPerThread + (t < extra ? 1 : 0);
        
        threadArgs[t].scene = &scene;
        threadArgs[t].camera = &camera;
        threadArgs[t].startY = startY;
        threadArgs[t].endY = endY;
        threadArgs[t].width = width;
        threadArgs[t].height = height;
        threadArgs[t].image = image;
        
        pthread_create(&threads[t], NULL, threadFunction, &threadArgs[t]);
        currentY = endY;
    }

    for (long t = 0; t < nThreads; ++t) {
        pthread_join(threads[t], NULL);
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

scene::Args parseArgs(int argc, char* argv[]) {
    scene::Args args;
    args.sceneFile = argv[1];
    args.isMultiThreaded = true;
    args.useBVH = true;
    args.enableBackFaceCulling = true;
    
    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], "-m") == 0) {
            args.isMultiThreaded = false;
        }
        if (strcmp(argv[i], "-b") == 0) {
            args.useBVH = false;
        }
        if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--no-cull") == 0) {
            args.enableBackFaceCulling = false;
        }
    }
    return args;
}


int main(int argc, char* argv[])
{

    scene::Args args = parseArgs(argc, argv);

    Scene scene;
    scene.loadSceneFromFile(args.sceneFile);
    scene.enableBackFaceCulling = args.enableBackFaceCulling;

    scene.getSummary();


    vector<vector<VectorFloatTriplet>> meshVertexNormals;
    vector<VectorFloatTriplet> triangleNormals;
    vector<vector<float>> cameraTriangleDeterminant;
    vector<vector<vector<float>>> cameraMeshDeterminant;

    precomputeMeshNormals(scene.meshes, meshVertexNormals, scene.vertices);
    precomputeTriangleNormals(scene.triangles, triangleNormals, scene.vertices);
    precomputeCameraTriangleDeterminant(scene, cameraTriangleDeterminant);
    precomputeCameraMeshDeterminant(scene, cameraMeshDeterminant);

    scene.cameraTriangleDeterminant = cameraTriangleDeterminant;
    scene.cameraMeshDeterminant = cameraMeshDeterminant;
    scene.meshVertexNormals = meshVertexNormals;
    scene.triangleNormals = triangleNormals;

    if (args.useBVH) {
        scene.buildBVH();
    }

    for (int i = 0; i < scene.cameras.size(); i++) {
        Camera camera = scene.cameras[i];
        scene.currentCameraIndex = i;
        int width = camera.imageResolution.x;
        int height = camera.imageResolution.y;
        unsigned char* image = new unsigned char[width * height * 3];
        if (args.isMultiThreaded) {
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
 