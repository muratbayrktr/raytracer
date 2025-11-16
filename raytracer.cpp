
#include <iostream>
#include "scene.h"
#include "utils.h"
#include <cmath>
#include "precompute.h"
#include <pthread.h>
#include <vector>
#include <chrono>
#include <unistd.h>
#include <fstream>
#include <atomic>
#define OUTPUT_PATH "../my_outputs_hw2/"
using namespace std;
using namespace scene;

atomic<int> g_pixelsProcessed(0);
int g_totalPixels = 0;

void printProgress(int processed, int total, std::chrono::time_point<std::chrono::high_resolution_clock> startTime) {
    double percentage = (100.0 * processed) / total;
    int barWidth = 50;
    int pos = barWidth * processed / total;
    
    auto now = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count();
    double pixelsPerMs = processed / (double)(elapsed + 1);
    int remaining = (total - processed) / (pixelsPerMs + 0.001f);
    
    std::cout << "\r[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(percentage) << "% (" << processed << "/" << total << ") ";
    std::cout << "ETA: " << remaining / 1000 << "s";
    
    // Print inline performance stats
    printPerfStatsInline();
    
    std::cout << "   " << std::flush;
}

VectorFloatTriplet __compute(Scene& scene, Camera& camera, int x, int y, int width, int height) {
    Ray ray = castRay(camera, x, y, width, height);
    Intersection intersection = intersect(scene, ray);
    VectorFloatTriplet pixelColor = computePixelColor(scene, ray, intersection);
    clamp(pixelColor, 0, 255);
    return pixelColor;
}

void printTimingStats() {
    std::cout << "\n=== Timing Statistics ===" << std::endl;
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
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
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
            
            int processed = ++g_pixelsProcessed;
            if (processed % 10000 == 0 || processed == g_totalPixels) {
                printPerfStats();
                printPerfStatsInline();
                printProgress(processed, g_totalPixels, args->startTime);
            }
        }
    }
    return NULL;
}

double multiThreadedRayTracing(Scene& scene, Camera& camera, int width, int height, unsigned char* image) {
    auto start = std::chrono::high_resolution_clock::now();

    g_totalPixels = width * height;
    g_pixelsProcessed = 0;

    long nThreads = sysconf(_SC_NPROCESSORS_ONLN);
    if (nThreads <= 0 || nThreads > 4) nThreads = 4;
    std::cout << "Number of threads: " << nThreads << std::endl;
    std::cout << "Rendering " << width << "x" << height << " (" << g_totalPixels << " pixels)..." << std::endl;

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
        threadArgs[t].startTime = start;
        
        pthread_create(&threads[t], NULL, threadFunction, &threadArgs[t]);
        currentY = endY;
    }

    for (long t = 0; t < nThreads; ++t) {
        pthread_join(threads[t], NULL);
    }
    
    std::cout << std::endl;
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Multi-threaded ray tracing time: " << duration.count() << " milliseconds" << std::endl;
    
    printPerfStats();
    
    return duration.count();
}

double singleThreadedRayTracing(Scene& scene, Camera& camera, int width, int height, unsigned char* image) {
    auto start = std::chrono::high_resolution_clock::now();
    
    g_totalPixels = width * height;
    g_pixelsProcessed = 0;
    
    std::cout << "Rendering " << width << "x" << height << " (" << g_totalPixels << " pixels)..." << std::endl;
    
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            VectorFloatTriplet pixelColor = __compute(scene, camera, x, y, width, height);
            image[(y * width + x) * 3] = (unsigned char) round(pixelColor.x);
            image[(y * width + x) * 3 + 1] = (unsigned char) round(pixelColor.y);
            image[(y * width + x) * 3 + 2] = (unsigned char) round(pixelColor.z);
            
            int processed = ++g_pixelsProcessed;
            if (processed % 10000 == 0 || processed == g_totalPixels) {
                printProgress(processed, g_totalPixels, start);
            }
        }
    }
    
    std::cout << std::endl;
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Single-threaded ray tracing time: " << duration.count() << " milliseconds" << std::endl;
    return duration.count();
}

scene::Args parseArgs(int argc, char* argv[]) {
    scene::Args args;
    args.sceneFile = argv[1];
    args.isMultiThreaded = true;
    args.useBVH = true;
    args.enableBackFaceCulling = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-m") == 0) {
            // Accept either -m 0/1 or just -m (legacy set to false)
            if (i + 1 < argc && (strcmp(argv[i + 1], "0") == 0 || strcmp(argv[i + 1], "1") == 0)) {
                args.isMultiThreaded = (atoi(argv[i + 1]) != 0);
                i++;
            } else {
                args.isMultiThreaded = false;
            }
        }
        else if (strcmp(argv[i], "-b") == 0) {
            // Accept either -b 0/1 or just -b (legacy set to true)
            if (i + 1 < argc && (strcmp(argv[i + 1], "0") == 0 || strcmp(argv[i + 1], "1") == 0)) {
                args.useBVH = (atoi(argv[i + 1]) != 0);
                i++;
            } else {
                args.useBVH = true;
            }
        }
        else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--no-cull") == 0) {
            // Accept either -c 0/1 or just -c (legacy set to false)
            if (i + 1 < argc && (strcmp(argv[i + 1], "0") == 0 || strcmp(argv[i + 1], "1") == 0)) {
                args.enableBackFaceCulling = (atoi(argv[i + 1]) != 0);
                i++;
            } else {
                args.enableBackFaceCulling = false;
            }
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
    vector<vector<double>> cameraTriangleDeterminant;
    vector<vector<vector<double>>> cameraMeshDeterminant;

    auto preprocessingStart = std::chrono::high_resolution_clock::now();
    
    // IMPORTANT: Build BVH first, as transformation precomputation needs it for world-space bounds
    if (args.useBVH) {
        scene.buildBVH();
    }
    
    scene.precomputeTransformations();

    precomputeMeshNormals(scene.meshes, meshVertexNormals, scene.vertices);
    precomputeTriangleNormals(scene.triangles, triangleNormals, scene.vertices);
    precomputeCameraTriangleDeterminant(scene, cameraTriangleDeterminant);
    precomputeCameraMeshDeterminant(scene, cameraMeshDeterminant);

    scene.cameraTriangleDeterminant = cameraTriangleDeterminant;
    scene.cameraMeshDeterminant = cameraMeshDeterminant;
    scene.meshVertexNormals = meshVertexNormals;
    scene.triangleNormals = triangleNormals;
    auto preprocessingEnd = std::chrono::high_resolution_clock::now();
    auto preprocessingTime = std::chrono::duration_cast<std::chrono::milliseconds>(preprocessingEnd - preprocessingStart);
    double preprocessingTimeMs = preprocessingTime.count();
    double totalTimeMs = preprocessingTimeMs;
    double renderTimeMs = 0;
    for (int i = 0; i < scene.cameras.size(); i++) {
        Camera camera = scene.cameras[i];
        scene.currentCameraIndex = i;
        int width = camera.imageResolution.x;
        int height = camera.imageResolution.y;
        unsigned char* image = new unsigned char[width * height * 3];
        if (args.isMultiThreaded) {
            renderTimeMs = multiThreadedRayTracing(scene, camera, width, height, image);
            totalTimeMs += renderTimeMs;
        } else {
            renderTimeMs = singleThreadedRayTracing(scene, camera, width, height, image);
            totalTimeMs += renderTimeMs;
        }
        string outputName = camera.imageName;
        scene.writePPM((OUTPUT_PATH + outputName).c_str(), image, width, height);
        delete[] image;

        json results = {
            {"sceneName", outputName},
            {"preprocessingTimeMs", preprocessingTimeMs},
            {"renderTimeMs", renderTimeMs},
            {"totalTimeMs", totalTimeMs},
            {"useBVH", args.useBVH},
            {"isMultiThreaded", args.isMultiThreaded},
            {"enableBackFaceCulling", args.enableBackFaceCulling}
        };
    
        // output to scene name file_results.json
        std::ofstream resultsFile((OUTPUT_PATH + outputName.substr(0, outputName.find_last_of('.')) + "_results.json").c_str());
        resultsFile << results.dump(4);
        resultsFile.close();
    }

    return 0;
}