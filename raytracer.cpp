
#include <iostream>
#include "scene.h"
#include "utils.h"
#include <cmath>
#include "precompute.h"
#include "overloads.h"
#include <pthread.h>
#include <vector>
#include <chrono>
#include <unistd.h>
#include <fstream>
#include <atomic>
// #include <SDL3/SDL.h> uncomment this if sdl3 is installed
#define OUTPUT_PATH "../my_outputs_hw4/"
#define JSON_OUTPUT_PATH "../my_outputs_hw4/benchmark/"
using namespace std;
using namespace scene;

atomic<int> g_pixelsProcessed(0);
int g_totalPixels = 0;

#define USE_GUI 0
#if USE_GUI
#include <SDL3/SDL.h>
#endif

#if USE_GUI
struct GuiState {
    bool initialized = false;
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;
    SDL_Texture* texture = nullptr;
    int width = 0;
    int height = 0;
};

static GuiState g_gui;
static bool g_useGui = false;
static bool g_guiInitTried = false;

static bool initGui(int width, int height) {
    // If we already tried and failed, don't spam attempts.
    if (g_guiInitTried && !g_gui.initialized) {
        return false;
    }
    if (g_gui.initialized && g_gui.width == width && g_gui.height == height) {
        return true;
    }

    if (!g_gui.initialized) {
        if (!SDL_Init(SDL_INIT_VIDEO)) {
            g_guiInitTried = true;
            std::cerr << "SDL_Init Error: " << SDL_GetError() << " (disabling GUI)" << std::endl;
            g_useGui = false;
            exit(1); // don't remove until we fix errors
            return false;
        }
    } else {
        if (g_gui.texture) {
            SDL_DestroyTexture(g_gui.texture);
            g_gui.texture = nullptr;
        }
        if (g_gui.renderer) {
            SDL_DestroyRenderer(g_gui.renderer);
            g_gui.renderer = nullptr;
        }
        if (g_gui.window) {
            SDL_DestroyWindow(g_gui.window);
            g_gui.window = nullptr;
        }
    }

    // SDL3: SDL_CreateWindow(title, w, h, flags)
    g_gui.window = SDL_CreateWindow("Raytracer", width, height, 0);
    if (!g_gui.window) {
        std::cerr << "SDL_CreateWindow Error: " << SDL_GetError() << std::endl;
        return false;
    }

    // SDL3: simplified renderer creation, no flags needed for our use
    g_gui.renderer = SDL_CreateRenderer(g_gui.window, nullptr);
    if (!g_gui.renderer) {
        std::cerr << "SDL_CreateRenderer Error: " << SDL_GetError() << std::endl;
        return false;
    }

    g_gui.texture = SDL_CreateTexture(g_gui.renderer,
                                      SDL_PIXELFORMAT_RGB24,
                                      SDL_TEXTUREACCESS_STREAMING,
                                      width,
                                      height);
    if (!g_gui.texture) {
        std::cerr << "SDL_CreateTexture Error: " << SDL_GetError() << std::endl;
        return false;
    }

    g_gui.width = width;
    g_gui.height = height;
    g_gui.initialized = true;
    return true;
}

static void updateGui(const unsigned char* image, int width, int height) {
    if (!g_useGui) return;
    if (!initGui(width, height)) return;

    SDL_Event e;
    while (SDL_PollEvent(&e)) {
        if (e.type == SDL_EVENT_QUIT) {
            // For now, just ignore and let rendering continue.
        }
    }

    SDL_UpdateTexture(g_gui.texture, nullptr, image, width * 3);
    SDL_RenderClear(g_gui.renderer);
    SDL_RenderTexture(g_gui.renderer, g_gui.texture, nullptr, nullptr);
    SDL_RenderPresent(g_gui.renderer);
}
#endif

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

VectorFloatTriplet __compute(Scene& scene, Camera& camera, double x, double y, int width, int height, double time, double random1, double random2) {
    Ray ray = castRay(camera, x, y, width, height, time, random1, random2);
    Intersection intersection = intersect(scene, ray);
    VectorFloatTriplet pixelColor = computePixelColor(scene, ray, intersection);
    clamp(pixelColor, 0, 255);
    return pixelColor;
}

void printTimingStats() {
    std::cout << "\n=== Timing Statistics ===" << std::endl;
}

// Thread argument structure for pthread (non-iterative sampling)
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

struct IterThreadArgs {
    Scene* scene;
    Camera* camera;
    int startY;
    int endY;
    int width;
    int height;
    VectorFloatTriplet* accum;
    int currentSample;
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
};

// Thread function for pthread
void* threadFunction(void* arg) {
    ThreadArgs* args = (ThreadArgs*)arg;
    Scene* scene = args->scene;
    Camera* camera = args->camera;
    
    for (int y = args->startY; y < args->endY; ++y) {
        for (int x = 0; x < args->width; ++x) {
            VectorFloatTriplet pixelColor = {0.0, 0.0, 0.0};
            int pixelIndex = y * args->width + x;
            int sampleIndex = pixelIndex * camera->numSamples;
            for (int k = 0; k < camera->numSamples; k++) {
                VectorFloatPenta sample = camera->samples[sampleIndex + k];
                double sx = x + sample.x;   // subpixel x
                double sy = y + sample.y;   // subpixel y
                // Use precomputed time and extra random dims per sample
                double sampleTime = sample.z;
                double r1 = sample.w;
                double r2 = sample.v;
                pixelColor = pixelColor + __compute(*scene, *camera, sx, sy, args->width, args->height, sampleTime, r1, r2);
            }
            pixelColor = pixelColor * (1.0 / camera->numSamples);
            args->image[pixelIndex * 3 + 0] = (unsigned char) round(pixelColor.x);
            args->image[pixelIndex * 3 + 1] = (unsigned char) round(pixelColor.y);
            args->image[pixelIndex * 3 + 2] = (unsigned char) round(pixelColor.z);
            
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

// Thread function for iterative sampling (one sample index over many pixels)
void* threadFunctionIterative(void* arg) {
    IterThreadArgs* args = (IterThreadArgs*)arg;
    Scene* scene = args->scene;
    Camera* camera = args->camera;
    int k = args->currentSample;
    
    for (int y = args->startY; y < args->endY; ++y) {
        for (int x = 0; x < args->width; ++x) {
            int pixelIndex = y * args->width + x;
            int sampleIndexBase = pixelIndex * camera->numSamples;
            VectorFloatPenta sample = camera->samples[sampleIndexBase + k];
            double sx = x + sample.x;
            double sy = y + sample.y;
            double sampleTime = sample.z;
            double r1 = sample.w;
            double r2 = sample.v;
            VectorFloatTriplet color =
                __compute(*scene, *camera, sx, sy, args->width, args->height, sampleTime, r1, r2);
            // Accumulate (no race: each pixel belongs to exactly one thread)
            args->accum[pixelIndex].x += color.x;
            args->accum[pixelIndex].y += color.y;
            args->accum[pixelIndex].z += color.z;
            
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

// Helper: build intermediate image from accumulation buffer and write canvas PNG
static void writeIterativeCanvas(Scene& scene,
                                 Camera& camera,
                                 VectorFloatTriplet* accum,
                                 int width,
                                 int height,
                                 int samplesSoFar,
                                 unsigned char* image) {
    double invSamples = 1.0 / std::max(1, samplesSoFar);
    for (int i = 0; i < width * height; ++i) {
        VectorFloatTriplet c = accum[i];
        c.x *= invSamples;
        c.y *= invSamples;
        c.z *= invSamples;
        clamp(c, 0, 255);
        image[i * 3 + 0] = (unsigned char) std::round(c.x);
        image[i * 3 + 1] = (unsigned char) std::round(c.y);
        image[i * 3 + 2] = (unsigned char) std::round(c.z);
    }
    
    std::string baseName = camera.imageName;
    size_t dotPos = baseName.find_last_of('.');
    if (dotPos != std::string::npos) {
        baseName = baseName.substr(0, dotPos);
    }
    std::string canvasName = OUTPUT_PATH + baseName + "_iterative.png";
    scene.writePPM(canvasName.c_str(), image, width, height);

    // Also update GUI window if enabled
#if USE_SDL3
    updateGui(image, width, height);
#endif
}

double multiThreadedRayTracing(Scene& scene, Camera& camera, int width, int height, unsigned char* image, bool iterativeSampling) {
    auto start = std::chrono::high_resolution_clock::now();

    int numSamples = camera.numSamples;
    if (iterativeSampling) {
        g_totalPixels = width * height * numSamples;
    } else {
        g_totalPixels = width * height;
    }
    g_pixelsProcessed = 0;

    long nThreads = sysconf(_SC_NPROCESSORS_ONLN);
    if (nThreads <= 0 || nThreads > 4) nThreads = 4;
    std::cout << "Number of threads: " << nThreads << std::endl;
    std::cout << "Rendering " << width << "x" << height << " (" << g_totalPixels << " pixels)..." << std::endl;

    if (!iterativeSampling) {
        // Original per-pixel sampling (all samples per pixel, then move on)
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
    } else {
        // Iterative, level-wise sampling over samples, with canvas update every 2 samples
        std::vector<pthread_t> threads(nThreads);
        std::vector<IterThreadArgs> threadArgs(nThreads);
        VectorFloatTriplet* accum = new VectorFloatTriplet[width * height];
        for (int i = 0; i < width * height; ++i) {
            accum[i].x = accum[i].y = accum[i].z = 0.0;
        }

        int rowsPerThread = height / nThreads;
        int extra = height % nThreads;

        for (int k = 0; k < numSamples; ++k) {
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
                threadArgs[t].accum = accum;
                threadArgs[t].currentSample = k;
                threadArgs[t].startTime = start;
                
                pthread_create(&threads[t], NULL, threadFunctionIterative, &threadArgs[t]);
                currentY = endY;
            }

            for (long t = 0; t < nThreads; ++t) {
                pthread_join(threads[t], NULL);
            }

            // Update canvas every 10th sample

            int samplesSoFar = k;
            if (samplesSoFar++ % 2 == 0) {
                writeIterativeCanvas(scene, camera, accum, width, height, samplesSoFar, image);
            }
        }

        // Final image after all samples
        writeIterativeCanvas(scene, camera, accum, width, height, numSamples, image);
        delete[] accum;
    }
    
    std::cout << std::endl;
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Multi-threaded ray tracing time: " << duration.count() << " milliseconds" << std::endl;
    
    printPerfStats();
    
    return duration.count();
}

double singleThreadedRayTracing(Scene& scene, Camera& camera, int width, int height, unsigned char* image, bool iterativeSampling) {
    auto start = std::chrono::high_resolution_clock::now();
    
    int numSamples = camera.numSamples;
    if (iterativeSampling) {
        g_totalPixels = width * height * numSamples;
    } else {
        g_totalPixels = width * height;
    }
    g_pixelsProcessed = 0;
    
    std::cout << "Rendering " << width << "x" << height << " (" << g_totalPixels << " pixels)..." << std::endl;
    
    if (!iterativeSampling) {
        // Original per-pixel sampling
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                VectorFloatTriplet pixelColor = {0.0, 0.0, 0.0};
                int pixelIndex = y * width + x;
                int sampleIndex = pixelIndex * camera.numSamples;
                for (int k = 0; k < camera.numSamples; k++) {
                    VectorFloatPenta sample = camera.samples[sampleIndex + k];
                    double sx = x + sample.x;   // subpixel x
                    double sy = y + sample.y;   // subpixel y
                    double sampleTime = sample.z;
                    double r1 = sample.w;
                    double r2 = sample.v;
                    pixelColor = pixelColor + __compute(scene, camera, sx, sy, width, height, sampleTime, r1, r2);
                }
                pixelColor = pixelColor * (1.0 / camera.numSamples);
                image[pixelIndex * 3 + 0] = (unsigned char) round(pixelColor.x);
                image[pixelIndex * 3 + 1] = (unsigned char) round(pixelColor.y);
                image[pixelIndex * 3 + 2] = (unsigned char) round(pixelColor.z);
                
                int processed = ++g_pixelsProcessed;
                if (processed % 10000 == 0 || processed == g_totalPixels) {
                    printProgress(processed, g_totalPixels, start);
                }
            }
        }
    } else {
        // Iterative, level-wise sampling
        VectorFloatTriplet* accum = new VectorFloatTriplet[width * height];
        for (int i = 0; i < width * height; ++i) {
            accum[i].x = accum[i].y = accum[i].z = 0.0;
        }

        for (int k = 0; k < numSamples; ++k) {
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    int pixelIndex = y * width + x;
                    int sampleIndexBase = pixelIndex * camera.numSamples;
                    VectorFloatPenta sample = camera.samples[sampleIndexBase + k];
                    double sx = x + sample.x;
                    double sy = y + sample.y;
                    double sampleTime = sample.z;
                    double r1 = sample.w;
                    double r2 = sample.v;
                    VectorFloatTriplet color =
                        __compute(scene, camera, sx, sy, width, height, sampleTime, r1, r2);
                    accum[pixelIndex].x += color.x;
                    accum[pixelIndex].y += color.y;
                    accum[pixelIndex].z += color.z;
                    
                    int processed = ++g_pixelsProcessed;
                    if (processed % 10000 == 0 || processed == g_totalPixels) {
                        printProgress(processed, g_totalPixels, start);
                    }
                }
            }

            int samplesSoFar = k + 1;
            if (samplesSoFar % 2 == 0) {
                writeIterativeCanvas(scene, camera, accum, width, height, samplesSoFar, image);
            }
        }

        // Final image
        writeIterativeCanvas(scene, camera, accum, width, height, numSamples, image);
        delete[] accum;
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
    args.iterativeSampling = false;

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
        else if (strcmp(argv[i], "-is") == 0 || strcmp(argv[i], "--iterative-sampling") == 0) {
            // Accept either -is 0/1 or just -is (default true)
            if (i + 1 < argc && (strcmp(argv[i + 1], "0") == 0 || strcmp(argv[i + 1], "1") == 0)) {
                args.iterativeSampling = (atoi(argv[i + 1]) != 0);
                i++;
            } else {
                args.iterativeSampling = true;
            }
        }
        else if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--gui") == 0) {
            // Accept either -g 0/1 or just -g (default true)
            if (i + 1 < argc && (strcmp(argv[i + 1], "0") == 0 || strcmp(argv[i + 1], "1") == 0)) {
                args.useGUI = (atoi(argv[i + 1]) != 0);
                i++;
            } else {
                args.useGUI = true;
            }
        }
    }
    return args;
}


int main(int argc, char* argv[])
{

    scene::Args args = parseArgs(argc, argv);

#if USE_GUI
    g_useGui = (args.iterativeSampling && args.useGUI);
#endif

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
        int numSamples = camera.numSamples;
        int width = camera.imageResolution.x;
        int height = camera.imageResolution.y;

        // Initialize GUI once per camera (before tracing starts), if requested.
#if USE_GUI
        if (g_useGui) {
            if (!initGui(width, height)) {
                std::cerr << "GUI initialization failed for camera " << i
                          << ", continuing without GUI." << std::endl;
                g_useGui = false;
            }
        }
#endif

        // numSamples is a perfect square (1, 4, 9, 16, etc.), total samples per pixel
        VectorFloatPenta* samples = new VectorFloatPenta[numSamples * width * height];
        precomputeSamples(numSamples, width, height, samples);
        camera.samples = samples;
        unsigned char* image = new unsigned char[width * height * 3];
        if (args.isMultiThreaded) {
            renderTimeMs = multiThreadedRayTracing(scene, camera, width, height, image, args.iterativeSampling);
            totalTimeMs += renderTimeMs;
        } else {
            renderTimeMs = singleThreadedRayTracing(scene, camera, width, height, image, args.iterativeSampling);
            totalTimeMs += renderTimeMs;
        }
        string outputName = camera.imageName;
        scene.writePPM((OUTPUT_PATH + outputName).c_str(), image, width, height);
        delete[] image;
        delete[] samples;

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
        std::ofstream resultsFile((JSON_OUTPUT_PATH + outputName.substr(0, outputName.find_last_of('.')) + "_results.json").c_str());
        resultsFile << results.dump(4);
        resultsFile.close();
    }

#if USE_GUI
    if (g_gui.initialized && g_useGui) {
        bool running = true;
        while (running) {
            SDL_Event e;
            while (SDL_PollEvent(&e)) {
                if (e.type == SDL_EVENT_QUIT) {
                    running = false;
                }
            }
            SDL_Delay(16);
        }
    }

    if (g_gui.initialized) {
        if (g_gui.texture) SDL_DestroyTexture(g_gui.texture);
        if (g_gui.renderer) SDL_DestroyRenderer(g_gui.renderer);
        if (g_gui.window) SDL_DestroyWindow(g_gui.window);
        SDL_Quit();
    }
#endif
    return 0;
}