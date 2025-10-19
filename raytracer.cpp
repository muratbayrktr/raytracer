
#include <iostream>
#include "scene.h"


void printScene(scene::Scene scene) {
    std::cout << "Scene:" << std::endl;
    std::cout << "BackgroundColor: " << scene.backgroundColor.x << " " << scene.backgroundColor.y << " " << scene.backgroundColor.z << std::endl;
    std::cout << "ShadowRayEpsilon: " << scene.shadowRayEpsilon << std::endl;
    std::cout << "IntersectionTestEpsilon: " << scene.intersectionTestEpsilon << std::endl;
    std::cout << "Cameras: " << scene.cameras.size() << std::endl;
    std::cout << "Lights: " << scene.pointLights.size() << std::endl;
    std::cout << "Materials: " << scene.materials.size() << std::endl;
    std::cout << "VertexData: " << scene.vertexData.data.size() << std::endl;
    std::cout << "Meshes: " << scene.meshes.size() << std::endl;
    std::cout << "Triangles: " << scene.triangles.size() << std::endl;
    std::cout << "Spheres: " << scene.spheres.size() << std::endl;
    std::cout << "Planes: " << scene.planes.size() << std::endl;
}

int main(int argc, char* argv[])
{
    scene::Scene scene;
    scene.loadSceneFromFile(argv[1]);

    // Render the scene
    printScene(scene);

    return 0;
}
