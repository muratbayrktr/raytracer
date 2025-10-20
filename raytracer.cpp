
#include <iostream>
#include "scene.h"


int main(int argc, char* argv[])
{
    scene::Scene scene;
    scene.loadSceneFromFile(argv[1]);

    // Render the scene
    scene.getSummary();

    return 0;
}
