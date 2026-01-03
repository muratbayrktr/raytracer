#include "hdr_io.h"
#define TINYEXR_USE_STB_ZLIB 1
#define TINYEXR_USE_MINIZ 0
#include "tinyexr.h"
#include "stb_image_write.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>
#include <chrono>

namespace hdr_io {

bool isHDRFile(const std::string& filename) {
    std::string lowerFilename = filename;
    std::transform(lowerFilename.begin(), lowerFilename.end(), lowerFilename.begin(), ::tolower);
    return (lowerFilename.length() >= 4 && lowerFilename.substr(lowerFilename.length() - 4) == ".exr") ||
           (lowerFilename.length() >= 4 && lowerFilename.substr(lowerFilename.length() - 4) == ".hdr");
}

bool writeEXR(const std::string& filename, const float* rgb, int width, int height) {
    // Convert RGB to RGBA (EXR expects RGBA)
    float* rgba = new float[width * height * 4];
    for (int i = 0; i < width * height; ++i) {
        rgba[i * 4 + 0] = rgb[i * 3 + 0];
        rgba[i * 4 + 1] = rgb[i * 3 + 1];
        rgba[i * 4 + 2] = rgb[i * 3 + 2];
        rgba[i * 4 + 3] = 1.0f;  // Alpha
    }
    
    const char* err = nullptr;
    int ret = SaveEXR(rgba, width, height, 4, 0, filename.c_str(), &err);
    
    delete[] rgba;
    
    if (ret != TINYEXR_SUCCESS) {
        if (err) {
            // Error message is already set by tinyexr
            FreeEXRErrorMessage(err);
        }
        return false;
    }
    
    return true;
}

bool writeHDR(const std::string& filename, const float* rgb, int width, int height) {
    int ret = stbi_write_hdr(filename.c_str(), width, height, 3, rgb);
    return ret != 0;
}

}

