#ifndef __HDR_IO__
#define __HDR_IO__

#include <string>

namespace hdr_io {
    // Check if filename is an HDR file (.exr or .hdr)
    bool isHDRFile(const std::string& filename);
    
    // Write EXR file using tinyexr
    bool writeEXR(const std::string& filename, const float* rgb, int width, int height);
    
    // Write HDR file using stbi_write_hdr
    bool writeHDR(const std::string& filename, const float* rgb, int width, int height);
}

#endif

