#ifndef __TONEMAP__
#define __TONEMAP__

#include "scene.h"
#include <vector>
#include <string>

namespace tonemap {
    // Compute luminance from RGB: Y = 0.2126*R + 0.7152*G + 0.0722*B
    double computeLuminance(const scene::VectorFloatTriplet& rgb);
    
    // Tone mapping operators
    scene::VectorFloatTriplet photographic(const scene::VectorFloatTriplet& rgb, double key, double burnOutPercent, const std::vector<double>& sortedLuminances);
    scene::VectorFloatTriplet filmic(const scene::VectorFloatTriplet& rgb, double key, double burnOutPercent, const std::vector<double>& sortedLuminances);
    scene::VectorFloatTriplet aces(const scene::VectorFloatTriplet& rgb, double key, double burnOutPercent, const std::vector<double>& sortedLuminances);
    
    // Apply saturation to RGB based on original RGB and compressed luminance
    scene::VectorFloatTriplet applySaturation(const scene::VectorFloatTriplet& rgb, const scene::VectorFloatTriplet& originalRgb, double saturation);
    
    // Apply gamma correction
    void applyGammaCorrection(scene::VectorFloatTriplet& rgb, double gamma);
    
    // Main tone mapping function that processes entire image
    void toneMapImage(float* hdrImage, unsigned char* ldrImage, int width, int height, const scene::TonemapSettings& settings);
}

#endif

