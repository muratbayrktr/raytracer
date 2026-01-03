#include "tonemap.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>
#include <fstream>
#include <map>
#include <chrono>

namespace tonemap {

double computeLuminance(const scene::VectorFloatTriplet& rgb) {
    return 0.2126 * rgb.x + 0.7152 * rgb.y + 0.0722 * rgb.z;
}

scene::VectorFloatTriplet photographic(const scene::VectorFloatTriplet& rgb, double luminanceScale, double burnOutPercent, const std::vector<double>& sortedLuminances) {
    double Yi = computeLuminance(rgb);
    if (Yi <= 0.0) return scene::VectorFloatTriplet{0.0, 0.0, 0.0};
    
    // Scale luminance by the precomputed luminanceScale (which is key / Lw)
    // This maps the scene's log-average luminance to the key value
    double L = luminanceScale * Yi;
    
    // Reinhard's simple operator: Ld = L / (1 + L)
    double Ld = L / (1.0 + L);
    
    // Handle burn-out if specified (extended Reinhard operator)
    if (burnOutPercent > 0.0 && !sortedLuminances.empty()) {
        int percentileIndex = static_cast<int>((1.0 - burnOutPercent / 100.0) * sortedLuminances.size());
        percentileIndex = std::max(0, std::min(percentileIndex, static_cast<int>(sortedLuminances.size() - 1)));
        double LwhiteValue = sortedLuminances[percentileIndex];
        if (LwhiteValue > 0.0) {
            double Lwhite = luminanceScale * LwhiteValue;
            // Extended Reinhard: Ld = L * (1 + L/Lwhite^2) / (1 + L)
            Ld = L * (1.0 + L / (Lwhite * Lwhite)) / (1.0 + L);
        }
    }
    
    // Scale RGB proportionally to preserve color
    double scale = Ld / Yi;
    return scene::VectorFloatTriplet{rgb.x * scale, rgb.y * scale, rgb.z * scale};
}

scene::VectorFloatTriplet filmic(const scene::VectorFloatTriplet& rgb, double luminanceScale, double burnOutPercent, const std::vector<double>& sortedLuminances) {
    // Filmic tone mapping using Uncharted 2 filmic curve
    double Yi = computeLuminance(rgb);
    if (Yi <= 0.0) return scene::VectorFloatTriplet{0.0, 0.0, 0.0};
    
    double L = luminanceScale * Yi;
    
    // Filmic curve: (x * (a * x + b)) / (x * (c * x + d) + e)
    double a = 2.51;
    double b = 0.03;
    double c = 2.43;
    double d = 0.59;
    double e = 0.14;
    
    double Ld = (L * (a * L + b)) / (L * (c * L + d) + e);
    Ld = std::max(0.0, std::min(1.0, Ld));  // Clamp to [0,1]
    
    double scale = Ld / Yi;
    return scene::VectorFloatTriplet{rgb.x * scale, rgb.y * scale, rgb.z * scale};
}

scene::VectorFloatTriplet aces(const scene::VectorFloatTriplet& rgb, double luminanceScale, double burnOutPercent, const std::vector<double>& sortedLuminances) {
    // ACES tone mapping approximation (Stephen Hill's fit)
    double Yi = computeLuminance(rgb);
    if (Yi <= 0.0) return scene::VectorFloatTriplet{0.0, 0.0, 0.0};
    
    double L = luminanceScale * Yi;
    
    // ACES curve approximation (same as filmic in this implementation)
    double a = 2.51;
    double b = 0.03;
    double c = 2.43;
    double d = 0.59;
    double e = 0.14;
    
    double Ld = (L * (a * L + b)) / (L * (c * L + d) + e);
    Ld = std::max(0.0, std::min(1.0, Ld));  // Clamp to [0,1]
    
    double scale = Ld / Yi;
    return scene::VectorFloatTriplet{rgb.x * scale, rgb.y * scale, rgb.z * scale};
}

scene::VectorFloatTriplet applySaturation(const scene::VectorFloatTriplet& rgb, const scene::VectorFloatTriplet& originalRgb, double saturation) {
    double Yi = computeLuminance(originalRgb);
    if (Yi <= 0.0) return rgb;
    
    // Equation from spec: Ro = Yo * (R/Yi)^s
    double s = saturation;
    double Yo = computeLuminance(rgb);
    
    double R = Yo * std::pow(std::max(0.0, originalRgb.x / Yi), s);
    double G = Yo * std::pow(std::max(0.0, originalRgb.y / Yi), s);
    double B = Yo * std::pow(std::max(0.0, originalRgb.z / Yi), s);
    
    return scene::VectorFloatTriplet{R, G, B};
}

void applyGammaCorrection(scene::VectorFloatTriplet& rgb, double gamma) {
    double invGamma = 1.0 / gamma;
    rgb.x = std::pow(std::max(0.0, std::min(1.0, rgb.x)), invGamma);
    rgb.y = std::pow(std::max(0.0, std::min(1.0, rgb.y)), invGamma);
    rgb.z = std::pow(std::max(0.0, std::min(1.0, rgb.z)), invGamma);
}

void toneMapImage(float* hdrImage, unsigned char* ldrImage, int width, int height, const scene::TonemapSettings& settings) {
    // Parse TMOOptions to get key and burnOutPercent
    std::istringstream iss(settings.tmoOptions);
    double key = 0.18;  // Default
    double burnOutPercent = 0.0;  // Default
    iss >> key;
    if (!iss.eof()) {
        iss >> burnOutPercent;
    }
    
    // Collect all luminances and compute log-average luminance (Reinhard's Lw)
    std::vector<double> luminances;
    luminances.reserve(width * height);
    double logSum = 0.0;
    int validPixels = 0;
    const double delta = 1e-6;  // Small value to avoid log(0)
    
    for (int i = 0; i < width * height; ++i) {
        scene::VectorFloatTriplet rgb{hdrImage[i * 3 + 0], hdrImage[i * 3 + 1], hdrImage[i * 3 + 2]};
        double lum = computeLuminance(rgb);
        if (lum > 0.0) {
            luminances.push_back(lum);
            logSum += std::log(delta + lum);
            validPixels++;
        }
    }
    
    // Compute log-average luminance (Lw) - Equation 1 from Reinhard
    double Lw = 1.0;  // Default fallback
    if (validPixels > 0) {
        Lw = std::exp(logSum / validPixels);
    }
    
    // Compute the luminance scaling factor: key / Lw
    // This normalizes the scene so that average luminance maps to 'key'
    double luminanceScale = (Lw > 0.0) ? (key / Lw) : key;
    
    // Sort luminances for percentile calculation
    std::vector<double> sortedLuminances = luminances;
    std::sort(sortedLuminances.begin(), sortedLuminances.end());
    
    // Apply tone mapping to each pixel
    for (int i = 0; i < width * height; ++i) {
        scene::VectorFloatTriplet originalRgb{hdrImage[i * 3 + 0], hdrImage[i * 3 + 1], hdrImage[i * 3 + 2]};
        scene::VectorFloatTriplet rgb = originalRgb;
        
        // Apply tone mapping operator
        if (settings.tmo == "Photographic") {
            rgb = photographic(originalRgb, luminanceScale, burnOutPercent, sortedLuminances);
        } else if (settings.tmo == "Filmic") {
            rgb = filmic(originalRgb, luminanceScale, burnOutPercent, sortedLuminances);
        } else if (settings.tmo == "ACES") {
            rgb = aces(originalRgb, luminanceScale, burnOutPercent, sortedLuminances);
        }
        
        // Apply saturation
        rgb = applySaturation(rgb, originalRgb, settings.saturation);
        
        // Clamp to [0,1]
        rgb.x = std::max(0.0, std::min(1.0, rgb.x));
        rgb.y = std::max(0.0, std::min(1.0, rgb.y));
        rgb.z = std::max(0.0, std::min(1.0, rgb.z));
        
        // Apply gamma correction
        applyGammaCorrection(rgb, settings.gamma);
        
        // Convert to 0-255 range
        ldrImage[i * 3 + 0] = static_cast<unsigned char>(std::round(rgb.x * 255.0));
        ldrImage[i * 3 + 1] = static_cast<unsigned char>(std::round(rgb.y * 255.0));
        ldrImage[i * 3 + 2] = static_cast<unsigned char>(std::round(rgb.z * 255.0));
    }
}

}


