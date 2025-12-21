#ifndef __TEXTURE__
#define __TEXTURE__

#include "scene.h"

namespace scene {
    // Forward declarations
    struct TextureMap;
    struct Image;
    struct Scene;
    struct Intersection;
    struct VectorFloatPair;
    struct VectorFloatTriplet;
    
    // Main texture sampling function
    VectorFloatTriplet sampleTexture(const TextureMap* textureMap,
                                     const VectorFloatPair& uv,
                                     const VectorFloatTriplet& position,
                                     const Scene* scene,
                                     bool applyNormalizer = true);
    
    // Texture type specific samplers
    VectorFloatTriplet sampleImageTexture(const TextureMap* textureMap,
                                          const VectorFloatPair& uv,
                                          const Scene* scene);
    
    VectorFloatTriplet samplePerlinNoise(const TextureMap* textureMap,
                                         const VectorFloatTriplet& position);
    
    VectorFloatTriplet sampleCheckerboard(const TextureMap* textureMap,
                                          const VectorFloatTriplet& position);
    
    // UV coordinate computation for different primitives
    VectorFloatPair computeUVCoordinates(const Intersection& intersection,
                                         const Scene& scene);
    
    // Interpolation methods
    VectorFloatTriplet bilinearInterpolation(const Image* image,
                                             double u, double v);
    
    VectorFloatTriplet trilinearInterpolation(const Image* image,
                                               double u, double v);
    
    // Perlin noise functions
    double perlinNoise(double x, double y, double z);
    double perlinNoise3D(double x, double y, double z);
    
    // Normal/Bump mapping support
    void computeTangentSpace(const VectorFloatTriplet& normal,
                             const VectorFloatTriplet& dpdu,
                             const VectorFloatTriplet& dpdv,
                             VectorFloatTriplet& tangent,
                             VectorFloatTriplet& bitangent);
    
    VectorFloatTriplet applyBumpMapping(const TextureMap* bumpMap,
                                        const VectorFloatPair& uv,
                                        const VectorFloatTriplet& position,
                                        const VectorFloatTriplet& geometricNormal,
                                        const Scene* scene,
                                        double bumpFactor);
    
    VectorFloatTriplet transformNormalFromTangentSpace(
        const VectorFloatTriplet& normalMapValue,
        const VectorFloatTriplet& tangent,
        const VectorFloatTriplet& bitangent,
        const VectorFloatTriplet& normal);
}

#endif

