#ifndef __BRDF__
#define __BRDF__

#include "scene.h"

namespace brdf {
    // Evaluate BRDF f_r(wi, wo) given material, BRDF settings, and vectors
    // Returns the BRDF value and fills in diffuse and specular components
    // wi: incoming light direction (from light to surface)
    // wo: outgoing view direction (from surface to camera)
    // N: surface normal
    scene::VectorFloatTriplet evaluateBRDF(const scene::Material& mat, 
                                     const scene::BRDF* brdf,
                                     const scene::VectorFloatTriplet& N,
                                     const scene::VectorFloatTriplet& wi,  // incoming light dir
                                     const scene::VectorFloatTriplet& wo,  // outgoing view dir
                                     scene::VectorFloatTriplet& diffuse,
                                     scene::VectorFloatTriplet& specular);
    
    // PDF for sampling this BRDF direction (for MIS)
    // Returns the probability density of sampling direction wi given wo and N
    double pdfBRDF(const scene::BRDF* brdf, 
                   const scene::VectorFloatTriplet& N,
                   const scene::VectorFloatTriplet& wi,
                   const scene::VectorFloatTriplet& wo);
}

#endif
