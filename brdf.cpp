#include "brdf.h"
#include "utils.h"
#include "overloads.h"
#include <cmath>

using namespace scene;

namespace brdf {
    
    VectorFloatTriplet evaluateBRDF(const Material& mat, 
                                     const BRDF* brdf,
                                     const VectorFloatTriplet& N,
                                     const VectorFloatTriplet& wi,
                                     const VectorFloatTriplet& wo,
                                     VectorFloatTriplet& diffuse,
                                     VectorFloatTriplet& specular) {
        
        // Default to OriginalBlinnPhong if no BRDF specified
        BRDFType brdfType = BRDFType::OriginalBlinnPhong;
        double exponent = mat.phongExponent;
        bool normalized = false;
        
        if (brdf) {
            brdfType = brdf->type;
            exponent = brdf->exponent;
            normalized = brdf->normalized;
        }
        
        // Ensure directions are normalized
        VectorFloatTriplet wi_norm = normalize(wi);
        VectorFloatTriplet wo_norm = normalize(wo);
        VectorFloatTriplet N_norm = normalize(N);
        
        // Compute cosine terms
        double NdotL = std::max(0.0, dotProduct(N_norm, wi_norm));
        double NdotV = std::max(0.0, dotProduct(N_norm, wo_norm));
        
        // Diffuse component (Lambertian)
        diffuse = mat.diffuseReflectance;
        
        // Specular component depends on BRDF type
        specular = VectorFloatTriplet{0.0, 0.0, 0.0};
        
        if (NdotL <= 0.0 || NdotV <= 0.0) {
            return VectorFloatTriplet{0.0, 0.0, 0.0};
        }
        
        switch (brdfType) {
            case BRDFType::OriginalBlinnPhong: {
                // f_s = k_s * (N.H)^p
                VectorFloatTriplet H = normalize(wi_norm + wo_norm);
                double NdotH = std::max(0.0, dotProduct(N_norm, H));
                double specFactor = std::pow(NdotH, exponent);
                specular = mat.specularReflectance * specFactor;
                break;
            }
            
            case BRDFType::OriginalPhong: {
                // f_s = k_s * (R.V)^p
                VectorFloatTriplet R = wi_norm - N_norm * (2.0 * dotProduct(wi_norm, N_norm));
                double RdotV = std::max(0.0, dotProduct(R, wo_norm));
                double specFactor = std::pow(RdotV, exponent);
                specular = mat.specularReflectance * specFactor;
                break;
            }
            
            case BRDFType::ModifiedBlinnPhong: {
                // f_s = k_s * (N.H)^p * (N.L)
                VectorFloatTriplet H = normalize(wi_norm + wo_norm);
                double NdotH = std::max(0.0, dotProduct(N_norm, H));
                double specFactor = std::pow(NdotH, exponent) * NdotL;
                
                if (normalized) {
                    // Normalized: f_s = k_s * (p+8)/(8*pi) * (N.H)^p
                    double normalization = (exponent + 8.0) / (8.0 * M_PI);
                    specFactor = std::pow(NdotH, exponent) * normalization;
                }
                
                specular = mat.specularReflectance * specFactor;
                break;
            }
            
            case BRDFType::ModifiedPhong: {
                // f_s = k_s * (R.V)^p * (N.L)
                VectorFloatTriplet R = wi_norm - N_norm * (2.0 * dotProduct(wi_norm, N_norm));
                double RdotV = std::max(0.0, dotProduct(R, wo_norm));
                double specFactor = std::pow(RdotV, exponent) * NdotL;
                
                if (normalized) {
                    // Normalized: f_s = k_s * (p+2)/(2*pi) * (R.V)^p
                    double normalization = (exponent + 2.0) / (2.0 * M_PI);
                    specFactor = std::pow(RdotV, exponent) * normalization;
                }
                
                specular = mat.specularReflectance * specFactor;
                break;
            }
            
            case BRDFType::TorranceSparrow: {
                // Microfacet BRDF - simplified version
                // Full implementation would require D (distribution), G (geometry), F (Fresnel)
                // For now, use a simplified version
                VectorFloatTriplet H = normalize(wi_norm + wo_norm);
                double NdotH = std::max(0.0, dotProduct(N_norm, H));
                double VdotH = std::max(0.0, dotProduct(wo_norm, H));
                
                // D term: Beckmann distribution
                double alpha = std::sqrt(2.0 / (exponent + 2.0));  // Roughness parameter
                double alpha2 = alpha * alpha;
                double cosThetaH2 = NdotH * NdotH;
                double tanThetaH2 = (1.0 - cosThetaH2) / (cosThetaH2 + 1e-10);
                double D = std::exp(-tanThetaH2 / alpha2) / (M_PI * alpha2 * cosThetaH2 * cosThetaH2);
                
                // G term: simplified geometric term (Smith's model approximation)
                double G = 1.0;  // Simplified - full implementation would compute G1 for each direction
                
                // F term: Fresnel (simplified Schlick's approximation)
                VectorFloatTriplet F = mat.specularReflectance;  // Simplified
                
                // Specular: D * G * F / (4 * (N.L) * (N.V))
                double denominator = 4.0 * NdotL * NdotV;
                if (denominator > 1e-10) {
                    double invDenom = 1.0 / denominator;
                    specular = (F * D * G) * invDenom;
                }
                
                // Diffuse with optional Fresnel term
                if (brdf && brdf->kdFresnel) {
                    // Use (1-F)*kd/pi instead of kd/pi
                    VectorFloatTriplet oneMinusF = VectorFloatTriplet{1.0, 1.0, 1.0} - F;
                    double invPi = 1.0 / M_PI;
                    diffuse = (oneMinusF * mat.diffuseReflectance) * invPi;
                } else {
                    double invPi = 1.0 / M_PI;
                    diffuse = mat.diffuseReflectance * invPi;
                }
                break;
            }
        }
        
        // Return total BRDF value
        return diffuse + specular;
    }
    
    double pdfBRDF(const BRDF* brdf, 
                   const VectorFloatTriplet& N,
                   const VectorFloatTriplet& wi,
                   const VectorFloatTriplet& wo) {
        
        // For now, assume cosine-weighted hemisphere sampling
        // This will be used for MIS weighting
        VectorFloatTriplet N_norm = normalize(N);
        VectorFloatTriplet wi_norm = normalize(wi);
        
        double cosTheta = std::max(0.0, dotProduct(N_norm, wi_norm));
        
        // Cosine-weighted PDF: cos(theta) / pi
        return cosTheta / M_PI;
    }
}
