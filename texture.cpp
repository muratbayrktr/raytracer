#include "texture.h"
#include "scene.h"
#include "utils.h"
#include "overloads.h"
#include <cmath>
#include <algorithm>

using namespace scene;

// Standard 512-entry permutation table for Perlin noise
static const int p[512] = {
    151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,
    8,99,37,240,21,10,23,190,6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,
    35,11,32,57,177,33,88,237,149,56,87,174,20,125,136,171,168,68,175,74,165,71,
    134,139,48,27,166,77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,
    55,46,245,40,244,102,143,54,65,25,63,161,1,216,80,73,209,76,132,187,208,89,
    18,169,200,196,135,130,116,188,159,86,164,100,109,198,173,186,3,64,52,217,226,
    250,124,123,5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,
    189,28,42,223,183,170,213,119,248,152,2,44,154,163,70,221,153,101,155,167,43,
    172,9,129,22,39,253,19,98,108,110,79,113,224,232,178,185,112,104,218,246,97,
    228,251,34,242,193,238,210,144,12,191,179,162,241,81,51,145,235,249,14,239,
    107,49,192,214,31,181,199,106,157,184,84,204,176,115,121,50,45,127,4,150,254,
    138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
    // repeat for wrapping
    151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,
    8,99,37,240,21,10,23,190,6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,
    35,11,32,57,177,33,88,237,149,56,87,174,20,125,136,171,168,68,175,74,165,71,
    134,139,48,27,166,77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,
    55,46,245,40,244,102,143,54,65,25,63,161,1,216,80,73,209,76,132,187,208,89,
    18,169,200,196,135,130,116,188,159,86,164,100,109,198,173,186,3,64,52,217,226,
    250,124,123,5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,
    189,28,42,223,183,170,213,119,248,152,2,44,154,163,70,221,153,101,155,167,43,
    172,9,129,22,39,253,19,98,108,110,79,113,224,232,178,185,112,104,218,246,97,
    228,251,34,242,193,238,210,144,12,191,179,162,241,81,51,145,235,249,14,239,
    107,49,192,214,31,181,199,106,157,184,84,204,176,115,121,50,45,127,4,150,254,
    138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
};

static double fade(double t) {
    return t * t * t * (t * (t * 6 - 15) + 10);
}

static double lerp(double a, double b, double t) {
    return a + t * (b - a);
}

static double grad(int hash, double x, double y, double z) {
    int h = hash & 15;
    double u = h < 8 ? x : y;
    double v = h < 4 ? y : (h == 12 || h == 14 ? x : z);
    return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
}

double scene::perlinNoise3D(double x, double y, double z) {
    int X = (int)floor(x) & 255;
    int Y = (int)floor(y) & 255;
    int Z = (int)floor(z) & 255;

    x -= floor(x);
    y -= floor(y);
    z -= floor(z);

    double u = fade(x);
    double v = fade(y);
    double w = fade(z);

    int A = p[X] + Y;
    int AA = p[A] + Z;
    int AB = p[A + 1] + Z;
    int B = p[X + 1] + Y;
    int BA = p[B] + Z;
    int BB = p[B + 1] + Z;

    return lerp(
        lerp(
            lerp(grad(p[AA], x, y, z),
                 grad(p[BA], x - 1, y, z), u),
            lerp(grad(p[AB], x, y - 1, z),
                 grad(p[BB], x - 1, y - 1, z), u), v),
        lerp(
            lerp(grad(p[AA + 1], x, y, z - 1),
                 grad(p[BA + 1], x - 1, y, z - 1), u),
            lerp(grad(p[AB + 1], x, y - 1, z - 1),
                 grad(p[BB + 1], x - 1, y - 1, z - 1), u), v), w);
}

double scene::perlinNoise(double x, double y, double z) {
    return perlinNoise3D(x, y, z);
}

VectorFloatTriplet scene::samplePerlinNoise(const TextureMap* textureMap,
                                            const VectorFloatTriplet& position) {
    if (!textureMap) {
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }

    double scaledX = position.x * textureMap->noiseScale;
    double scaledY = position.y * textureMap->noiseScale;
    double scaledZ = position.z * textureMap->noiseScale;

    double noiseValue = 0.0;
    double amplitude = 1.0;

    for (int k = 0; k < textureMap->numOctaves; k++) {
        double freq = pow(2.0, k);
        noiseValue += amplitude * perlinNoise3D(scaledX * freq, scaledY * freq, scaledZ * freq);
        amplitude *= 0.5;
    }

    if (textureMap->noiseConversion == NoiseConversion::AbsVal) {
        noiseValue = std::abs(noiseValue);
    } else {
        // Convert Perlin noise to [0,1]
        noiseValue = (noiseValue + 1.0) * 0.5;
        noiseValue = std::max(0.0, std::min(1.0, noiseValue));
    }

    return VectorFloatTriplet{noiseValue, noiseValue, noiseValue};
}

VectorFloatTriplet scene::sampleCheckerboard(const TextureMap* textureMap,
                                             const VectorFloatTriplet& position) {
    if (!textureMap) {
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }

    double x = (position.x + textureMap->offset.x) * textureMap->scale;
    double y = (position.y + textureMap->offset.y) * textureMap->scale;
    double z = (position.z + textureMap->offset.z) * textureMap->scale;

    bool xBool = ((int)floor(x)) % 2 != 0;
    bool yBool = ((int)floor(y)) % 2 != 0;
    bool zBool = ((int)floor(z)) % 2 != 0;

    bool xorXY = xBool != yBool;

    return (xorXY != zBool) ? textureMap->blackColor : textureMap->whiteColor;
}

VectorFloatTriplet scene::bilinearInterpolation(const Image* image,
                                                double u, double v) {
    if (!image || (!image->data && !image->hdrData)) {
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }

    u = std::max(0.0, std::min(1.0, u));
    v = std::max(0.0, std::min(1.0, v));

    double x = u * (image->width - 1);
    double y = v * (image->height - 1);

    int x0 = (int)floor(x);
    int y0 = (int)floor(y);
    int x1 = std::min(x0 + 1, image->width - 1);
    int y1 = std::min(y0 + 1, image->height - 1);

    double fx = x - x0;
    double fy = y - y0;

    int idx00 = (y0 * image->width + x0) * image->channels;
    int idx10 = (y0 * image->width + x1) * image->channels;
    int idx01 = (y1 * image->width + x0) * image->channels;
    int idx11 = (y1 * image->width + x1) * image->channels;

    VectorFloatTriplet c00, c10, c01, c11;

    if (image->isHDR && image->hdrData) {
        // HDR data is already in float format
        if (image->channels >= 3) {
            c00 = VectorFloatTriplet{image->hdrData[idx00], image->hdrData[idx00 + 1], image->hdrData[idx00 + 2]};
            c10 = VectorFloatTriplet{image->hdrData[idx10], image->hdrData[idx10 + 1], image->hdrData[idx10 + 2]};
            c01 = VectorFloatTriplet{image->hdrData[idx01], image->hdrData[idx01 + 1], image->hdrData[idx01 + 2]};
            c11 = VectorFloatTriplet{image->hdrData[idx11], image->hdrData[idx11 + 1], image->hdrData[idx11 + 2]};
        } else if (image->channels == 1) {
            double g00 = image->hdrData[idx00];
            double g10 = image->hdrData[idx10];
            double g01 = image->hdrData[idx01];
            double g11 = image->hdrData[idx11];
            c00 = VectorFloatTriplet{g00, g00, g00};
            c10 = VectorFloatTriplet{g10, g10, g10};
            c01 = VectorFloatTriplet{g01, g01, g01};
            c11 = VectorFloatTriplet{g11, g11, g11};
        } else {
            return VectorFloatTriplet{0.0, 0.0, 0.0};
        }
    } else if (image->data) {
        // LDR data needs normalization
        if (image->channels >= 3) {
            c00 = VectorFloatTriplet{image->data[idx00] / 255.0, image->data[idx00 + 1] / 255.0, image->data[idx00 + 2] / 255.0};
            c10 = VectorFloatTriplet{image->data[idx10] / 255.0, image->data[idx10 + 1] / 255.0, image->data[idx10 + 2] / 255.0};
            c01 = VectorFloatTriplet{image->data[idx01] / 255.0, image->data[idx01 + 1] / 255.0, image->data[idx01 + 2] / 255.0};
            c11 = VectorFloatTriplet{image->data[idx11] / 255.0, image->data[idx11 + 1] / 255.0, image->data[idx11 + 2] / 255.0};
        } else if (image->channels == 1) {
            double g00 = image->data[idx00] / 255.0;
            double g10 = image->data[idx10] / 255.0;
            double g01 = image->data[idx01] / 255.0;
            double g11 = image->data[idx11] / 255.0;
            c00 = VectorFloatTriplet{g00, g00, g00};
            c10 = VectorFloatTriplet{g10, g10, g10};
            c01 = VectorFloatTriplet{g01, g01, g01};
            c11 = VectorFloatTriplet{g11, g11, g11};
        } else {
            return VectorFloatTriplet{0.0, 0.0, 0.0};
        }
    } else {
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }

    VectorFloatTriplet c0 = c00 * (1.0 - fx) + c10 * fx;
    VectorFloatTriplet c1 = c01 * (1.0 - fx) + c11 * fx;
    VectorFloatTriplet result = c0 * (1.0 - fy) + c1 * fy;
    return result;
}

VectorFloatTriplet scene::trilinearInterpolation(const Image* image, double u, double v) {
    // Approximates trilinear by blending bilinear with a 2x-downsampled version
    if (!image || (!image->data && !image->hdrData)) {
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }

    double mipLevel = 0.5;
    VectorFloatTriplet color0 = bilinearInterpolation(image, u, v);

    int w2 = image->width / 2;
    int h2 = image->height / 2;
    if (w2 < 1) w2 = 1;
    if (h2 < 1) h2 = 1;

    double u2 = u * 0.5;
    double v2 = v * 0.5;
    int x = (int)(u2 * (w2 - 1));
    int y = (int)(v2 * (h2 - 1));
    x = std::max(0, std::min(w2 - 1, x));
    y = std::max(0, std::min(h2 - 1, y));

    int idx = ((y * 2) * image->width + (x * 2)) * image->channels;
    VectorFloatTriplet color1;
    if (image->channels >= 3 && idx + 2 < image->width * image->height * image->channels) {
        if (image->isHDR && image->hdrData) {
            color1 = VectorFloatTriplet{
                image->hdrData[idx],
                image->hdrData[idx + 1],
                image->hdrData[idx + 2]
            };
        } else if (image->data) {
            color1 = VectorFloatTriplet{
                image->data[idx] / 255.0,
                image->data[idx + 1] / 255.0,
                image->data[idx + 2] / 255.0
            };
        } else {
            color1 = color0;
        }
    } else {
        color1 = color0;
    }

    double blend = mipLevel - floor(mipLevel);
    return color0 * (1.0 - blend) + color1 * blend;
}

VectorFloatTriplet scene::sampleImageTexture(const TextureMap* textureMap,
                                             const VectorFloatPair& uv,
                                             const Scene* scene) {
    if (!textureMap || !scene) {
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }
    const Image* image = scene->getImageById(textureMap->imageId);
    if (!image || (!image->data && !image->hdrData)) {
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }

    double u = uv.x - floor(uv.x);
    double v = uv.y - floor(uv.y);

    // Clamp after wrapping to handle numerical errors
    u = std::max(0.0, std::min(1.0, u));
    v = std::max(0.0, std::min(1.0, v));

    VectorFloatTriplet result;
    switch (textureMap->interpolation) {
        case InterpolationMode::Nearest: {
            int x = (int)(u * (image->width - 1));
            int y = (int)(v * (image->height - 1));
            x = std::max(0, std::min(image->width - 1, x));
            y = std::max(0, std::min(image->height - 1, y));
            int idx = (y * image->width + x) * image->channels;
            if (image->isHDR && image->hdrData) {
                // HDR data is already in float format
                if (image->channels >= 3) {
                    result = VectorFloatTriplet{
                        image->hdrData[idx],
                        image->hdrData[idx + 1],
                        image->hdrData[idx + 2]
                    };
                } else if (image->channels == 1) {
                    double gray = image->hdrData[idx];
                    result = VectorFloatTriplet{gray, gray, gray};
                } else {
                    result = VectorFloatTriplet{0.0, 0.0, 0.0};
                }
            } else if (image->data) {
                // LDR data needs normalization
                if (image->channels >= 3) {
                    result = VectorFloatTriplet{
                        image->data[idx] / 255.0,
                        image->data[idx + 1] / 255.0,
                        image->data[idx + 2] / 255.0
                    };
                } else if (image->channels == 1) {
                    double gray = image->data[idx] / 255.0;
                    result = VectorFloatTriplet{gray, gray, gray};
                } else {
                    result = VectorFloatTriplet{0.0, 0.0, 0.0};
                }
            } else {
                result = VectorFloatTriplet{0.0, 0.0, 0.0};
            }
            break;
        }
        case InterpolationMode::Bilinear:
            result = bilinearInterpolation(image, u, v);
            break;
        case InterpolationMode::Trilinear:
            result = trilinearInterpolation(image, u, v);
            break;
        default:
            result = bilinearInterpolation(image, u, v);
            break;
    }
    return result;
}

VectorFloatPair scene::computeUVCoordinates(const Intersection& intersection,
                                            const Scene& scene) {
    VectorFloatPair uv{0.0, 0.0};
    switch (intersection.kind) {
        case Intersection::Kind::Triangle: {
            // Barycentric interpolation for triangles
            if (intersection.faceIndex >= 0 &&
                intersection.faceIndex < (int)scene.triangles.size()) {
                const Triangle& tri = scene.triangles[intersection.faceIndex];
                if (tri.indices.x < (int)scene.texCoords.size() &&
                    tri.indices.y < (int)scene.texCoords.size() &&
                    tri.indices.z < (int)scene.texCoords.size()) {
                    double alpha = 1.0 - intersection.beta - intersection.gamma;
                    uv.x = alpha * scene.texCoords[tri.indices.x].x +
                           intersection.beta * scene.texCoords[tri.indices.y].x +
                           intersection.gamma * scene.texCoords[tri.indices.z].x;
                    uv.y = alpha * scene.texCoords[tri.indices.x].y +
                           intersection.beta * scene.texCoords[tri.indices.y].y +
                           intersection.gamma * scene.texCoords[tri.indices.z].y;
                }
            }
            break;
        }
        case Intersection::Kind::Mesh: {
            int meshIdx = -1;
            if (intersection.containerIndex >= 0) {
                meshIdx = intersection.containerIndex;
            } else if (intersection.containerIndex == -2) {
                // Try matching mesh by material
                for (size_t i = 0; i < scene.meshes.size(); i++) {
                    if (intersection.faceIndex >= 0 &&
                        intersection.faceIndex < (int)scene.meshes[i].faces.size() &&
                        scene.meshes[i].material == intersection.material) {
                        meshIdx = i;
                        break;
                    }
                }
                // Fallback: match by face index
                if (meshIdx == -1) {
                    for (size_t i = 0; i < scene.meshes.size(); i++) {
                        if (intersection.faceIndex >= 0 &&
                            intersection.faceIndex < (int)scene.meshes[i].faces.size()) {
                            meshIdx = i;
                            break;
                        }
                    }
                }
            }
            if (meshIdx >= 0 && meshIdx < (int)scene.meshes.size() &&
                intersection.faceIndex >= 0) {
                const Mesh& mesh = scene.meshes[meshIdx];
                if (intersection.faceIndex < (int)mesh.faces.size()) {
                    const VectorIntTriplet& texFace = mesh.texCoordIndices.empty()
                        ? mesh.faces[intersection.faceIndex]
                        : mesh.texCoordIndices[intersection.faceIndex];

                    if (texFace.x >= 0 && texFace.y >= 0 && texFace.z >= 0 &&
                        texFace.x < (int)scene.texCoords.size() &&
                        texFace.y < (int)scene.texCoords.size() &&
                        texFace.z < (int)scene.texCoords.size()) {
                        double alpha = 1.0 - intersection.beta - intersection.gamma;
                        uv.x = alpha * scene.texCoords[texFace.x].x +
                               intersection.beta * scene.texCoords[texFace.y].x +
                               intersection.gamma * scene.texCoords[texFace.z].x;
                        uv.y = alpha * scene.texCoords[texFace.x].y +
                               intersection.beta * scene.texCoords[texFace.y].y +
                               intersection.gamma * scene.texCoords[texFace.z].y;
                    }
                }
            }
            break;
        }
        case Intersection::Kind::Sphere: {
            // Spherical mapping
            VectorFloatTriplet normal = normalize(intersection.geometricNormal);
            double u = atan2(normal.z, normal.x) / (2.0 * M_PI) + 0.5;
            double v = acos(normal.y) / M_PI;
            uv.x = u;
            uv.y = v;
            break;
        }
        case Intersection::Kind::Plane: {
            // Planar mapping (uses x and z)
            VectorFloatTriplet point = intersection.point;
            uv.x = point.x * 0.1;
            uv.y = point.z * 0.1;
            uv.x = uv.x - floor(uv.x);
            uv.y = uv.y - floor(uv.y);
            break;
        }
        default:
            uv = VectorFloatPair{0.0, 0.0};
            break;
    }
    return uv;
}

void scene::computeTangentSpace(
    const VectorFloatTriplet& normal,
    const VectorFloatTriplet& dpdu,
    const VectorFloatTriplet& dpdv,
    VectorFloatTriplet& tangent,
    VectorFloatTriplet& bitangent) {
    tangent = normalize(dpdu);
    double dot = dotProduct(tangent, normal);
    tangent = normalize(tangent - normal * dot);
    bitangent = crossProduct(normal, tangent);
}

VectorFloatTriplet scene::transformNormalFromTangentSpace(
    const VectorFloatTriplet& normalMapValue,
    const VectorFloatTriplet& tangent,
    const VectorFloatTriplet& bitangent,
    const VectorFloatTriplet& normal) {
    // Convert from [0,1] to [-1,1], then from tangent to world space
    VectorFloatTriplet n = (normalMapValue * 2.0) - VectorFloatTriplet{1.0, 1.0, 1.0};
    n = normalize(n);
    VectorFloatTriplet result;
    result.x = tangent.x * n.x + bitangent.x * n.y + normal.x * n.z;
    result.y = tangent.y * n.x + bitangent.y * n.y + normal.y * n.z;
    result.z = tangent.z * n.x + bitangent.z * n.y + normal.z * n.z;
    return normalize(result);
}

VectorFloatTriplet scene::applyBumpMapping(const TextureMap* bumpMap,
                                           const VectorFloatPair& uv,
                                           const VectorFloatTriplet& position,
                                           const VectorFloatTriplet& geometricNormal,
                                           const Scene* scene,
                                           double bumpFactor) {
    if (!bumpMap || !scene)
        return geometricNormal;

    VectorFloatTriplet tangent, bitangent;
    VectorFloatTriplet up{0.0, 1.0, 0.0};
    VectorFloatTriplet right{1.0, 0.0, 0.0};

    if (std::abs(dotProduct(geometricNormal, up)) < 0.999)
        tangent = normalize(crossProduct(geometricNormal, up));
    else
        tangent = normalize(crossProduct(geometricNormal, right));
    bitangent = normalize(crossProduct(geometricNormal, tangent));

    VectorFloatTriplet height = sampleTexture(bumpMap, uv, position, scene, false);
    double h = (height.x + height.y + height.z) / 3.0;

    double du = 0.0, dv = 0.0;

    if (bumpMap->type == "image") {
        // Compute bump map gradient in UV space (image textures)
        double epsilon = 1.0 / 512.0;

        VectorFloatPair uvU = uv; uvU.x += epsilon; uvU.x = uvU.x - floor(uvU.x);
        VectorFloatPair uvV = uv; uvV.y += epsilon; uvV.y = uvV.y - floor(uvV.y);

        VectorFloatTriplet hU = sampleTexture(bumpMap, uvU, position, scene, false);
        VectorFloatTriplet hV = sampleTexture(bumpMap, uvV, position, scene, false);
        double hu = (hU.x + hU.y + hU.z) / 3.0;
        double hv = (hV.x + hV.y + hV.z) / 3.0;

        du = (hu - h) / epsilon;
        dv = (hv - h) / epsilon;
    } else {
        // Compute gradient in position space (procedural textures)
        double eps = 0.001;
        VectorFloatTriplet posU = position + tangent * eps;
        VectorFloatTriplet posV = position + bitangent * eps;
        VectorFloatTriplet hU = sampleTexture(bumpMap, uv, posU, scene, false);
        VectorFloatTriplet hV = sampleTexture(bumpMap, uv, posV, scene, false);
        double hu = (hU.x + hU.y + hU.z) / 3.0;
        double hv = (hV.x + hV.y + hV.z) / 3.0;
        du = (hu - h) / eps;
        dv = (hv - h) / eps;
    }

    VectorFloatTriplet perturbed = geometricNormal - tangent * (du * bumpFactor) - bitangent * (dv * bumpFactor);
    return normalize(perturbed);
}

VectorFloatTriplet scene::sampleTexture(const TextureMap* textureMap,
                                        const VectorFloatPair& uv,
                                        const VectorFloatTriplet& position,
                                        const Scene* scene,
                                        bool applyNormalizer) {
    if (!textureMap || !scene)
        return VectorFloatTriplet{0.0, 0.0, 0.0};

    VectorFloatTriplet textureValue;
    if (textureMap->type == "image") {
        textureValue = sampleImageTexture(textureMap, uv, scene);
        if (applyNormalizer && textureMap->normalizer > 0.0) {
            // Check if the image is HDR - HDR values are already in proper radiance units
            const Image* image = scene->getImageById(textureMap->imageId);
            if (image && image->isHDR) {
                // For HDR textures, values are already in proper radiance units
                // Only apply normalizer if it's not 1 (default), to allow user scaling
                if (textureMap->normalizer != 1.0) {
                    textureValue = textureValue * (1.0 / textureMap->normalizer);
                }
                // Otherwise, keep HDR values as-is
            } else {
                // For LDR textures, scale from [0,1] to [0,255/normalizer]
                double normalizerCorrection = 255.0 / textureMap->normalizer;
                textureValue = textureValue * normalizerCorrection;
            }
        }
    } else if (textureMap->type == "perlin") {
        textureValue = samplePerlinNoise(textureMap, position);
    } else if (textureMap->type == "checkerboard") {
        textureValue = sampleCheckerboard(textureMap, position);
    } else {
        return VectorFloatTriplet{0.0, 0.0, 0.0};
    }
    return textureValue;
}
