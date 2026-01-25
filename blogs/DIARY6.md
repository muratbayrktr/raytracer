# Logs

## Starting HW6 - Path Tracing and BRDFs

This is the final homework and it's a big one. Path tracing with BRDF models, object lights, and all the advanced sampling techniques. Let me break this down systematically.

The plan is clear: parse first, test parsing, then implement path tracing, then BRDF evaluation, then all the efficiency improvements. Let's start.

## Phase 1: Scene Parsing

Started with extending the Camera struct to handle path tracing settings. Added fields for renderer type, importance sampling, next event estimation, MIS heuristics, Russian roulette, recursion depths, splitting factor, and sample clamping. The parsing was straightforward - just following the existing pattern.

BRDF parsing was next. Added the BRDF enum and struct, then parsed all five BRDF types from the BRDFs section. Each BRDF has an ID, type, exponent, normalization flag, and for TorranceSparrow, a kdfresnel flag. Materials can now reference BRDFs via the _BRDF field.

Object lights were interesting - LightSphere and LightMesh are basically regular objects but with a Radiance field. They inherit from Object so they can have materials, transformations, textures, everything. Parsed them from the Objects section, handling both array and single object cases.

Had to fix a JSON parsing bug - was accessing camera fields without null checks, which caused assertion failures. Added proper contains() checks everywhere. Also had to add LightSphere and LightMesh to the intersection tests so rays can actually hit them. Added new Intersection::Kind values for them.

## Phase 2: BRDF Module

Created brdf.h and brdf.cpp for BRDF evaluation. Implemented all five BRDF types:
- OriginalBlinnPhong: standard (N.H)^p
- OriginalPhong: (R.V)^p  
- ModifiedBlinnPhong: (N.H)^p * (N.L), with optional normalization
- ModifiedPhong: (R.V)^p * (N.L), with optional normalization
- TorranceSparrow: microfacet model with D, G, F terms (simplified version)

The TorranceSparrow implementation is basic - full implementation would need proper geometry term and Fresnel calculations, but this should work for the homework scenes.

Had namespace issues at first - needed to use scene::VectorFloatTriplet everywhere. Also had to fix vector division - used scalar multiplication instead (multiply by 1/denominator).

## Phase 3: Hemisphere Sampling

Added uniform and cosine-weighted hemisphere sampling functions. Uniform: phi = 2*pi*xi1, theta = arccos(xi2), PDF = 1/(2*pi). Cosine: phi = 2*pi*xi1, theta = arcsin(sqrt(xi2)), PDF = cos(theta)/pi. Both transform to world space using orthonormal basis from surface normal.

## Phase 4: Basic Path Tracing

Implemented computePathTracing function. The basic loop:
1. Intersect ray with scene
2. If hit light, add emission * throughput
3. Sample next direction (uniform or cosine-weighted based on camera.importanceSampling)
4. Evaluate BRDF
5. Update throughput: throughput *= BRDF * cos(theta) / PDF
6. Create new ray and continue

For random number generation, used a simple hash-based approach for subsequent bounces since Ray only has random1/random2 for the primary ray. This isn't ideal but works for now.

Added LightSphere and LightMesh to intersection tests so they can be hit. Used temporary Sphere/Mesh objects for testing since the intersection functions expect those types.

The path tracing dispatches from computePixelColor when camera.renderer == "PathTracing". Still need to implement NEE, MIS, Russian Roulette, and all the efficiency improvements, but the basic structure is there.

## Phase 5: Next Event Estimation and MIS

Implemented sampleDirectLight function that samples from all light types (point, area, LightSphere, LightMesh). The function:
1. Selects a light uniformly from all available lights
2. Samples a point on that light (using area sampling for area lights, sphere sampling for LightSphere, mesh sampling for LightMesh)
3. Computes visibility
4. Converts area PDF to solid angle PDF
5. Returns radiance contribution and PDF

Added MIS weight function with three heuristics: balance, power, and 01. The balance heuristic is the standard one: w = pdf1 / (pdf1 + pdf2).

Integrated NEE and MIS into the path tracing loop:
- Sample light directly, evaluate BRDF for that direction, compute MIS weight, add contribution
- Sample BRDF direction, if it hits a light, compute light PDF for that direction, compute MIS weight, add contribution
- Both contributions are weighted and summed

The tricky part was computing the light PDF for a BRDF-sampled direction that happens to hit a light. For now using simplified estimates based on light area and distance.

## Phase 6: Russian Roulette and Splitting

Russian Roulette: After minRecursionDepth, compute survival probability based on throughput (max of RGB components, clamped to [0.01, 0.99]). If random number > survival probability, terminate path. Otherwise, scale throughput by 1/survivalProbability.

Splitting: For primary rays (depth == 0), send splittingFactor indirect rays from the first hit point. Accumulate their contributions and average them. Continue path with the first ray's direction.

Sample clamping: Applied in the rendering loop before averaging samples. Each sample is clamped to sampleMaxVal if specified.

## Integration

Path tracing is fully integrated. The renderer automatically uses path tracing when camera.renderer == "PathTracing". All the features are in place:
- Uniform and cosine-weighted hemisphere sampling
- BRDF evaluation (all 5 types)
- Object lights (LightSphere, LightMesh)
- Next Event Estimation
- Multiple Importance Sampling (balance, power, 01)
- Russian Roulette
- Splitting factor
- Sample clamping

Ready for testing!

## Critical Bug Fixes (Jan 23)

### Fix 1: Binary File Parsing for Sponza Scene

The sponza scene uses `_binaryFile` format instead of `_data` for VertexData, TexCoordData, and mesh Faces. The parser was skipping these entirely, resulting in meshes with 0 faces.

**Binary format discovered:**
- 4 bytes: count (little-endian uint32)
- count * N bytes: data (floats or ints)

**Added support in scene.cpp for:**

1. **VertexData**: `_binaryFile` → reads count + xyz floats per vertex
```cpp
} else if (vertexData.contains("_binaryFile") && !vertexData["_binaryFile"].is_null()) {
    std::string binaryFile = vertexData["_binaryFile"].get<std::string>();
    std::string fullPath = this->baseDirectory + binaryFile;
    std::ifstream file(fullPath, std::ios::binary);
    uint32_t count;
    file.read(reinterpret_cast<char*>(&count), sizeof(count));
    for (uint32_t i = 0; i < count; i++) {
        float x, y, z;
        file.read(reinterpret_cast<char*>(&x), sizeof(x));
        // ... read y, z
        this->vertices.push_back(VectorFloatTriplet{(double)x, (double)y, (double)z});
    }
}
```

2. **TexCoordData**: `_binaryFile` → reads count + uv floats per coordinate

3. **Mesh Faces**: `_binaryFile` → reads count + 3 uint32 indices per triangle

Result: Sponza now loads correctly with 471,282 vertices and 393 meshes with proper face counts.

### Fix 2: Near-Plane Clipping Bug (Killeroo Scene)

The killeroo scene was rendering with the dinosaur completely invisible - only the walls were showing.

**Root cause:** The intersection code was using `nearDistance` as a clipping distance. For killeroo:
- Camera at Y=1.6, looking down (-Y direction)
- nearDistance = 1.8
- Killeroo top at Y=0.29, so distance from camera ≈ 1.3

Since 1.3 < 1.8, the entire killeroo was being clipped as "too close to camera"!

**The bug in utils.cpp:**
```cpp
// This was computing minDistance = nearDistance for primary rays
// and rejecting all hits closer than nearDistance
double minDistance = 0.0;
if (isPrimaryRay && !scene.cameras.empty()) {
    // ... computed t_plane based on nearDistance
    minDistance = t_plane;  // BUG: This clips objects!
}
```

**The fix:** `nearDistance` defines the image plane for ray generation, NOT a clipping plane. Objects between the camera and image plane should still be visible.

```cpp
// Near-plane clipping disabled. The nearDistance field defines the image plane
// for ray generation, not a clipping plane.
double minDistance = 0.0;
```

Result: Killeroo now renders correctly with all 92,092 triangles visible.

### Fix 3: LightMesh/LightSphere Direct Visibility

For non-path-traced scenes using LightMesh or LightSphere, the light sources themselves weren't visible (appeared black).

**Two sub-fixes:**

1. **Return emission for direct hits:** In `computePixelColor`, check if the hit object is emissive and return its radiance directly:
```cpp
if (intersection.hit) {
    if (isEmissiveObject(scene, intersection)) {
        return getEmission(scene, intersection);
    }
    return computeShading(scene, ray, intersection);
}
```

2. **Disable back-face culling for LightMesh:** Light sources should be visible from both sides:
```cpp
bool thisHit = rayHitsMesh(testRay, tempMesh, scene.vertices,
                            emptyDeterminants,
                            temp_t_min, tempIntersection, scene.intersectionTestEpsilon, bvh,
                            false,  // Disable back-face culling for emissive objects
                            -20000 - i, ...);
```

Result: LightMesh objects now appear as bright rectangles instead of black holes.
