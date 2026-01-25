# My Raytracer Journey – Fall 2025 CENG 795 HW6 Blog

> **Murat Bayraktar – 2448199**

This is HW6 — path tracing with BRDFs, object lights, and all the advanced sampling techniques. I thought hw5 would be the last one due to time constraints but it looks like I had to do one more. The bugs, the fixes, the chaos, and the screenshots. This time, I had to be more selective, due to time constraints I couldn't handle all the cases perfectly, but I got the core ideas working.

# TL;DR

### What I implemented
- Path tracing with BRDF evaluation (all 5 types)
- Object lights (LightSphere, LightMesh)
- Hemisphere sampling (uniform and cosine-weighted)
- Next Event Estimation (NEE)
- Multiple Importance Sampling (MIS) with balance heuristic
- Russian Roulette
- Splitting factor
- Sample clamping

### What I fixed
- Binary file parsing for Sponza scene (was loading 0 faces)
- Near-plane clipping bug (was clipping killeroo completely)
- LightMesh/LightSphere direct visibility (were appearing black)

### What still needs work
- Some cornellbox scenes don't match ground truth perfectly
- TorranceSparrow BRDF is simplified (needs proper geometry term)
- Some advanced MIS heuristics not fully implemented
- Time constraints meant I couldn't polish everything

# The Path Tracing Journey

HW6 started with a fundamental shift: moving from Whitted-style ray tracing to path tracing. For five homeworks, my renderer had been living in a world where light bounced once or twice, and then we gave up. Now we're entering the real world where light bounces forever (or until Russian Roulette kills it), and I had to convince my code that this was okay.

The basic path tracing loop is straightforward:
1. Intersect ray with scene
2. If hit light, add emission * throughput
3. Sample next direction (uniform or cosine-weighted)
4. Evaluate BRDF
5. Update throughput: throughput *= BRDF * cos(theta) / PDF
6. Create new ray and continue

But then you add NEE, MIS, Russian Roulette, splitting, clamping... and suddenly it's not so straightforward anymore.

# BRDF Evaluation: The Five Types

I implemented all five BRDF types:
- **OriginalBlinnPhong**: standard (N.H)^p
- **OriginalPhong**: (R.V)^p  
- **ModifiedBlinnPhong**: (N.H)^p * (N.L), with optional normalization
- **ModifiedPhong**: (R.V)^p * (N.L), with optional normalization
- **TorranceSparrow**: microfacet model with D, G, F terms (simplified version)

The TorranceSparrow implementation is basic — full implementation would need proper geometry term and Fresnel calculations, but this should work for the homework scenes. I got the idea, but time constraints meant I couldn't fully polish it.

# The Binary File Parsing Nightmare

The Sponza scene was rendering completely empty — no geometry at all. Turns out Sponza uses `_binaryFile` format instead of `_data` for VertexData, TexCoordData, and mesh Faces. The parser was skipping these entirely, resulting in meshes with 0 faces.

**Binary format discovered:**
- 4 bytes: count (little-endian uint32)
- count * N bytes: data (floats or ints)

I added support for reading binary files in `scene.cpp`:
- **VertexData**: `_binaryFile` : reads count + xyz floats per vertex
- **TexCoordData**: `_binaryFile` : reads count + uv floats per coordinate
- **Mesh Faces**: `_binaryFile` : reads count + 3 uint32 indices per triangle

Result: Sponza now loads correctly with 471,282 vertices and 393 meshes with proper face counts. But the rendering still has issues — more on that later.

# The Near-Plane Clipping Bug

The killeroo scene was rendering with the dinosaur completely invisible — only the walls were showing. Not "dark" — completely invisible. Like the dinosaur had been Thanos-snapped out of existence.

**Root cause:** The intersection code was using `nearDistance` as a clipping distance. For killeroo:
- Camera at Y=1.6, looking down (-Y direction)
- nearDistance = 1.8
- Killeroo top at Y=0.29, so distance from camera ≈ 1.3

Since 1.3 < 1.8, the entire killeroo was being clipped as "too close to camera"!

**The fix:** `nearDistance` defines the image plane for ray generation, NOT a clipping plane. Objects between the camera and image plane should still be visible. I removed the clipping logic, and killeroo now renders correctly with all 92,092 triangles visible.

# Object Lights: When Lights Are Also Objects

Object lights (LightSphere and LightMesh) are interesting — they're basically regular objects but with a Radiance field. They inherit from Object so they can have materials, transformations, textures, everything.

The tricky part was making them visible. For non-path-traced scenes using LightMesh or LightSphere, the light sources themselves weren't visible (appeared black).

**Two sub-fixes:**

1. **Return emission for direct hits:** In `computePixelColor`, check if the hit object is emissive and return its radiance directly.

2. **Disable back-face culling for LightMesh:** Light sources should be visible from both sides.

Result: LightMesh objects now appear as bright rectangles instead of black holes.

# Next Event Estimation and MIS

Next Event Estimation (NEE) samples lights directly, which reduces variance significantly. The function:
1. Selects a light uniformly from all available lights
2. Samples a point on that light (using area sampling for area lights, sphere sampling for LightSphere, mesh sampling for LightMesh)
3. Computes visibility
4. Converts area PDF to solid angle PDF
5. Returns radiance contribution and PDF

Multiple Importance Sampling (MIS) combines light sampling and BRDF sampling. I implemented the balance heuristic: w = pdf1 / (pdf1 + pdf2). The tricky part was computing the light PDF for a BRDF-sampled direction that happens to hit a light. For now using simplified estimates based on light area and distance.

# Russian Roulette and Splitting

Russian Roulette: After minRecursionDepth, compute survival probability based on throughput (max of RGB components, clamped to [0.01, 0.99]). If random number > survival probability, terminate path. Otherwise, scale throughput by 1/survivalProbability.

Splitting: For primary rays (depth == 0), send splittingFactor indirect rays from the first hit point. Accumulate their contributions and average them. Continue path with the first ray's direction.

Sample clamping: Applied in the rendering loop before averaging samples. Each sample is clamped to sampleMaxVal if specified.

# Results Gallery

Some of my renders from HW6:

| Scene | Image |
| --- | --- |
| Diffuse Cornell Box (Default) | ![](./outputs_hw6/diffuse_cornell_box_default_phot.png) |
| Diffuse Cornell Box (Importance + NEE + MIS) | ![](./outputs_hw6/diffuse_cornell_box_importance_nee_mis_balance_phot.png) |
| Cornell Box Jaroslav (Glossy) | ![](./outputs_hw6/cornellbox_jaroslav_glossy_phot.png) |

The diffuse cornell box scenes work reasonably well. The basic path tracing produces noisy but correct results, and adding importance sampling, NEE, and MIS reduces the noise significantly.

# What Didn't Work (Time Constraints)

Unfortunately, time constraints meant I couldn't handle all cases perfectly. Here are some scenes that still have issues:

| Scene | Issue |
| --- | --- |
| Killeroo BlinnPhong | ![](./outputs_hw6/killeroo_blinnphong_phot.png) |
| Sponza Path | ![](./outputs_hw6/sponza_path_phot.png) |

The killeroo scene renders, but the lighting doesn't quite match the ground truth. The BRDF evaluation might be off, or there could be issues with the sampling. The Sponza scene loads correctly now, but the path tracing produces results that don't match the expected output. I tried something, got the idea somehow, but couldn't completely polish it and it has some issues still.

Some of the more advanced cornellbox scenes (like the ones with sphere lights or prism lights) also don't match perfectly. The core path tracing works, but the edge cases and optimizations need more work.

# Final Thoughts

This was a challenging homework. Path tracing is conceptually simple but getting all the details right is hard. I got the core ideas working:
- Basic path tracing loop
- BRDF evaluation
- Object lights
- NEE and MIS
- Russian Roulette and splitting

But time constraints meant I couldn't polish everything. Some scenes work well, some don't. The diffuse cornell box scenes are close to ground truth, but the more complex scenes (killeroo, sponza) still have issues.

Looking back, I'm happy I got the basic path tracing working. The renderer can now handle global illumination, which is a huge step forward. But there's still work to be done on the advanced features and edge cases.

It's been a wild ride.
