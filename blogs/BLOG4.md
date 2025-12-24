# My Raytracer Journey – Fall 2025 CENG 795 HW4 Blog

> **Murat Bayraktar – 2448199**

This is HW4 — the one where I thought "how hard can texturing be?" and then spent two days discovering that the answer is "very." Turns out wrapping images around 3D objects involves approximately 47 different ways to mess up color spaces, UV coordinates, and epsilon values.

---

# TL;DR

### What I implemented
- Image texturing with nearest, bilinear, and trilinear interpolation
- Bump mapping
- Normal mapping
- Procedural textures (Perlin noise, checkerboard)
- Background textures (`replace_background`, `replace_all`)
- UV coordinate handling with face offsets
- Normalizer field support

### What I fixed
- Bump mapping epsilon too large (0.01 → 0.001)
- Earth texture rotation (`atan2` sign was wrong)
- Grayscale bilinear interpolation sampling only one corner
- Bump gradient calculation (forgot to divide by epsilon)
- Tangent space edge case when normal is parallel to up vector
- Background textures showing black ([0,1] vs [0,255] confusion)
- Face offset parsing with large/negative offsets
- Empty VertexData causing crashes
- UV clamping breaking tiled textures
- Procedural bump mapping not working (UV vs position space)
- Normalizer field not being applied
- `replace_all` textures completely white
- Normalizer incorrectly applied to bump maps

### What I learned
- Color spaces will haunt your dreams
- Epsilon is never what you think it should be
- Every texture bug looks like a different bug



# The Sun Rises from the West

The sphere UV mapping scene has two Earth globes. One should show Asia, the other should show America. Mine showed Asia on both, but the second one was also stretched and glitchy.

The hack? The UV formula for spheres:

```cpp
// Wrong
float u = atan2(z, x) / (2.0 * M_PI) + 0.5;

// Correct
float u = -atan2(z, x) / (2.0 * M_PI) + 0.5;
```

That minus sign. That one minus sign. Hours of debugging for a minus sign.

| Before (Asia twice, second one glitchy) | After (Asia and America) |
| --- | --- |
| ![](./faulty_hw4/sphere_nobump_bump_rotation_problem_and_bump_glitchy.png) | ![](./outputs_hw4/sphere_nobump_bump.png) |

Also had to add proper UV wrapping with `u = u - floor(u)` to handle values outside [0,1].



# The Walls That Vanished

The Killeroo scene has a dinosaur in a Cornell box with bump-mapped walls. My version had... a dinosaur floating in the void.

Turns out `sampleTexture()` returns colors in [0,1] range, but the renderer expects [0,255] for the final output. When I returned the background texture as [0,1], the values got clamped to basically black.

| Before (The void) | After (Actual walls) |
| --- | --- |
| ![](./faulty_hw4/killeroo_bump_walls_no_bump_problem.png) | ![](./outputs_hw4/killeroo_bump_walls.png) |

The fix in `computePixelColor()`:
```cpp
// Scale background texture by 255 before returning
return backgroundTextureColor * 255.0;
```


# The VeachAjar Mystery

The Veach ajar scene is a classic test scene with teapots, a door, and a painting on the wall. My version had:
- A black door (no wood texture)
- A black painting (no image)
- Dark table
- Overall sadness

This was actually multiple bugs stacked on top of each other:

1. **Normalizer field not applied**: Wood textures had `Normalizer: 30` but I was dividing by 255 always
2. **`replace_background` parsing failures**: String/number type mismatches
3. **Background still black**: Even after fixing parsing, the [0,1] → [0,255] issue struck again

| Before (Everything dark and missing) | After (The scene actually renders) |
| --- | --- |
| ![](./faulty_hw4/VeachAjar_dark_and_missing.png) | ![](./outputs_hw4/VeachAjar.png) |

The normalizer fix was interesting. Some textures use values outside [0,255]:
```cpp
// For Normalizer=255: pixel/255 → [0,1]
// For Normalizer=30:  pixel/30 → values up to ~8.5
// For Normalizer=1:   pixel/1 → direct values [0,255]
float normalizerCorrection = 255.0 / normalizer;
```

# The Noisy Checkerboard

The plane_bilinear scene is just a checkered floor extending to the horizon. The ground truth shows clean moiré aliasing patterns. Mine looked like someone threw random pixels at the screen.

The issue: even with `numSamples=1`, I was using random subpixel offsets.

```cpp
// Wrong - random offset for each pixel
float ksi_1 = uniform_random(0, 1);
float ksi_2 = uniform_random(0, 1);

// Correct for single sample - use pixel center
if (numSamples == 1) {
    ksi_1 = 0.5f;
    ksi_2 = 0.5f;
}
```

| Before (Random noise) | After (Clean-ish aliasing) |
| --- | --- |
| ![](./faulty_hw4/plane_bilinear_dark_and_light_problematic.png) | ![](./outputs_hw4/plane_bilinear.png) |

The ground truth has that beautiful moiré pattern because it's deterministically sampling pixel centers. My noisy version was sampling random positions within each pixel, turning predictable aliasing into unpredictable chaos.

# Ellipsoids Gone Stretchy

The ellipsoids scene has various textured ellipsoids. The standing one in the middle should have a nice checkerboard pattern, but mine had stretched vertical stripes.

The issue was UV clamping vs. wrapping. When texture coordinates go beyond [0,1] (like `TexCoordData` with values like "0 4 4 4 4 0" for tiling), clamping collapses everything to the edge texels.

```cpp
// Wrong - clamp to [0,1]
u = std::clamp(u, 0.0f, 1.0f);

// Correct - wrap for tiling
u = u - floor(u);
v = v - floor(v);
```

| Before (Stretched stripes) | After (Proper checkerboard) |
| --- | --- |
| ![](./faulty_hw4/ellipsoids_texture_wrong_stretch.png) | ![](./outputs_hw4/ellipsoids_texture.png) |



# Galactica Static

### Face Offset Parsing

Scenes with `_vertexOffset` and `_textureOffset` in faces needed special handling. The formula `(raw - 1) + (offset - 1)` failed for large/negative offsets. Changed to `raw + offset - 1`.

In the below faulty case I was also having some issues with mapping the background as well. So it looks really awkward on the left.


| Before (Trash) | After |
| --- | --- |
| ![](./faulty_hw4/galactica_static.png) | ![](./outputs_hw4/galactica_static.png) |



## More Small Bug Fixes

### Grayscale Bilinear Bug
For 1-channel images, I was only sampling one corner and using it for all four. Fixed by properly sampling all four corners.

### Tangent Space Edge Case
Cross product with up vector (0,1,0) fails when normal is parallel to it. Added check:
```cpp
if (fabs(dotProduct(normal, up)) > 0.999) {
    tangent = normalize(crossProduct(normal, right));  // Use right vector instead
}
```

### Procedural Bump Mapping
Perlin and checkerboard procedural textures ignore UV coordinates — they only depend on position. So computing gradients by offsetting UVs gave zero gradients. Fixed by computing gradients in position space instead.


# Results Gallery

Some of my favorite renders from HW4:

| Scene | Image |
| --- | --- |
| Wood Box | ![](./outputs_hw4/wood_box_all.png) |
| Perlin Bump Cube | ![](./outputs_hw4/cube_perlin_bump.png) |
| Brick Wall Normal Map | ![](./outputs_hw4/brickwall_with_normalmap.png) |
| Galactica Dynamic | ![](./outputs_hw4/galactica_dynamic.png) |
| Dragon | ![](./outputs_hw4/dragon_new.png) |
| Tap (from HW3, now with textures) | ![](./outputs_hw4/mytap_final.png) |



# Benchmark Results

All results with BVH enabled and multi-threading.

| Scene | Pre (ms) | Render (ms) | Total (ms) | BVH | MT | Final Image |
| --- | --- | --- | --- | --- | --- | --- |
| galactica_dynamic | 0 | 56,419 | 56,419 | ✅ | ✅ | ![](./outputs_hw4/galactica_dynamic.png) |
| mytap_final | 13 | 21,363 | 21,376 | ✅ | ✅ | ![](./outputs_hw4/mytap_final.png) |
| sphere_normal | 0 | 8,128 | 8,128 | ✅ | ✅ | ![](./outputs_hw4/sphere_normal.png) |
| VeachAjar | 41 | 4,116 | 4,157 | ✅ | ✅ | ![](./outputs_hw4/VeachAjar.png) |
| killeroo_bump_walls | 25 | 3,844 | 3,869 | ✅ | ✅ | ![](./outputs_hw4/killeroo_bump_walls.png) |
| galactica_static | 0 | 631 | 631 | ✅ | ✅ | ![](./outputs_hw4/galactica_static.png) |
| sphere_perlin_bump | 0 | 161 | 161 | ✅ | ✅ | ![](./outputs_hw4/sphere_perlin_bump.png) |
| bump_mapping_transformed | 0 | 154 | 154 | ✅ | ✅ | ![](./outputs_hw4/bump_mapping_transformed.png) |
| sphere_perlin_scale | 0 | 150 | 150 | ✅ | ✅ | ![](./outputs_hw4/sphere_perlin_scale.png) |
| sphere_perlin | 0 | 147 | 147 | ✅ | ✅ | ![](./outputs_hw4/sphere_perlin.png) |
| wood_box_all | 0 | 134 | 134 | ✅ | ✅ | ![](./outputs_hw4/wood_box_all.png) |
| brickwall_with_normalmap | 0 | 115 | 115 | ✅ | ✅ | ![](./outputs_hw4/brickwall_with_normalmap.png) |
| cube_wall_normal | 0 | 112 | 112 | ✅ | ✅ | ![](./outputs_hw4/cube_wall_normal.png) |
| cube_perlin_bump | 0 | 110 | 110 | ✅ | ✅ | ![](./outputs_hw4/cube_perlin_bump.png) |
| wood_box | 0 | 109 | 109 | ✅ | ✅ | ![](./outputs_hw4/wood_box.png) |
| sphere_nobump_bump | 0 | 107 | 107 | ✅ | ✅ | ![](./outputs_hw4/sphere_nobump_bump.png) |
| wood_box_no_specular | 0 | 107 | 107 | ✅ | ✅ | ![](./outputs_hw4/wood_box_no_specular.png) |
| ellipsoids_texture | 0 | 103 | 103 | ✅ | ✅ | ![](./outputs_hw4/ellipsoids_texture.png) |
| cube_perlin | 0 | 102 | 102 | ✅ | ✅ | ![](./outputs_hw4/cube_perlin.png) |
| sphere_nobump_justbump | 0 | 101 | 101 | ✅ | ✅ | ![](./outputs_hw4/sphere_nobump_justbump.png) |
| cube_wall | 0 | 100 | 100 | ✅ | ✅ | ![](./outputs_hw4/cube_wall.png) |
| sphere_nearest_trilinear | 0 | 100 | 100 | ✅ | ✅ | ![](./outputs_hw4/sphere_nearest_trilinear.png) |
| cube_waves | 0 | 99 | 99 | ✅ | ✅ | ![](./outputs_hw4/cube_waves.png) |
| cube_cushion | 0 | 99 | 99 | ✅ | ✅ | ![](./outputs_hw4/cube_cushion.png) |
| sphere_nearest_bilinear | 0 | 98 | 98 | ✅ | ✅ | ![](./outputs_hw4/sphere_nearest_bilinear.png) |
| plane_trilinear | 0 | 59 | 59 | ✅ | ✅ | ![](./outputs_hw4/plane_trilinear.png) |
| plane_bilinear | 0 | 56 | 56 | ✅ | ✅ | ![](./outputs_hw4/plane_bilinear.png) |
| plane_nearest | 0 | 56 | 56 | ✅ | ✅ | ![](./outputs_hw4/plane_nearest.png) |


# What I Learned

HW4 taught me that:

- Texturing has more edge cases than actual rendering code
- The difference between [0,1] and [0,255] will ruin your week
- Epsilon values are never obvious — always too big or too small
- A single minus sign can rotate your entire planet

*But honestly, after four homeworks, the raytracer feels complete. It can handle meshes, transformations, materials, lights, sampling, motion blur, depth of field, area lights, and now textures. Looking back at the creepy bunny from HW1, we've come a long way.*

Until the next homework.

