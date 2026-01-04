# Logs

## Jan 2 Friday Morning 

I am starting this homework a little bit of late because of the holiday. It looks like I have one final battle with this raytracer. I will kick start with the parsing stuff to get up and running quickly.

The homework is about HDR rendering and advanced lights. Three new light types: directional, spot, and environment lights. Plus tone mapping operators and HDR file support. This is the last one, let's make it count.

Started with the easy stuff - adding the structs for the new lights. DirectionalLight is straightforward - just direction and radiance, no distance attenuation. SpotLight needs position, direction, intensity, and two angles (coverage and falloff). SphericalDirectionalLight is the interesting one - it uses an HDR environment map image. Added them to scene.h following the same pattern as PointLight and AreaLight.

Parsing was also straightforward. Just followed the existing pattern for PointLight/AreaLight parsing. The JSON structure is pretty clear from the spec. One thing I noticed - SphericalDirectionalLight has optional `_type` (latlong or probe) and `Sampler` (uniform or cosine), so I added defaults for those.

## Jan 2 Friday Afternoon

Now for the tone mapping stuff. Added `TonemapSettings` struct to Camera. The spec says it can be a list, so I made it a vector. Each entry has TMO (Photographic/Filmic/ACES), TMOOptions (key and burn-out percent as a string), Saturation, Gamma, and Extension.

Parsing the tone mapping was a bit tricky because TMOOptions is a string like "0.18 1" that I need to parse into two doubles. But it's working now.

Next up: HDR image loading. I already have `tinyexr.h` included, so that's good. Modified the Image struct to have both `data` (for LDR) and `hdrData` (for HDR), plus an `isHDR` flag. The loading code now checks the file extension - if it's `.exr`, use `LoadEXR`, if it's `.hdr`, use `stbi_loadf`, otherwise use regular `stbi_load`. The destructor needed updating to handle both cases - EXR uses `free()` while HDR uses `stbi_image_free()`.

## Jan 2 Friday Evening

Created the tone mapping module (`tonemap.h`/`tonemap.cpp`). The photographic operator follows the Reinhard paper - Equations 1-4. The burn-out percentile calculation requires sorting all luminances first, which is a bit expensive but necessary. Filmic and ACES are similar curve-based operators. The saturation application uses the formula from the spec: `Ro = Yo * (R/Yi)^s`.

Also created `hdr_io.h`/`hdr_io.cpp` for writing HDR files. `writeEXR` uses `SaveEXR` from tinyexr, and `writeHDR` uses `stbi_write_hdr`. EXR expects RGBA so I convert RGB to RGBA by adding alpha=1.0.

Now the big change - modifying the rendering pipeline to use float buffers instead of clamping to [0,255]. Changed `unsigned char* image` to `float* hdrImage` everywhere. The thread functions now accumulate float colors directly. After rendering, I check if the output filename is HDR (.exr or .hdr), and if so, save it directly. If tone mapping is specified, apply it and save PNG. Otherwise, convert to LDR with clamping for backward compatibility.

The `__compute` function was clamping colors to [0,255] which I removed - we want HDR values now!

## Jan 3 Saturday Morning

Implemented the three new light types in `computeShading()`. Directional lights are simple - no distance attenuation, just check if surface faces the light, shadow test, then add diffuse + specular. Pretty straightforward.

Spot lights are more interesting. Need to calculate the angle between the light direction and the vector to the intersection point. If it's outside the coverage angle, skip. If it's between falloff and coverage, apply the spot attenuation formula from the spec:

```cpp
double spotAttenuation = std::pow((cosAlpha - cosCoverage) / (cosFalloff - cosCoverage), 4.0);
```

Then multiply by distance attenuation like point lights. The formula in the spec uses degrees, so I convert to radians.

Environment lights are the most complex. Need to sample a direction from the upper hemisphere - either uniform or cosine-weighted. For cosine sampling, I use inversion sampling: `cosTheta = sqrt(r1)`, `sinTheta = sqrt(1 - r1)`, `phi = 2*PI*r2`. Then transform to world space using the surface normal and an orthonormal basis.

Convert the direction to UV coordinates - latlong uses `atan2` and `acos`, probe uses a different formula. Then sample from the HDR environment map using bilinear interpolation. The tricky part is that HDR images use `float*` data, not `unsigned char*`, so I had to modify the bilinear interpolation to handle both cases. Actually, I just wrote the lookup inline in the environment light code since it's HDR-specific.

The PDF calculation is important - for uniform it's `1/(2*PI)`, for cosine it's `cos(theta)/PI`. Then divide the radiance by PDF to get the correct estimate.

## Jan 3 Saturday Afternoon

Added light transformations. Directional lights transform the direction vector, spot lights transform both position and direction. SphericalDirectionalLight doesn't need transformations since it uses image lookup.

Updated the Makefile - the new source files (tonemap.cpp, hdr_io.cpp) are automatically included with `*.cpp`, so that's fine. Just added `-I.` to be safe.

Found a bug - in the single-threaded rendering path, I was still using `image` instead of `hdrImage` in a couple places. Fixed that. Also the iterative canvas writing needed updating to handle float buffers.

One thing I'm worried about - the environment light sampling. The cosine-weighted hemisphere sampling should be correct, but I want to test it with actual scenes to make sure the illumination looks right. The UV conversion formulas from the spec should be correct though.

## Jan 3 Saturday Evening

Everything compiles! No linter errors. The implementation is complete. I've added:

- Three new light types (directional, spot, environment)
- HDR image loading (EXR and HDR formats)
- HDR image writing
- Three tone mapping operators (Photographic, Filmic, ACES)
- Modified rendering pipeline to use float buffers

Now I need to test it with actual scenes. The environment maps should create some really nice lighting effects. And the tone mapping should make HDR renders look good when converted to PNG. This is the final homework - let's see if everything works as expected!

## Jan 3 Saturday Night - Bug Fixes

Ran into several issues during testing that needed fixing:

### 1. Tone Mapping - Log-Average Luminance

The photographic tone mapping was way too bright! The issue was I was using `key * Yi` directly instead of properly normalizing by the scene's log-average luminance. According to Reinhard's paper, you need to:

1. Compute log-average luminance: `Lw = exp(sum(log(δ + Li)) / N)`
2. Scale each pixel by `(key / Lw) * Yi`, not just `key * Yi`

This makes the average luminance of the scene map to the key value (typically 0.18 = middle gray), which is essential for proper HDR→LDR conversion.

### 2. Environment Light Monte Carlo Estimator

The environment light formula had wrong scaling. For cosine-weighted sampling with physically correct Lambertian BRDF (ρ/π):
- PDF = cos(θ)/π
- Integrand = L * (ρ/π) * cos(θ)
- Estimator = Integrand/PDF = L * ρ (the π terms cancel!)

So the correct formula is just `envContribution = radiance * diffuseReflectance` - no π multiplier needed.

### 3. Latlong Environment Map UV Mapping

This one took several iterations to get right! The environment maps were showing the wrong parts of the scene - horizontally flipped and rotated. The correct formula for latlong mapping is:

```cpp
u = 0.5 + std::atan2(direction.x, -direction.z) / (2.0 * M_PI);
v = std::acos(clampedY) / M_PI;
```

Where:
- `atan2(x, -z)` computes azimuth angle (horizontal rotation around Y-axis)
- `acos(y)` computes polar angle from top pole
- The `0.5` offset centers forward direction (-Z) at panorama center

After these fixes, both mirror_sphere_env and glass_sphere_env scenes match the ground truth perfectly!

## Jan 4 Sunday Evening - Bug Fixes

Found a critical bug when testing the `sphere_point_hdr_texture` scene - it was rendering completely black! Turned out the texture sampling functions weren't handling HDR images at all.

### 1. HDR Texture Sampling

The problem was in `sampleImageTexture`, `bilinearInterpolation`, and `trilinearInterpolation` - they all had this check:

```cpp
if (!image || !image->data) {
    return VectorFloatTriplet{0.0, 0.0, 0.0};
}
```

For HDR images, the data is stored in `image->hdrData`, not `image->data`, so the check was failing and returning black. Also, the sampling code was trying to access `image->data[idx] / 255.0` which doesn't make sense for HDR - those values are already floats in proper radiance units.

Fixed by checking for both LDR and HDR data:

```cpp
if (!image || (!image->data && !image->hdrData)) {
    return VectorFloatTriplet{0.0, 0.0, 0.0};
}
```

Then branching based on `image->isHDR` - for HDR, use `hdrData` directly (no division by 255), for LDR, use `data` with normalization. Updated all three interpolation functions the same way.

### 2. Normalizer for HDR Textures

There was another subtle issue - the normalizer correction was applying `255.0 / normalizer` to all textures. This makes sense for LDR (scaling from [0,1] to [0,255/normalizer]), but HDR values are already in proper radiance units and shouldn't get multiplied by 255.

Fixed by checking if the image is HDR before applying the normalizer. For HDR textures, only apply the normalizer division if it's not 1.0 (to allow user scaling), but skip the 255x multiplication. For LDR, keep the existing behavior.

Now `sphere_point_hdr_texture` renders correctly with the grace cathedral environment map showing on the sphere instead of being black!
