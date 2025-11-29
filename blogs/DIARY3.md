# Logs 

## 2025-11-25 Evening

-Insert GTA San Andreas Beginning scene meme here- "Ah... Here we go again.". Since I am getting used to it, this time I won't wander around and jump directly into writing the parser verbatim. I remember the discussions well from the lectures; though, I might need to revisit. Yet the idea is simple we got some configs for our sampling methods and we will deliver multi-sampling via those parameters and probably the pseudo codes from the lecture slides. I don't expect this one to be tough. 

Here is what I have so far:
- Suggested method is `jittered sampling`
- `NumSamples`: Inside the `Camera` object. Hw states: "This element is defined as a perfect square (e.g., 1, 2, 4, 16, etc.)"
```
i = 0
for y = 0 to nRows – 1
    for x = 0 to nCols – 1
        𝜉_1=uniform_random(0, 1)
        𝜉_2=uniform_random(0, 1)
        samples[i].x = (x + 𝜉_1) / nCols
        samples[i].y = (y + 𝜉_2) / nRows
        i = i + 1 
```
- `ApertureSize`: Inside the `Camera` object. This may not exists, so I might use a flag for this one to conditionally enable depth-of-field thing. Lens is assumed flat square, this is easier. I won't have to deal with square root etc. to adjust a good sampling distribution on the spherical one. Also one notable thing: "Note that the aperture’s center position is the camera origin. The aperture sits on the w = 0 plane."
- `FocusDistance`: Also inside the `Camera`.
- `AreaLight`: Defined as square
    - Position: The center point of the area light
    - Normal: The surface normal of the light
    - Size: The edge length of the area light. Its area will be equal to the square of this value
    - Radiance: Radiance of the light

    There is some workload here. I gotta be careful with calculations.

- `MotionBlur`: This one simulates the blurry look of fast moving objects. It's defined for objects using the "MotionBlur" element in the JSON (like in dragon dynamic.json). For simplicity, we're only dealing with translational movements - no rotations.

MotionBlur simulates the blur caused by objects moving during the camera’s exposure. In the JSON scene description (e.g., dragon_dynamic.json), objects can have a MotionBlur element that specifies a translation vector describing how much the object moves over the entire shutter interval.

To implement this, each ray carries a time parameter t ∈ [0, 1], which is a normalized shutter time:
- t = 0: shutter opens → object is at its base transformed position.
- t = 1: shutter closes → object is translated by the full MotionBlur vector.

For a ray with time t, the object’s transform becomes:
```
T(t) = T_base * translate(t * motionBlurVector)
```
Each primary sample (camera ray) draws a random t uniformly in [0, 1], and all secondary rays spawned from that primary ray reuse the same t. This ensures temporal coherence for a given pixel path and produces a smooth stochastic motion blur, as described in distribution_ray_tracing_explained.pdf.

- `Roughness`: This is a scalar parameter in the Material that controls how glossy or diffuse specular interactions appear for mirrors, conductors, and dielectrics. It determines how far the reflected (and, for dielectrics, refracted) direction can deviate from the ideal perfect direction.
- Small roughness : directions are tightly clustered around the perfect reflection/refraction => sharp, polished highlights.
- Large roughness : directions are widely spread => blurry reflections and softer highlights.

Implementation-wise, I first compute the ideal reflection or refraction direction, then sample a random direction in a lobe around it whose spread is controlled by the roughness value. For dielectrics, the ideal refracted direction is computed via Snell’s law, and the outgoing direction is sampled around that. This follows the stochastic sampling approach in distribution_ray_tracing_explained.pdf.


I handled the parsing and also the initial sampling. It didn't work out at first but then I noticed I was only passing the jitter as if they were the subpixels.

```cpp
for (int y = args->startY; ...) {
    for (int x = 0; x ...) {
        ...
        for (int k = 0; k < camera->numSamples; k++) {
            VectorFloatTriplet jitter = camera->samples[sampleIndex + k];
            double sx = x + jitter.x;   // subpixel x 
            double sy = y + jitter.y;   // subpixel y
            pixelColor = pixelColor + __compute(..., sx, sy, ...); // previously I was using jitter.x and jitter.y which created nonsense xd
        }
        pixelColor = pixelColor * (1.0 / camera->numSamples);
        ...rest
        
    }
}
```

## 2025-11-26

Keeping this as a diary was really good idea because now I am able to instantly start working on things as I am continuing from where I left off by recalling from previous days logs. Now that I am done with sampling part. I am moving on to aperture and focus distance.

It just took 6-7 lines of code and I wasn't expecting this and I think it works just fine. 

After I implement the sampling the speed kinda dropped but I will take care of it later not a big problem for now.

I now tried some scenes with brushed metal and area light. I thought area light is working okay-ish but it doesn't create a brightness on the ceiling. Similarly conductor and dielectrics are kinda off when considering the area light. Something is not adding up but I'll find it.


## 2025-11-28

Today I finally tracked down the infamous "green glass" bug. Transformed dielectric objects were rendering with this sickly green tint instead of proper glass refraction. Classic.

The culprit? A sneaky normal-flipping optimization that was "helping" everywhere except where it mattered most. In two places in `utils.cpp`, I had code like:

```cpp
// Ensure normals face the camera (critical for mirror reflections and lighting)
if (dotProduct(worldGeomNormal, ray.direction) > 0.0) {
    worldGeomNormal = -worldGeomNormal;
    worldShadingNormal = -worldShadingNormal;
}
```

Looks innocent, right? Even has a helpful comment. But here's the thing - dielectrics use this check to determine if we're entering or exiting the glass:

```cpp
bool entering = dotProduct(ray.direction, N) < 0.0;
```

When you always flip the normal to face the camera, `dotProduct(ray.direction, N)` is ALWAYS negative. So `entering` is ALWAYS true. The refracted ray traveling inside the glass thinks it's still entering, uses the wrong refractive indices, computes garbage directions, and boom - green artifacts everywhere.

The fix was simple once I understood the problem:
1. Removed the automatic normal flipping in both `rayHitsMesh` (for transformed objects) and the post-processing section (for smooth shaded meshes)
2. Added the normal flip specifically in `computeShading()` but ONLY for non-dielectric materials

Also found a bonus bug while I was at it - `rayHitsSphere` was using `min(t1, t2)` for sphere intersection which completely breaks when the ray starts inside the sphere (like, you know, when you're inside a glass sphere after refraction). Fixed that too by properly selecting the closest positive t value.

Moral of the story: what helps opaque materials can absolutely destroy transparent ones. Dielectrics are special snowflakes and they need their original normals to know which side they're on.

#### Mesh Scaling Bug

Today I tracked down why scaled meshes were “disappearing” from the render while still casting shadows.
Root cause: for transformed meshes I was using the global world-space t_min as the intersection cutoff in object space. When the mesh was scaled (especially down), the object-space hit distance became larger than world t_min, so all primary ray hits were rejected.

Fix: for transformed meshes I introduced a separate local t_min in object space (initialized to inf) and kept the incoming t_min only as original_t_min in world space. After intersecting in object space and back-transforming the hit point, I compare the resulting worldDistance with original_t_min and update the global t_min only if this hit is closer. This decouples object-space intersection distances from world-space pruning and prevents scaled meshes from disappearing.