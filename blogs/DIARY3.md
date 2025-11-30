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
- t = 0: shutter opens => object is at its base transformed position.
- t = 1: shutter closes => object is translated by the full MotionBlur vector.

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

## 2025-11-29

Today I implemented motion blur. The parsing was already done - each object has `motionBlur` and `hasMotionBlur` fields, and the `Ray` struct already had a `time` field that was just sitting there unused.

The idea is straightforward: motion blur is a world-space translation that happens AFTER all other transformations. At time t=0, the object is at its transformed position. At t=1, it's translated by the full motionBlur vector. Each sample gets a random time in [0,1], and all secondary rays (reflections, refractions, shadows) keep the same time as their parent ray.

Implementation approach - instead of actually moving objects (expensive), I offset the ray in the opposite direction:
1. Object at time t is at position P + t*motionBlur
2. Equivalently, test ray with origin shifted by -t*motionBlur against object at P
3. If hit, shift hit point back by +t*motionBlur

The changes were pretty localized:
- `castRay` now accepts a time parameter and stores it in the ray
- `reflect` and `refract` propagate time to secondary rays
- Shadow rays in `isInShadow` and `computeShading` also propagate time
- Added two helper functions: `applyMotionBlurToRay` and `correctHitPointForMotionBlur`
- In `intersect`, for each object type (planes, triangles, meshes, spheres, mesh instances) I create an offset ray if the object has motion blur, do the normal intersection test, then correct the hit point if there was a hit
- In `raytracer.cpp`, each sample generates `sampleTime = uniform_random(0.0, 1.0)` and passes it to `__compute`

One subtlety: since motion blur is just a translation, the distance along the ray is unchanged. So t_min comparisons work correctly even when different objects have different motion blur offsets.

The AABB early rejection for mesh instances uses the offset ray but against the static bounding box - this is conservative (might miss some early rejections) but correct. Could optimize later by expanding the AABB by the motion blur extent.

On top of that, I finally fixed an annoying near-plane / image-plane bug. The all-intersection path and BVH traversal were happily reporting hits that were technically between the camera origin and the image plane when `camera.nearDistance` was small, which meant some objects could “pop” in front of the film even though they should have been clipped. The first attempt at fixing this used `camera.nearDistance` directly as a radial cutoff (distance from the camera), which behaved like clipping with a sphere and produced that weird “growing circle” effect when I changed nearDistance.

The real fix was to compute, for each primary ray, where that specific ray actually crosses the image plane and use that as the per-ray minimum distance. In `intersect()` I now do:

```cpp
double minDistance = 0.0;
bool isPrimaryRay = (ray.depth == 0 && !ray.shadowRay && !ray.reflectionRay && !ray.refractionRay);
if (isPrimaryRay && !scene.cameras.empty()) {
    const Camera& cam = scene.cameras[scene.currentCameraIndex];
    VectorFloatTriplet n = normalize(cam.gaze);
    VectorFloatTriplet e = cam.position;
    VectorFloatTriplet o = ray.origin;
    VectorFloatTriplet d = normalize(ray.direction);

    double denom = dotProduct(d, n);
    if (std::fabs(denom) > 1e-9) {
        double numer = cam.nearDistance - dotProduct(o - e, n);
        double t_plane = numer / denom;
        if (t_plane > 0.0) {
            minDistance = t_plane;
        }
    }
}
```

Then I thread this `minDistance` into `rayHitsPlane`, `rayHitsSphere`, `rayHitsTriangle`, `rayHitsMesh`, and `MeshBVH::traverse`. For untransformed primitives/meshes, I reject hits with `t < minDistance` (so anything before the image plane along that ray is ignored). For transformed ones (where local-space t doesn’t match world-space distance), I let the code compute the world-space hit point first and then discard any candidate whose world-space distance along the ray is less than `minDistance` before updating `t_min`. Secondary rays (shadows, reflections, refractions) still see everything from their own origin since they don’t use the camera’s image plane for clipping.

On top of that I added an “actual raytracer-style” progressive preview for multi-sampling. Instead of doing all samples per pixel in one go, I can now iterate **level-wise** over the sample index: for each sample `k`, I sweep all pixels, accumulate into a floating-point `accum` buffer, and periodically (every few samples) normalize and dump the current state. The precomputed `VectorFloatPenta` per pixel/sample (jitter, time, and extra random dims) made this trivial to wire up – no RNG cost in the hot path. I also wrapped a small SDL-based GUI around the iterative path so that when it’s enabled, the renderer opens a window at the image resolution and re-blits the current `accum` buffer each time I write the `_iterative.png` snapshot. Net effect: I can literally watch the image converge sample-by-sample like a “real” progressive raytracer, while still keeping the old one-shot sampling mode as the default.

## 2025-11-30

Today I finally understood that my “tap + water inside the tap” disaster was not some mysterious motion blur glitch, not a smoothing bug, not even a transform bug – it was pure, old-fashioned intersection ordering chaos in concave, layered geometry. Once I wrote down the sequence of hits for a primary ray (tap outer => water outer => water inner => tap inner => exit…) it clicked: different parts of the code were comparing different notions of “distance” (local-space t here, world-space t there, `original_t_min` vs `local_t_min`), and sometimes an inner surface was winning over an outer one just because it was competing in the wrong metric. The symptoms matched perfectly: water somehow “in front” of the tap, ghosty transparency, and objects disappearing or popping when I introduced transforms or BVH.

The fix was to make the intersection pipeline brutally consistent: **one ray, one distance metric, everywhere**. Every primitive or BVH traversal now produces a candidate hit point, I transform that to world space if needed, compute `distance_world = length(hitPoint_world - originalRay.origin)`, and only then do I let it compete against the global `t_min`. Local/object-space t is now just an internal parameter to find the point, never the thing we compare for “closest hit” once transforms enter the picture. For transformed triangles and meshes (including instances), I use a separate `local_t_min` and `localIntersection` in object space, back-transform the winner, compute its world-space distance, apply near-plane clipping in world coordinates, and only if `distance_world < t_min` do I overwrite the global intersection and `t_min`. No more resetting t_min to “original” after testing candidates, no more mixing local and world distances.

In practice this means outermost layers always win unless there is a real geometric hole: the tap’s outer shell wins over the water, the water wins over its own inner surface, and the inner tap walls never magically float in front. Motion blur stays happy because it just offsets rays and then recomputes the distance from the original origin after correcting the hit point. It’s surprisingly satisfying to see how a single rule – “compare only world-space distance from the original ray” – cleans up a whole zoo of concave mesh bugs and makes the renderer feel a lot more like what real ray tracers actually do under the hood.