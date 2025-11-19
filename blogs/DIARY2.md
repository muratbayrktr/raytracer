# Raytracer with Transformations and BVH


# Logs

## 2025-11-10

It's been a while since the assignment has been announced. I started this one a little late and unfortunately I have 2 bugs to fix before starting the actual homework. I have a week or so and I plan to finish these two main bugs in couple of days so at least after friday to sunday I will have full dedicated time just for transformations.

## 2025-11-11

I finally got around to fixing and improving the BVH implementation that I mentioned was buggy at the end of the first homework. The issues were pretty interesting - some scenes were working fine (like bunny) but others were completely broken (cornellbox meshes disappearing, chinese dragon looking like scattered points). Classic case of "works on my machine" but actually "works on some scenes" lol. 

Funny how, I always thought BVH was the problem but it turns out the actual issue is at the intersection of multiple edge cases. 
For instance, chinese dragon had too small faces and I was using the default epsilon of 1e-6. :D I haven't noticed this until I debugged the calculations with pen and paper. On the other hand, I rewrote the bvh implementation similar to the one before but I noticed the trees were somehow unbalanced. Fixed it by playing around the splitAxis. 

Also one performance improvement I made was replacing the `std::vector<int>` stack with a fixed-size array `int stack[64]`. This eliminates dynamic allocations during traversal and should be faster. The stack depth rarely exceeds 64 levels for typical meshes, so this is a safe optimization. Plus it's more cache-friendly.

Now my bvh is waay faster.

## 2025-11-12

Dielectrics. I had problems with them in the previous homework. I did it on the last day of the homework 1 and couldn't get it running properly. I noticed culling was messing up with my calculations in refraction. Since we move inside the object, the second hit always from the inside and that's why they were off. I noticed this after trying to render the mirror scene in the second homework. Because mirrors were kinda dark at some angles.

Also ambient was not added on dielectric cases, it's assignment was after the return statements of dielectric and conductor's.

## 2025-11-14

I started implementing transformations. Started with scene parsing. Handled most of it. I think I am good to use them. 

Simple scene is working fine however, mesh instances are not working properly. Also I noticed that when there is a mesh instance, the code runs extremely slow. I think something is off I added performance profiling and the time spent in tracing and the number of calls. Will continue tomorrow.

```
=== Performance Stats ===
Intersect calls: 695184
BVH traversals: 65536
Triangle tests: 775184
Transformed mesh calls: 695184
World bounds rejects: 629648
World bounds accepts: 65536
Ray transforms: 65536
World bounds reject rate: 90.5727%

=== Timing Breakdown (ms) ===
Time in intersect(): 99.074 ms
Time in rayHitsMesh(): 462.72 ms
  - BVH traverse: 34.1818 ms
  - Ray transform: 12.4943 ms
  - AABB tests: 18.4531 ms
  - Triangle tests: 0 ms
  - Back transform: 1.77951 ms
=========================
```

## 2025-11-15

After performance profiling I found out below issues:

For the simple scene:
- Triangle tests: 10.6ms (775,184 tests)
- BVH traverse: 36.5ms
- Transforms: 11.7ms
- Total tracked: ~84ms out of 506ms in rayHitsMesh
For the large scene:
- 172,364 triangle tests for just 10,000 pixels
- ETA: ~72 minutes (4334 seconds)
- Only 27% world bounds rejection rate (vs 90% in simple scene!)

This means the rayHitsMesh timing only accounts for 131ms, but intersect() takes 46 seconds. That means something ELSE in the intersect() function is eating all the time.

```
Time in intersect(): 102594 ms       (102 SECONDS)
  - Planes: 0.408382 ms              (negligible)
  - Triangles: 0.568624 ms           (negligible)
  - Meshes: 123.744 ms               (0.12 seconds)
  - Spheres: 0.387979 ms             (negligible)
  - Instances: 102465 ms             (102 SECONDS)
  - PostProcess: 0.4962 ms           (negligible)
```
and of course this is because I was COPYING a vector of 1.8 million faces for every ray that passes the world-space bounds test xd. Moving forward... lessons learned... 

I noticed a bug in lookAt camera type, I will put an example for it in the blog it's a funny one but it seems because of it I couldn't properly render dragon metal.

```cpp
// I forgot the conversion xd
float fovY = parseSingleValue<float>(cameraData["FovY"]) * M_PI / 180.0f; // convert degrees to radians
```

It's been a while since transformations and instancing are (not) working properly :-). But I found out issues regarding bvh side and precomputed determinants. I was just using the precomputed determinants for transformed meshes and that was causing issues.

What's more is I had to first handle BVH and then precompute transformation because precomputation needs it for world-space bounds. 

I was able to get the marching dragons run somehow but the color was off. This was mainly because MeshInstances in JSON don't have a Material field, but the parser tried to read it, causing a crash. Instead it should have used the default value.

I thing this is the last bug that I will have solved today. Two berserkers was pain in the... neck. The second berserker has transformation s2 = "0.75 -0.75 0.75" (negative Y scale), which is a reflection. This flips the triangle winding order, causing normals to point inward; thus, everything renders dark. Now it's working... ish. Because some parts of the berserker guys are still dark and I will try to fix it at last.

## 2025-11-16

Today I wrote the `benchmark.py` and the `runner.py` scripts. My cpp code generates a metadata.json file for every camera in every scene. Like this:
```json
{
    "enableBackFaceCulling": false,
    "isMultiThreaded": true,
    "preprocessingTimeMs": 615.0,
    "renderTimeMs": 2208.0,
    "sceneName": "dragon_metal.png",
    "totalTimeMs": 2823.0,
    "useBVH": true
}
```
This pave the way for creating the benchmark for my blog. And the runner code helped me convert frames into a video to create the same outputs just like the ground truth.

PS: I tried to fix berserkers but I couldn't get it working properly. I still don't understand why some parts are still dark.
