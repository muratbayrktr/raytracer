# Raytracer

To get the project up and running run

```
make all
```

after that use following to produce the output images from JSON scene files:

```
./raytracer ../inputs/scene.json
```

where `scene.json` is the JSON scene file you want to render.

# Logs 

## Until 2025-10-21

I decided to include my journey of implementing the raytracer in this README.md file to edit it later. I also plan to keep a casual language to write this so I can be more relaxed and write this as if I'm talking to a friend.

Up until this point I was able to parse the scene file and print the scene summary. I was also able to write an image by assigning a constant color to each pixel, (Woah! 2 for loops for writing an image? I'm a genius!).

I tried to keep the `scene.h` as simple as possible; however, I noticed I misunderstood the faces and vertex data. Vertex data is a list of vertices, and faces is a list of indices to the vertices. I was trying to parse the faces as a list of vertices, which is wrong. 

Along the way I also noticed that the `scene.cpp` file was getting too big, so I decided to split it into smaller files. I created a `utils.h` file to contain the utility functions that I used in the `scene.cpp` file. I also created a `raytracer.cpp` file to contain the main function of the raytracer.

I'm planning to get this up and running for only `bunny.json` scene file for now and then implement other details later. I haven't implemented the other camera type as well. So there's still a lot to do.
My next TODO follows as:
```
TODO: I have to do the following to compute pixelColor:
1. Cast a ray from the camera to the pixel
2. Check if the ray intersects with any object, I'll prolly need to write separate function for different objects
3. If it intersects, compute the color of the pixel
Somewhere around here I gotta apply shading
4. If it doesn't intersect, return the background color

```
That's all for now.

## 2025-10-23 

I am currently working on the `castRay` function. My plan is to implement this function from 1st slide 15th page of the lecture notes **verbatim**. Let's see how it goes.

## 2025-10-24

Alright, time to actually implement that `castRay` function. Followed the lecture notes pretty closely - slide 15 from the first deck. Added the Ray and Intersection structs to `scene.h` to keep track of rays and where they hit stuff.

Also I am glad from the beginning I am using operator overloading for vector operations, that saved me a lot of lines. Nothing fancy, just addition, subtraction, multiplication, cross product, dot product, normalization - the usual. These are gonna be super useful for all the ray-object intersection calculations.

Got the basic camera ray casting working! For each pixel, I'm computing the ray direction using the camera's position, gaze direction, and the near plane.

## 2025-10-25

Today was all about fixing stupid bugs and refactoring. 

First issue: I realized I was mixing up 0-based and 1-based indexing for the vertex data. The JSON files use 1-based indexing, but C++ arrays are 0-based. Had to add `-1` everywhere when accessing vertices. Fixed this for faces, triangles, spheres, and planes.

Also did a major refactor of the mesh normal computation and the PPM writing function. The PPM writer was kinda messy before, so I cleaned it up.

Then I improved the memory management in the main raytracer loop and enhanced the pixel color computation. 

These all happened because I was first testing ray-plane intersection and when I output the image I got the plane on the wrong half of the image. So while trying couple of fixes I did some refactoring.


## 2025-10-26

Alright, so I've been making some solid progress today. Got the basic ray-sphere and ray-plane intersections working. The math from the lecture slides was actually pretty straightforward to implement once I got my head around it.

So here's a fun one, spent way too long debugging why my spheres were looking weird. Turns out I had a sign error in the discriminant calculation for ray-sphere intersection. The formula should have `- r²` but I accidentally wrote `+ r²`.

The correct formula is:
```
D = (d·(o-c))² - (d·d)·((o-c)·(o-c) - r²)
```

Not:
```
D = (d·(o-c))² - (d·d)·((o-c)·(o-c) + r²)
```

lol.

Okay so this was the big one. I implemented the `computeShading` function with ambient, diffuse, and specular components. But when I rendered the scene, I was only getting the ambient color. Everything looked super dark and flat, no matter what I did.

My diffuse and specular calculations were basically being skipped. Facepalm moment for sure.

Fixed it by removing that check entirely (it was redundant anyway since I was already using `max(0.0f, ...)` for the dot products).

I also cleaned up the sphere intersection code - simplified it to just take `min(t1, t2)` instead of checking both intersections separately. Because of this I wasted couple of hours debugging why my spheres are creepy.

Enabled triangle intersection in the intersect function too, though I still need to actually implement `rayHitsTriangle`. 

```
TODO:
1. Implement ray-triangle intersection (I'll use barycentric coordinates)
2. Implement ray-mesh intersection (basically just loop through triangle faces)
3. Shadow rays (check if point is in shadow before adding diffuse/specular)
4. Maybe add some reflection/refraction
5. Test with more complex scenes
```

The spheres with plane scene is rendering correctly now with proper lighting and specular highlights. Pretty happy with it now. 

## 2025-10-26 Evening & Night

I am more of a zone-out person. So although I have made little progress last week, the biggest progress was today which I sat more than 8-9 hours straight. For the last couple of hours I've been working on implementing the ray-triangle intersection.

There were 2 main methods which I could have used: 
1- Lengthy method
2- Barycentric coordinates method

I thought using barycentric coordinates will allow me to early quit with some checks. Moreover, while calculating the determinant for `t` I noticed we were only using the origin and the triangle vertices. **Which means** I don't need to know the ray direction so origin will be constant and I will be able to calculate the determinant pre-tracing.

But now that I am writing what I have done, I noticed I haven't implemented the reflection and refraction yet. So while reflecting the ray I will be sending will be a new ray with a new origin and direction meaning the moment I need to calculate the reflected ray triangle intersection if I use the precomputed determinant I will be screwed. I am glad I am writing this down right now otherwise I wouldn't have noticed this. 

Back to the calculations below is a summary tbh, I wrote these down on the comments of the code, so since I am moving them here I will just delete from there.
```
    The idea is to go from the area of the triangle to the barycentric coordinates of the intersection point.
    This will allow me to early quit with some checks as well.

    We can sit down and write cramer's rule function or determinant function to solve the equation.
    HOWEVER, I want speed and I won't be calculating the bigger determinants or I won't need cramer for other things.
    Instead, I will just plug in the formula for alpha, beta and gamma in the derived form.

    Taken from lecture slides:
    beta = | ax-ox ax-cx dx |
           | ay-oy ay-cy dy |
           | az-oz az-cz dz |
           -------------------
                  |area|

    gamma = | ax-bx ax-ox dx |
           | ay-by ay-oy dy |
           | az-bz az-oz dz |
           -------------------
                  |area|

    alpha = 1 - beta - gamma

    determinantBeta = (ax - ox) * [(ay-cy)*dz - (az-cz)*dy] 
                    + (ay - oy) * [(az-cz)*dx - (ax-cx)*dz]
                    + (az - oz) * [(ax-cx)*dy - (ay-cy)*dx]
    

    determinantGamma = (ax - bx) * [(ay-oy)*dz - (az-oz)*dy] 
                    + (ay - by) * [(az-oz)*dx - (ax-ox)*dz]
                    + (az - bz) * [(ax-ox)*dy - (ay-oy)*dx]
                    
    beta = determinantBeta / |area|
    gamma = determinantGamma / |area|
    alpha = 1 - beta - gamma

    For the t calculation we have (replacing first column with o-a):

    determinantT = | ax-bx ax-cx ax-ox |
                   | ay-by ay-cy ay-oy |
                   | az-bz az-cz az-oz |

    determinantT = (ax-bx) * [(ay-cy)*(az-oz) - (az-cz)*(ay-oy)]
                 + (ay-by) * [(az-cz)*(ax-ox) - (ax-cx)*(az-oz)]
                 + (az-bz) * [(ax-cx)*(ay-oy) - (ay-cy)*(ax-ox)]

    t = determinantT / |area|
```

Now that `utils.cpp` is bloated with precomputations, main logic and operator overloads I will just separate them into smaller files.

```
overloads.h - for operator overloads
precomputations.h - for precomputations
utils.h - for the main logic
utils.cpp - for the implementation of the main logic
```

The simple.json scene file is almost rendering correctly now. I’ll take a break and come back later tonight to figure out why it’s only nearly correct and then continue with the other scenes. It seems that some parsing issues are still present — probably because I skipped some of the additional guidelines at the beginning to get my raytracer running faster. I guess I’m paying for that shortcut now.

## 2025-10-27

I fixed the **almost correct** rendering of the simple.json scene file. It was because of the intersection test epsilon that I had been using for the `beta` and `gamma` checks. Now moving on to the other scenes.

I still need to implement several things which I put below in order of priority, I am starting with the parsing issue first then comes the shadows.
- [x] Fix rendering bug in meshes
   It looks like as soon as I find the intersection in the mesh check I was returning true and that caused my creepy looking bunny :D I fixed it in this commit.
- [x] Parsing issue with other scenes
       This was a silly mistake. Some materials have conditional fields and I was trying to parse them everytime. Now I am conditionally parsing them with the `contains()` function.
- [x] Shadow rays
       Update: I implemented the shadow rays and it works! It was really straightforward considering the fact that I just used the `intersect()` function to check with the new shadow ray. Also I updated the `Ray` struct to include the shadow ray flag. And now I am using it inside the `rayHitsTriangle` to whether to calculate the determinant on the fly or not. 
- [x] Reflection
       Update: I implemented the reflection. Once I get the layout right, adding these features on top was easy for me tbh. My raytracer now works with the `spheres_mirror.json` scene file and it mostly looks correct. However, there is one tiny difference that the big ball where it touches the plane doesn't have enough whiteness as the ground truth. Therefore, I gotta take alook at this after I finish this completely. There might be some tiny floating point issues OR it could be something more serious.
- [ ] Refraction
- [ ] Performance optimizations -- I am thinking of threading, if I have time I will experiment with the acceleration structures.