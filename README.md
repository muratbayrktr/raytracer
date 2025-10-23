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
