all:
	g++ *.cpp -O3 -o raytracer -std=c++11 -lm -lpthread 

clean:
	rm raytracer

clear:
	rm *.ppm
