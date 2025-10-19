all:
	g++ *.cpp -O3 -o raytracer -std=c++11 -lm -I/opt/homebrew/include

clean:
	rm raytracer

clear:
	rm *.ppm
