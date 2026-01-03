all:
	g++ *.cpp -O3 -o raytracer -std=c++11 -lm -lpthread -w -I.

gui:
	g++ *.cpp -O3 -o raytracer -std=c++11 -lm -lpthread -w `pkg-config --cflags --libs sdl3`

clean:
	rm raytracer

clear:
	rm *.ppm
