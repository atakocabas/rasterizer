all:
	g++ *.cpp -o rasterizer -std=c++11
clean:
	rm -rf *.ppm ./rasterizer