all:
	g++ *.cpp -o  rasterizer -std=c++11 -O3
clean:
	rm -rf *.ppm ./rasterizer