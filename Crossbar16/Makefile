CXX = g++
CFLAGS = -std=c++11 -O2 -Wall -pthread -lpthread 

CrossbarSHD: crossbar_shd.cpp ctpl.h
	$(CXX) $(CFLAGS) -o CrossbarSHD crossbar_shd.cpp ctpl.h

clean:
	rm CrossbarSHD 
