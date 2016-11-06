CXX = g++
CFLAGS = -Wall -std=c++11

Crossbar16: readfasta.cpp compare16.o
	$(CXX) $(CFLAGS) -o Crossbar16 readfasta.cpp

compare16.o: compare16.cpp
	$(CXX) $(CFLAGS) -o compare16.o -c compare16.cpp

clean:
	rm Crossbar16 *.o
