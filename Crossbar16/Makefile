CXX = g++
CFLAGS = -std=c++11 -O2 -pthread -lpthread -lboost_thread

Crossbar16: readfasta.cpp compare16.cpp compare4.cpp
	$(CXX) $(CFLAGS) -o Crossbar16 readfasta.cpp compare16.cpp compare4.cpp

clean:
	rm Crossbar16 
