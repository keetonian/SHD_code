CC=gcc
CFLAGS = -c -O3 -Wall -msse -msse2 -I . -I/usr/local/boost_1_62_0 -msse4.2
LDFLAGS = -lz -lm 
SOURCES = baseFAST.c CommandLineParser.c Common.c HashTable.c MrFAST.c Output.c Reads.c RefGenome.c
OBJECTS = $(SOURCES:.c=.o) ../SHD/vector_filter.o ../SHD/bit_convert.o ../SHD/mask.o ../SHD/popcount.o ../SHD/print.o 
EXECUTABLE = mrfast


all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

.c.o:
	$(CC) $(CFLAGS) $< -o $@ 
clean:
	rm -f *.o *~ \#* mrfast
