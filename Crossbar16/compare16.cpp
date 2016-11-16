#include <vector>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

/**
 * Function to do the 16-bit comparison.
 * ref: pointer to reference genome (large compared to read). Format: 16-bit 1-hot encodings.
 * read: pointer to read sequence. Again, 16-bit 1-hot encodings.
 * shift: multiplex distance
 * threshold: minimum value allowed.
 */
static std::vector<long> compare16(vector<unsigned short>* ref, vector<unsigned short> read_in, vector<unsigned short> inverse, int shift, int threshold) {

  // Compare the read against every spot in the reference
  // Return an array of locations in the reference that are viable
  // Alternatively could return an array of strings of viable comparison spots.

  // Set up space for results.
  vector<long> results(0);

  // Set up variables
  int i = 0;
  long k = 0;
  unsigned int read_size = read_in.size();
  vector<unsigned short> read(0);
  vector<unsigned short> inv(0);
  read.reserve(read_size);
  inv.reserve(read_size);
  unsigned int ref_size = ref->size();

  for(i=0; i < read_size; i++) {
    // Prepare the read parameters
    unsigned short a = read_in.at(i);
    unsigned short b = inverse.at(i);
    int j;
    for(j=0; j <= shift; j++) {
      if(i-j>=0){
        a = a | read_in.at(i-j);
        b = b | inverse.at(i-j);
      }
      if((unsigned)(i+j)<read_size){
        a = a | read_in.at(i+j);
        b = b | inverse.at(i+j);
      }
    }
    inv.push_back(b);
    read.push_back(a);
  }

  // Compare read strand with entire reference genome.
  for(k = 0; k <= ref_size - read_size; k++) {
    int result = 0;
    int error = 0;
    int error_inv = 0;
    for(i = 0; i < read_size; i++) {
      // Prepare the read parameters
      unsigned short a = read.at(i);
      unsigned short b = inv.at(i);
      unsigned short r = ref->at(k+i);

      // Do the comparison
      if(!(a & r)){
        error++;
      }
      if(!(b & r)) {
        error_inv++;
      }
      if(error > threshold && error_inv > threshold) {
        break;
      }

    }
    //If the result is over or equal to the threshold
    if(error <= threshold){
			unsigned int kk = 0;
			printf("%lu\tn\t%d\t", k, error);
			for(; kk < read_size; kk+=2) {
        unsigned short c1 = ref->at(k+kk);
        char c2 = 0x20;
        char c3 = 0x20;
        switch(c1) {
          case 0x8000: c2 = 'A'; c3 = 'A'; break;
          case 0x4000: c2 = 'A'; c3 = 'T'; break;
          case 0x2000: c2 = 'A'; c3 = 'C'; break;
          case 0x1000: c2 = 'A'; c3 = 'G'; break;
          case 0x0800: c2 = 'T'; c3 = 'A'; break;
          case 0x0400: c2 = 'T'; c3 = 'T'; break;
          case 0x0200: c2 = 'T'; c3 = 'C'; break;
          case 0x0100: c2 = 'T'; c3 = 'G'; break;
          case 0x0080: c2 = 'C'; c3 = 'A'; break;
          case 0x0040: c2 = 'C'; c3 = 'T'; break;
          case 0x0020: c2 = 'C'; c3 = 'C'; break;
          case 0x0010: c2 = 'C'; c3 = 'G'; break;
          case 0x0008: c2 = 'G'; c3 = 'A'; break;
          case 0x0004: c2 = 'G'; c3 = 'T'; break;
          case 0x0002: c2 = 'G'; c3 = 'C'; break;
          case 0x0001: c2 = 'G'; c3 = 'G'; break;
          default: c2 = ' '; c3 = ' '; break;
        }
				printf("%c%c", c2, c3);
			}
			printf("\n");
      results.push_back(k);
		}
		
		if(error_inv <= threshold){
			unsigned int kk = 0;
			printf("%lu\tri\t%d\t", k, error_inv);
			for(; kk < read_size; kk+=2) {
        unsigned short c1 = ref->at(k+kk);
        char c2 = 0x20;
        char c3 = 0x20;
        switch(c1) {
          case 0x8000: c2 = 'A'; c3 = 'A'; break;
          case 0x4000: c2 = 'A'; c3 = 'T'; break;
          case 0x2000: c2 = 'A'; c3 = 'C'; break;
          case 0x1000: c2 = 'A'; c3 = 'G'; break;
          case 0x0800: c2 = 'T'; c3 = 'A'; break;
          case 0x0400: c2 = 'T'; c3 = 'T'; break;
          case 0x0200: c2 = 'T'; c3 = 'C'; break;
          case 0x0100: c2 = 'T'; c3 = 'G'; break;
          case 0x0080: c2 = 'C'; c3 = 'A'; break;
          case 0x0040: c2 = 'C'; c3 = 'T'; break;
          case 0x0020: c2 = 'C'; c3 = 'C'; break;
          case 0x0010: c2 = 'C'; c3 = 'G'; break;
          case 0x0008: c2 = 'G'; c3 = 'A'; break;
          case 0x0004: c2 = 'G'; c3 = 'T'; break;
          case 0x0002: c2 = 'G'; c3 = 'C'; break;
          case 0x0001: c2 = 'G'; c3 = 'G'; break;
          default: c2 = ' '; c3 = ' '; break;
        }
				printf("%c%c", c2, c3);
			}
			printf("\n");
      results.push_back(k);
		}
  }

  return results;
}
