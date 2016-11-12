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
static std::vector<long> compare4(vector<unsigned char>* ref, vector<unsigned char> read, int shift, int threshold) {

  // Compare the read against every spot in the reference
  // Return an array of locations in the reference that are viable
  // Alternatively could return an array of strings of viable comparison spots.

  // Set up space for results.
  vector<long> results(0);

  // Set up variables
  int i = 0;
  long k = 0;
  unsigned int read_size = read.size();
  unsigned int ref_size = ref->size();

  // Compare read strand with entire reference genome.
  for(k = 0; k <= ref_size - read_size; k++) {
    int result = 0;
    int error = 0;

    for(i = 0; (unsigned)i < read_size; i++) {

      // Prepare the read parameters
      unsigned char a = read.at(i);
      int j;
      for(j = 0; j <= shift; j++) {
        if(i-j>=0)
          a = a | read.at(i-j);
        if((unsigned)(i+j)<read_size)
          a = a | read.at(i+j);
      }

      // Do the comparison
      if(a & ref->at(k+i))
        result++;
      else
        error++;
      if(error > threshold)
        break;

    }
    //If the result is over or equal to the threshold
    if(result >= (int)read_size - threshold){
      unsigned int kk = 0;
      printf("%lu:\n", k);
      for(; kk < read_size; kk++) {
        unsigned char c1 = ref->at(k+kk);
        char c2 = 0x20;
        switch(c1) {
          case 0x8: c2 = 'A'; break;
          case 0x4: c2 = 'T'; break;
          case 0x2: c2 = 'C'; break;
          case 0x1: c2 = 'G'; break;
          case 0x0: c2 = 'N'; break;
          default: c2 = ' '; break;
        printf("%c", c2);
        }
      }
      printf("\n");
      results.push_back(k);
		}
  }

  return results;
}
