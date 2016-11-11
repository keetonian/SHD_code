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
static std::vector<long> compare16(vector<unsigned short>* ref, vector<unsigned short>* read, int shift, int threshold) {

  // Compare the read against every spot in the reference
  // Return an array of locations in the reference that are viable
  // Alternatively could return an array of strings of viable comparison spots.

  // Set up space for results.
  vector<long> results(0);

  // Set up variables
  int i = 0;
  long k = 0;
  unsigned int read_size = read->size();
  unsigned int ref_size = ref->size();

  // Compare read strand with entire reference genome.
  for(k = 0; k <= ref_size - read_size; k++) {
    int result = 0;
    for(i = 0; (unsigned)i < read_size; i++) {

      // Prepare the read parameters
      unsigned short a = read->at(i);
      int j;
      for(j = 0; j <= shift; j++) {
        if(i-j>=0)
          a = a | read->at(i-j);
        if((unsigned)(i+j)<read_size)
          a = a | read->at(i+j);
      }

      // Do the comparison
      if(a & ref->at(k+i))
        result++;

    }
    //If the result is over or equal to the threshold
    if(result >= (int)read_size - threshold)
      results.push_back(k);
  }

  return results;
}
