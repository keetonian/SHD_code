#include <stdio.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>

#include "ctpl.h"

using namespace std;

char* ParseFile(char** argv, int i);
void PrepareReference16(char* ref_file, vector<unsigned short> * ref);
void PrepareReference4(char* ref_file, vector<unsigned char> * ref);
unsigned short ConvertCharacters16(char char1, char char2);
unsigned short ConvertInverseCharacters16(char char1, char char2);
void read_compare_func16(int id, string header, string read_line, vector<unsigned short> * ref, int shift, int threshold);
void read_compare_func4(int id, string header, string read_line, vector<unsigned char> * ref, int shift, int threshold);
void CompareRead16(char* read_file, char* reference_file, int shift, int threshold);
void CompareRead4(char* read_file, char* reference_file, int shift, int threshold);
unsigned char ConvertCharacter4(char char1);
unsigned char ConvertCharacterInverse4(char char1);

string compare16(vector<unsigned short>* ref, vector<unsigned short> read_in, vector<unsigned short> inverse, int shift, int threshold, int* num_matches);
string compare4(vector<unsigned char>* ref, vector<unsigned char> read_in, vector<unsigned char> inverse, int shift, int threshold, int* num_matches);

static int report_total_matches = 0;

int main(int argc, char** argv)
{
  if(argc < 2) {
    cerr << "Error: no arguments specified.\n";
    return 1;
  }

  char* reference_file = 0;
  char* read_file = 0;

  // Set some default parameters
  int shift = 0;
  int encoding = 0;
  int threshold = 20;
  int i = 1;
  for(; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
        case 'g': //reference genome
          {
            reference_file = ParseFile(argv, i);
            if(reference_file == 0) {
              cerr << "Error parsing reference file.\n";
              return 1;
            }
            // Debugging output, remove later.
            cout << "Reference: " << reference_file << "\n";
          }
          break;
        case '4': encoding = 0;
                  break;
        case '1': encoding = 1;
                  break;
        case 'm': report_total_matches = 1;
                  break;
        case 'h': //help
                  {
                    cout << "This program compares a read sequence with a reference genome.\n";
                    cout << "\t-g: reference genome path\n";
                    cout << "\t-r: read sequence path\n";
                    cout << "\t-h: display this help file\n";
                    cout << "\t-t: tolerated error threshold\n";
                    cout << "\t\t0: only exact matches (no errors tolerated)\n";
                    cout << "\t-s: shift distance\n";
                    cout << "\t-4: use 4-bit encodings\n";
                    cout << "\t-16: use 16-bit encodings\n";
                    cout << "\t-m: turn on total match reporting\n";
                    // Using the help flag will only print this help dialog.
                    return 0;
                  }
                  break;
        case 'r': //read file
                  {
                    read_file = ParseFile(argv, i);
                    if(read_file == 0) {
                      cerr << "Error parsing read file.\n";
                      return 1;
                    }
                    // Debugging output, remove later.
                    cout << "Read: " << read_file << "\n";
                  }
                  break;
        case 't': //threshold (read size - threshold) 0 means no errors tolerated
                  {
                    if(argv[i][2] != 0)
                      threshold = atoi(argv[i]+2);
                    else
                      threshold = atoi(argv[i+1]);
                  }
                  break;
        case 's': //shift distance
                  {
                    if(argv[i][2] != 0)
                      shift = atoi(argv[i]+2);
                    else
                      shift = atoi(argv[i+1]);
                  }
                  break;
        default: //unknown flag
                  {
                    cerr << "Error: illegal flag: " << argv[i] << "\n";
                    cerr << "Use -h to print help dialog.\n";
                    return 1;
                  }
                  break;
      }
    }
  }

  cout<<"shift: "<< shift<<endl;
  cout<<"threshold: "<< threshold<<endl;

  if(!reference_file || !read_file) {
    cerr<< "Need to specify a reference and read file.\n";
    return 1;
  } 

  /* It would be worth it to test this code to make sure it:
   *    Compares each step of the reference only 1 time (no more, no less)
   *    Accurately returns the number of matches.
   */

  if(encoding){
    cout<<"16 Bit encodings.\n";
    CompareRead16(read_file, reference_file, shift, threshold);
  }
  else{
    cout<<"4 Bit encodings.\n";
    CompareRead4(read_file, reference_file, shift, threshold);
  }

  // This doesn't yet test if multiple of the same flag are used, or other cases.
  // The program still needs to check the file extensions and types of the read and ref given.

  return 0;
}

/*
 * Parses a file path from the command line arguments
 */
char* ParseFile(char** argv, int i) {
  char* file_name = 0;
  if(argv[i][2] == 0) {
    //Means that the file should be in the next argument.
    file_name = argv[i+1];
    if(file_name[0] == '-') {
      return 0;
    }
  }
  else {
    //Means that the file should be part of the current argument.
    file_name = argv[i] + 2;
  }
  return file_name;
}

/*
 * Compares each read in the read_file to the reference genome according to the shift and threshold
 * 16 bit version
 */
void CompareRead16(char* read_file, char* reference_file, int shift, int threshold) {
  std::vector<unsigned short> ref(0);
  ref.reserve(3500000000);
  PrepareReference16(reference_file, &ref);
  ctpl::thread_pool p(8);

  // Set up and read from the file with the read sequences
  string read_line;
  ifstream file(read_file);
  if(file.is_open()) {

    // Step through the read sequences
    while(getline(file, read_line)) {
      if(read_line.size() > 2 && read_line[0] == '@') {
        char c = read_line[0];
        string s = "";
        int cindex = 0;
        while(c != ' ') {
          s+=c;
          cindex++;
          c = read_line[cindex];
        }
        s+='\t';

        getline(file, read_line);
        p.push(read_compare_func16, s, read_line, &ref, shift, threshold);
      }
    }
  }
  else {
    cerr << "Unable to open reference file \n";
    exit(1);
  }
} 

/*
 * Helper function for CompareRead16. This function is passed into thread arguments
 */
void read_compare_func16(int id, string header, string read_line, vector<unsigned short> * ref, int shift, int threshold) {
  int j = 0;
  // Save read sequence
  header += read_line + '\n';

  std::vector<unsigned short> readv(0);
  std::vector<unsigned short> read_inverse(0);
  int read_size = read_line.size();
  for(; j < read_size-1; j++) {
    readv.push_back(ConvertCharacters16(read_line[j], read_line[j+1]));
    read_inverse.push_back(ConvertInverseCharacters16(read_line[read_size-(j+1)], read_line[read_size-(j+2)]));
  }

  int num_matches = 0;
  header += compare16(ref, readv, read_inverse, shift, threshold, &num_matches);

  if(report_total_matches)
    cout<<header<<"TOTAL MATCHES: " << num_matches << '\n' << endl;
  else
    cout<<header<<endl;
}

/*
 * Compares each read in the read_file to the reference genome according o the shift and threshold
 * 4 bit version
 */
void CompareRead4(char* read_file, char* reference_file, int shift, int threshold) {
  std::vector<unsigned char> ref(0);
  ref.reserve(3500000000);
  PrepareReference4(reference_file, &ref);
  ctpl::thread_pool p(8);

  // Set up and read from the file with the read sequences
  string read_line;
  ifstream file(read_file);
  if(file.is_open()) {

    // Step through the read sequences
    while(getline(file, read_line)) {
      if(read_line.size() > 2 && read_line[0] == '@') {
        char c = read_line[0];
        string s = "";
        int cindex = 0;
        while(c != ' ') {
          s+=c;
          cindex++;
          c = read_line[cindex];
        }
        s+='\t';

        getline(file, read_line);
        p.push(read_compare_func4, s, read_line, &ref, shift, threshold);
      }
    }
  }
  else {
    cerr << "Unable to open reference file \n";
    exit(1);
  }
} 

/*
 * Helper function for CompareRead4. THis function is passed into threads to be run
 */
void read_compare_func4(int id, string header, string read_line, vector<unsigned char> * ref, int shift, int threshold) {
  int j = 0;
  // Save read sequence
  header += read_line + '\n';

  std::vector<unsigned char> readv(0);
  std::vector<unsigned char> read_inverse(0);
  int read_size = read_line.size();
  for(; j < read_size-1; j++) {
    readv.push_back(ConvertCharacter4(read_line[j]));
    read_inverse.push_back(ConvertCharacterInverse4(read_line[read_size-(j+1)]));
  }

  int num_matches = 0;
  header += compare4(ref, readv, read_inverse, shift, threshold, &num_matches);

  if(report_total_matches)
    cout<<header<<"TOTAL MATCHES: " << num_matches << '\n' << endl;
  else
    cout<<header<<endl;
}

/*
 * Opens the reference file, converts it to 16-bit encodings.
 * This method is very memory-hungry. It attemps to load the entire reference.
 * AA: 0x8000 
 * AT: 0x4000
 * AC: 0x2000
 * AG: 0x1000
 * TA: 0x0800
 * TT: 0x0400
 * TC: 0x0200
 * TG: 0x0100
 * CA: 0x0080
 * CT: 0x0040
 * CC: 0x0020
 * CG: 0x0010
 * GA: 0x0008
 * GT: 0x0004
 * GC: 0x0002
 * GG: 0x0001
 */
void PrepareReference16(char* ref_file, vector<unsigned short> * ref){
  string line;
  ifstream reference(ref_file);
  if(reference.is_open()) {
    char previous = 0;
    // What if we tried buffers instead of lines?
    while(getline(reference, line)) {
      if(line.size() ==0)
        continue;
      if (line[0] != '>'){
        unsigned int i = 0;
        if(previous)
          ref->push_back(ConvertCharacters16(previous, line[0]));
        for(; i < line.size()-1; i++) {
          ref->push_back(ConvertCharacters16(line[i], line[i+1]));
        }
        previous = line[line.size()-1];
      }
      else {

      }
    }

    reference.close();
  }
  else {
    cerr << "Unable to open reference file\n";
    exit(1);
  }
}

/*
 * Converts 2 characters to their 16-bit encodings
 * */
unsigned short ConvertCharacters16(char char1, char char2) {
  if(char1 == 'A') {
    switch(char2) {
      case 'A': return 0x8000;
      case 'T': return 0x4000;
      case 'C': return 0x2000;
      case 'G': return 0x1000;
      default: return 0;
    }
  }
  else if(char1 == 'T') {
    switch(char2) {
      case 'A': return 0x0800;
      case 'T': return 0x0400;
      case 'C': return 0x0200;
      case 'G': return 0x0100;
      default: return 0;
    }
  }
  else if(char1 == 'C') {
    switch(char2) {
      case 'A': return 0x0080;
      case 'T': return 0x0040;
      case 'C': return 0x0020;
      case 'G': return 0x0010;
      default: return 0;
    }
  }
  else if(char1 == 'G') {
    switch(char2) {
      case 'A': return 0x0008;
      case 'T': return 0x0004;
      case 'C': return 0x0002;
      case 'G': return 0x0001;
      default: return 0;
    }
  }
  else
    return 0;
}

/*
 * Converts 2 characters to their 16-bit encodings
 * */
unsigned short ConvertInverseCharacters16(char char1, char char2) {
  if(char1 == 'T') {
    switch(char2) {
      case 'T': return 0x8000;
      case 'A': return 0x4000;
      case 'G': return 0x2000;
      case 'C': return 0x1000;
      default: return 0;
    }
  }
  else if(char1 == 'A') {
    switch(char2) {
      case 'T': return 0x0800;
      case 'A': return 0x0400;
      case 'G': return 0x0200;
      case 'C': return 0x0100;
      default: return 0;
    }
  }
  else if(char1 == 'G') {
    switch(char2) {
      case 'T': return 0x0080;
      case 'A': return 0x0040;
      case 'G': return 0x0020;
      case 'C': return 0x0010;
      default: return 0;
    }
  }
  else if(char1 == 'C') {
    switch(char2) {
      case 'T': return 0x0008;
      case 'A': return 0x0004;
      case 'G': return 0x0002;
      case 'C': return 0x0001;
      default: return 0;
    }
  }
  else
    return 0;
}

/*
 * Converts the reference into 4 bit codes
 */
void PrepareReference4(char* ref_file, vector<unsigned char> * ref){
  string line;
  ifstream reference(ref_file);
  if(reference.is_open()) {

    // What if we tried buffers instead of lines?
    while(getline(reference, line)) {
      if(line.size() ==0)
        continue;
      if (line[0] != '>'){
        unsigned int i = 0;
        for(; i < line.size(); i++) {
          ref->push_back(ConvertCharacter4(line[i]));
        }
      }
      else {

      }
    }

    reference.close();
  }
  else {
    cerr << "Unable to open reference file\n";
    exit(1);
  }
}


unsigned char ConvertCharacters4(char char1, char char2) {
  return (ConvertCharacter4(char1) << 4) + ConvertCharacter4(char2);
}

unsigned char ConvertCharacter4(char char1) {
  switch(char1) {
    case 'A': return 0x8;
    case 'T': return 0x4;
    case 'C': return 0x2;
    case 'G': return 0x1;
    default: return 0;
  }
}

unsigned char ConvertCharacterInverse4(char char1) {
  switch(char1) {
    case 'T': return 0x8;
    case 'A': return 0x4;
    case 'G': return 0x2;
    case 'C': return 0x1;
    default: return 0;
  }
}



/**
 * Function to do the 16-bit comparison.
 * ref: pointer to reference genome (large compared to read). Format: 16-bit 1-hot encodings.
 * read: pointer to read sequence. Again, 16-bit 1-hot encodings.
 * shift: multiplex distance
 * threshold: minimum value allowed.
 */
string compare16(vector<unsigned short>* ref, vector<unsigned short> read_in, vector<unsigned short> inverse, int shift, int threshold, int* num_matches) {

  // Set up space for results.
  stringstream s;
  int matches = 0;


  // Set up variables
  unsigned int i = 0;
  unsigned int k = 0;
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
    unsigned int j;
    for(j=0; j <= (unsigned)shift; j++) {
      if(i >= j){
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
    int error = 0;
    int error_inv = 0;
    for(i = 0; i < read_size; i++) {
      // Prepare the read parameters
      unsigned short a = read.at(i);
      unsigned short b = inv.at(i);
      unsigned short r = ref->at(k+i);

      // Do the comparison
      error+=(!(a&r));
      error_inv+=(!(b&r));
      if(error > threshold && error_inv > threshold) {
        break;
      }

    }
    //If the result is over or equal to the threshold
    if(error <= threshold){
      unsigned int kk = 0;
      s << k << "\tn\t" << error << '\t';
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
        s << c2 << c3;
      }
      s<<"\n";
      matches += 1;
    }

    if(error_inv <= threshold){
      unsigned int kk = 0;
      s << k << "\tri\t" << error_inv << '\t';
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
        s << c2 << c3;
      }
      s << "\n";
      matches += 1;
    }
  }

  (*num_matches) = matches;

  return s.str();
}


/**
 * Function to do the 4-bit comparison.
 * ref: pointer to reference genome (large compared to read). Format: 16-bit 1-hot encodings.
 * read: pointer to read sequence. Again, 4-bit 1-hot encodings.
 * shift: multiplex distance
 * threshold: minimum value allowed.
 */
string compare4(vector<unsigned char>* ref, vector<unsigned char> read_in, vector<unsigned char> inverse, int shift, int threshold, int* num_matches) {

  // Set up space for results.
  stringstream s;
  int matches = 0;

  // Set up variables
  unsigned int i = 0;
  unsigned int k = 0;
  unsigned int read_size = read_in.size();
  unsigned int ref_size = ref->size();
  vector<unsigned char> read(0);
  vector<unsigned char> inv(0);
  inv.reserve(read_size);
  read.reserve(read_size);

  // Prepare the forward and backwards reads
  for(i = 0; (unsigned)i < read_size; i++) {
    unsigned char a = read_in.at(i);
    unsigned char b = inverse.at(i);
    unsigned int j; 
    for(j= 0; j <= (unsigned)shift; j++) {
      if(i >= j){
        a = a | read_in.at(i-j);
        b = b | inverse.at(i-j);
      }
      if((i+j)<read_size){
        a = a | read_in.at(i+j);
        b = b | inverse.at(i+j);
      }
    }
    read.push_back(a);
    inv.push_back(b);
  }

  // Compare read strand with entire reference genome.
  for(k = 0; k <= ref_size - read_size; k++) {
    int error = 0;
    int error_inv = 0;

    for(i = 0; i < read_size; i++) {
      // Prepare the read parameters
      unsigned char a = read.at(i);
      unsigned char b = inv.at(i);
      unsigned char r = ref->at(k+i);

      // Do the comparison
      error += (!(a&r));
      error_inv += (!(b&r));
      if(error > threshold && error_inv > threshold) {
        break;
      }
    }

    //If the result is over or equal to the threshold
    if(error <= threshold){
      unsigned int kk = 0;
      s << k << "\tn\t" << error << '\t';
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
        }
        s << c2;
      }
      s << "\n";
      matches += 1;
    }

    //If the result is over or equal to the threshold
    if(error_inv <= threshold){
      unsigned int kk = 0;
      s << k << "\tn\t" << error_inv << '\t';
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
        }
        s << c2;
      }
      s << "\n";
      matches += 1;
    }
  }

  (*num_matches) = matches;

  return s.str();
}
