#include <stdio.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <fstream>

#include "compare16.cpp"

using namespace std;

char* ParseFile(char** argv, int i);
int PrepareData(char* ref_file, char* read_file);
static std::vector<unsigned short> * ref = new vector<unsigned short>();
static std::vector<unsigned short> * read = new vector<unsigned short>();

int main(int argc, char** argv)
{
  if(argc < 2) {
    cerr << "Error: no arguments specified.\n";
    return 1;
  }

  char* reference_file;
  char* read_file;
  int shift = 0;
  int threshold = 126;
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
        case 'h': //help
          {
            cout << "This program compares a read sequence with a reference genome.\n";
            cout << "\t-g: reference genome path\n";
            cout << "\t-r: read sequence path\n";
            cout << "\t-h: display this help file\n";
            cout << "\t-t: tolerated error threshold\n";
            cout << "\t-s: shift distance\n";
            // Using the help flag will only print this help dialog.
            return 0;
          }
          break;
        case 'r': //read
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
        case 't': //threshold
          {
            if(argv[i][2] != 0)
              threshold = atoi(argv[i]+2);
            else
              threshold = atoi(argv[i+1]);
            printf("threshold: %d\n", threshold);
          }
          break;
        case 's': //shift distance
          {
            if(argv[i][2] != 0)
              shift = atoi(argv[i]+2);
            else
              shift = atoi(argv[i+1]);
            printf("shift: %d\n", shift);
          }
          break;
        default: //unknown flag
          {
            cerr << "Error: illegal flag: " << argv[i] << "\n";
            return 1;
          }
          break;
      }
    }
  }
  // This doesn't yet test if multiple of the same flag are used, or other cases.


  // If the program gets to this point, it should have a file path to both the read and the ref.
  // The program still needs to check the file extensions and types of the read and ref given.
  // This is where the rest of the work will be done.
  //
  // Change the input from string to array of unsigned short:
  PrepareData(reference_file, read_file);
  std::vector<long> *positions = compare16(ref, read, shift, threshold);
  
  for(i = 0; (unsigned)i < positions->size(); i++)
    cout << "Element: " << i << ":" << positions->at(i) <<endl;

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

int PrepareData(char* ref_file, char* read_file){
  // Need to fully implement this method. Right now using dummy data:
  for(int i = 0; i < 100000000; i++) {
    ref->push_back(1<<i);
  }

  for(int j = 0; j < 5; j++) {
    read->push_back(1<<j);
  }

  return 0;
}

/*
 * Opens the reference file, converts it to 16-bit encodings.
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
int PrepareReference16(char* ref_file, vector<long> * chromosomes){
  FILE * reference;
  if((reference = fopen(ref_file,"r")) == NULL){
    printf("Reference file does not exist.\n");
    exit(1);
  }
  else
    printf("Reference opened.\n");

  // Do the conversions here.
  // Store converted values into ref_file
  // Store locations of the starts of chromosomes into chromosomes

  if(reference) {
    if(fclose(reference)){
      printf("Reference file not closed.\n");
    }
  }
  return 0;
}
