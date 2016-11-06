#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>

#include "readfasta.hpp"
#include "compare16.hpp"

using namespace std;

int main(int argc, char** argv)
{
  if(argc < 2) {
    cerr << "Error: no arguments specified.\n";
    return 1;
  }

  char* reference_file;
  char* read_file;

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
