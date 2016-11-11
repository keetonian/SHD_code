#include <stdio.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <fstream>

#include "compare16.cpp"

using namespace std;

char* ParseFile(char** argv, int i);
void PrepareReference16(char* ref_file, vector<long> * chromosomes);
unsigned short ConvertCharacters16(char char1, char char2);
static std::vector<unsigned short> * ref = new vector<unsigned short>(0);
static std::vector<unsigned short> * read = new vector<unsigned short>(0);

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
  int threshold = 90;
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
            cout << "\t\t0: only exact matches (no errors tolerated)
            cout << "\t-s: shift distance\n";
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

  if(!reference_file || !read_file) {
    cerr<< "Need to specify a reference and read file.\n";
    return 1;
  }

  // ### It would be good to split the rest of main into another function ### \\
  // It would also be worth it to test this code to make sure it:
  //    Compares each step of the reference only 1 time (no more, no less)
  //    Accurately returns the number of matches.

  // Not used yet.
  vector<long> chromosomes();

  // Set up and read from the file with the read sequences
  string read_line;
  ifstream file(read_file);
  if(file.is_open()) {
    long linec = 0;

    // Set up and read from the reference genome file
    string ref_line, ref_line2;
    ifstream reference(reference_file);
    if(reference.is_open()){
      char previous = 0;

      // Step through the read sequences
      while(getline(file, read_line)) {
        if(((linec - 1) % 4) == 0){
          int j = 0;

          // Print the current read sequence
          printf("%s\n", read_line.c_str());

          // Convert the read into the 16 bit encodings
          std::vector<unsigned short> readv(0);
          if(previous)
            readv.push_back(ConvertCharacters16(previous, read_line[0]));
          for(; j < read_line.size()-1; j++) {
            readv.push_back(ConvertCharacters16(read_line[j], read_line[j+1]));
          }
          previous = read_line[read_line.size() -1];

          // Step through the reference and do comparison
          getline(reference, ref_line); // first line is chromosomal information
          getline(reference, ref_line);
          int matches = 0;
          unsigned int loc = 1;
          while(getline(reference, ref_line2)) {
            if(ref_line2[0] == '>') // contains chromosomal information. skip.
              continue;
            string line = ref_line + ref_line2;
            
            //Convert the reference into the 16 bit encodings.
            vector<unsigned short> refv(0);
            int i = 0;
            for(; i < string.size()-1; i++) {
              refv.push_back(ConvertCharacters16(line[i], line[i+1]));
            }
            ref_line = ref_line2;
            std::vector<long> * v = compare16(&refv, &readv, shift, threshold);
            
            // Loc only tracks how many lines of the reference genome have been read
            // It does not count chromosomal information lines.
            loc++;

            // If a match or matches were found, print the entire string.
            // This will later turn to printing either a location or just the (read sized) reference
            if(v->size()){
              printf("%s\n",line.c_str());
              printf("%d\n", loc);
            }

            // Increment the total matches
            matches += v->size();
          }

          // Return to the beginning of the reference file
          reference.seekg(0, reference.beg);

          // Print total number of matches
          printf("TOTAL MATCHES: %d\n", matches);
        }
        linec++;
      }
    }
    else {
      cerr << "Unable to open reference file \n";
      file.close();
      return 1;
    }
  }
  else {
    cerr << "Unable to open reference file \n";
    return 1;
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
void PrepareReference16(char* ref_file, vector<long> * chromosomes){
  string line;
  ifstream reference(ref_file);
  if(reference.is_open()) {
    char previous = 0;

    // What if we tried buffers instead of lines?
    while(getline(reference, line)) {
      if (line[0] != '>'){
        int i = 0;
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


// 4-bit encodings. We need to decide how we would like to do these.
unsigned char ConvertCharacter4(char char1);

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

