#include <stdio.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <fstream>

#include "compare16.cpp"

using namespace std;

char* ParseFile(char** argv, int i);
int PrepareData(char* ref_file, char* read_file);
void PrepareReference16(char* ref_file, vector<long> * chromosomes);
unsigned short ConvertCharacters16(char char1, char char2);
static std::vector<unsigned short> * ref = new vector<unsigned short>(400);
static std::vector<unsigned short> * read = new vector<unsigned short>(200);

int main(int argc, char** argv)
{
  if(argc < 2) {
    cerr << "Error: no arguments specified.\n";
    return 1;
  }

  char* reference_file = 0;
  char* read_file = 0;
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

  if(!reference_file || !read_file) {
    cerr<< "Need to specify a reference and read file.\n";
    return 1;
  }

  vector<long> chromosomes();

  string read_line;
  ifstream file(read_file);
  if(file.is_open()) {
    long linec = 0;
    string ref_line, ref_line2;
    ifstream reference(reference_file);
    if(reference.is_open()){
      char previous = 0;

      // #### NOTE: IT DOESN'T SEEM TO BE READING THIS YET ####

      while(getline(file, read_line)) {
        if(((linec - 1) % 4) == 0){
          // Convert the read
          int j = 0;
          vector<unsigned short> readv(read_line.size());
          if(previous)
            readv.push_back(ConvertCharacters16(previous, read_line[0]));
          for(; j < read_line.size()-1; j++) {
            readv.push_back(ConvertCharacters16(read_line[j], read_line[j+1]));
          }
          previous = read_line[read_line.size() -1];

          // Step through the reference and do comparison
          getline(reference, ref_line);
          getline(reference, ref_line);
          int size = ref_line.size() << 1;
          while(getline(reference, ref_line2)) {
            if(ref_line2[0] == '>')
              continue;
            string line = ref_line + ref_line2;
            vector<unsigned short> refv(size);
            int i = 0;
            for(; i < size-1; i++) {
              refv.push_back(ConvertCharacters16(line[i], line[i+1]));
            }
            
            ref_line = ref_line2;
            cout<<compare16(&refv, &readv, shift, threshold) << endl;
        }
        linec++;
        reference.seekg(0, reference.beg);
      }
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
      

  //PrepareReference16(reference_file, chromosomes);
  printf("Reference size:%lu\n", ref->size());
  // This doesn't yet test if multiple of the same flag are used, or other cases.


  // If the program gets to this point, it should have a file path to both the read and the ref.
  // The program still needs to check the file extensions and types of the read and ref given.
  // This is where the rest of the work will be done.
  //
  // Change the input from string to array of unsigned short:
  // PrepareData(reference_file, read_file);
  std::vector<long> *positions = compare16(ref, read, shift, threshold);
  
  for(i = 0; (unsigned)i < positions->size(); i++)
    cout << "Element " << i << ":" << positions->at(i) <<endl;

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
void PrepareReference16(char* ref_file, vector<long> * chromosomes){
  string line;
  ifstream reference(ref_file);
  if(reference.is_open()) {
    char previous = 0;
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

