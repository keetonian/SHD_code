#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <stdlib.h>
#include "vector_filter.h"

using namespace std;
void helper_func(int id, vector<char> * ref, string line, string header, int error);


int main(int argc, char** argv) {
  if(argc != 4){
    cout<<"Usage: ./prepareinput <error> <read_file> <mrfast_file>"<<endl;
    cerr<<"Wrong number of arguments"<<endl;
    return 0;
  }

  //TODO: use a read file with mrfast file

  int error = atoi(argv[1]);
  char* ref_file = argv[2];
  char* read_file = argv[3];

  vector<char> ref(0);
  ref.reserve(32000);

  string line;
  ifstream reference(ref_file);
  if(reference.is_open()) {
    while(getline(reference, line)) {
      if(line.size() == 0)
        continue;
      if(line[0] != '>') {
        int i = 0;
        for(; i < line.size(); i++) {
          ref.push_back(line[i]);
        }
      }
    }
    reference.close();
  }
  else {
    cerr << "Unable to open reference file" << endl;
    exit(1);
  }

  long ref_size = ref.size();
  ifstream read(read_file);
  if(read.is_open()) {
    while(getline(read, line)) {
      if(line.size() > 2 && line[0] == '@' && line[1] == 'E') {
        char c = line[0];
        string s = "";
        int cindex = 0;
        while(c != ' ') {
          s += c;
          cindex++;
          c = line[cindex];
        }
        getline(read, line);
        helper_func(1, &ref, line, s, error);
        //p.push(helper_func, &ref, line, s, error);
      }
    }
    read.close();
  }
  else {
    cerr<<"Error in reading read file"<<endl;
  }
  return 0;
}


void helper_func(int id, vector<char> * ref, string line, string header, int error) {

  char read_t[128] __aligned__;
  char ref_t[128] __aligned__;
  char init_all_NULL[128] = "";
  strncpy(read_t, init_all_NULL, 128);
  strncpy(ref_t, init_all_NULL, 128);

  long ref_size = ref->size();
  header += line + '\n';
  int read_size = line.size();
  strncpy(read_t, line.c_str(), read_size);
  long i = 0;
  for(; i < ref_size-(read_size+1); i++) {
    //cout << line << endl;
    for(int j = 0; j < read_size; j++) {
      ref_t[j] = ref->at(i+j);
    }
    ref_t[read_size] = '\0';
    if(bit_vec_filter_sse1(read_t, ref_t, read_size, error)) {
      header += to_string(i) + '\n';
    }
  }
  cout << header << endl;
}
