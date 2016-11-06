long* compare(char* ref, char* read) {
  // Compare the read against every spot in the reference
  // Return an array of locations in the reference that are viable
  // Alternatively could return an array of strings of viable comparison spots.
  //

  int i = 0;
  int read_length = 128; //Find a more dynamic way of determining this.
  while(*(ref + read_length) != 0){
    for(i = 0; i < read_length; i++) {
      // Compare the genomes

    }
    ref++;
  }
  return 0;
}
