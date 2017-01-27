TODO:
Modify SHD for our algorithms
1. Modify SHD for the 4-bit encoding algorithm that uses an error threshold at math.pow(2, log2(read_size)-1). In this case, a larger number means more accuracy, whereas the traditional SHD algorithm has a lower number mean greater accuracy.
2. Modify SHD for the 16-bit encoding algorithm.

Graph:
  Filters:
    Speed for GateKeeper, SHD, Ours
    # reads for GateKeeper, SHD, Ours

TODO:
Crossbar algorithm: output in SAM format (or at least report by chromosome+index)
Consolidate tools
Clean up repository
Write comparison code that takes the actual matches and checks for false pos/neg matches for each type of result
Sort matches, make them easier to access
Generate graphs, visuals of data
Continue reading papers on the subject, prepare a presentation.

