# Simple DNA Sequence Read Mapper
### This project implements a short-read mapping program in Python that aligns sequencing reads to one or more reference sequences using a k-mer based seed and extend approach.
### The program takes a reference FASTA file and reads FASTA file as input, indexes the reference into k-mers and then maps each read by identifying matching k-mers in the reference.
### It also computes the mismatches(Hamming distance) between the read and potential reference regions. The program reports, for each read, the best mapping locations and mismatch counts, calculates coverage for each reference sequence, and produces a histogram showing the number of best mapping locations per read.

# Program Execution 
### The program is executed from the command line using Python and the argparse module for argument parsing. It requires three input paths and a k-mer size:
***python readmapper.py -k <k-mer size> <reference_file> <reads_file> <output_directory>***

***For example: python readmapper.py -k 5 ref.fasta reads.fasta results/***
### The program prints progress messages to the terminal and produces three output files in the specified output directory. 
* output1.txt – list of best mapping locations (read ID, reference ID, start position, mismatches).
* output2.txt – coverage values for each reference sequence.
* mapping_histogram.pdf – histogram of the number of best mapping locations per read
