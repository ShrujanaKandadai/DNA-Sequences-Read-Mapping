#Submitted by Shrujana Kandadai
import argparse 
import os
import matplotlib
from matplotlib import pyplot as plt

def load_reference_fasta(filepath):
    '''Reads a reference FASTA file and returns a dictionary {ref_id: sequence}.'''
    lines= []
    ref_sequences = {}
    with open(filepath, 'r') as f:
        for line in f:
            lines.append(line.strip())
        
        for i in range(len(lines)):
            if lines[i][0] == '>':
                ref_id = lines[i][1:]
                sequence = lines[i+1]
                ref_sequences[ref_id] = sequence

    return ref_sequences

def build_ref_kmer_index(ref_sequences,k):
    """Builds a dictionary {kmer: [(ref_id, position), ...]} for the reference."""
    read_kmers = {}
    t1 = ()
    for key, value in ref_sequences.items():
        for i in range(len(value)-k+1):
            kmer = value[i:i+k]
            if kmer in read_kmers:
                t1 = (key, i)
                read_kmers[kmer].append(t1)
            else:
                read_kmers[kmer] = [(key, i)]
    return read_kmers

def load_reads_fasta(file):
    """Reads a reads FASTA file and returns a dictionary {read_id: sequence}."""
    sequences= []
    reads = {}
    with open(file, 'r') as f:
        for line in f:
            sequences.append(line.strip())
        for i in range(len(sequences)):
            if sequences[i][0] == '>':
                key = sequences[i][1:]
                value = sequences[i+1]
                reads[key] = value
    return reads
     
def build_read_kmers(reads, k):
    """Builds a dictionary {read_id: [(kmer, position), ...]} for all reads."""
    read_kmers = {}
    for key, value in reads.items():
        read_kmers[key] = []
        for i in range(len(value)-k+1):
            kmer = value[i:i+k]
            read_kmers[key].append((kmer, i))
    return read_kmers        

def hamming_distance(seq1, seq2):
    """Computes the number of mismatched bases between two equal-length sequences."""
    mismatches = 0
    # Both sequences are length 'r'
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mismatches += 1
    return mismatches

def map_reads(read_kmers, ref_sequences, ref_kmers, reads, k, out1_path):
    """Maps reads to reference sequences using k-mer matches and mismatch counting."""
    mapping_results = {} #Dictionary to store mapping results for each read
    
    with open(out1_path, 'w') as f_out:

        #Looping through each read and its list of (k-mer,position) pairs
        for read_id, kmer_list in read_kmers.items():
          read_seq = reads[read_id] #getting the full read 
          min_mismatches = float('inf') ## Initializing to track lowest mismatch count
          best_locations = [] #to store positions in reference with minimum mismatches
          checked_positions = set() #prevent re-checking the same locations in reference
          
          #iterating through each k-mer in the read 
          for kmer, read_pos in kmer_list:
              # Only proceeding if the k-mer exists in reference k-mer index
              if kmer in ref_kmers:
                   # Each reference k-mer may appear in multiple references and positions
                  for ref_id, ref_pos in ref_kmers[kmer]:
                      # Compute candidate alignment start position in the reference
                      start_pos = ref_pos - read_pos

                       # Skip if read would extend beyond reference boundaries
                      if start_pos < 0 or start_pos + len(read_seq)> len(ref_sequences[ref_id]):
                          continue
                      
                      # Skip if this (reference, start position) combination was already tested
                      if (ref_id, start_pos) in checked_positions:
                          continue
                      checked_positions.add((ref_id, start_pos))

                       # Extract corresponding segment from reference to compare with read
                      ref_segment = ref_sequences[ref_id][start_pos:start_pos + len(read_seq)]
                      
                      # Compute number of mismatched bases (Hamming distance)
                      mismatches = hamming_distance(read_seq, ref_segment)

                      # Update best match if fewer mismatches are found
                      if mismatches < min_mismatches:
                          min_mismatches = mismatches 
                          best_locations = [(ref_id, start_pos, mismatches)]

                      # If tie, keep both positions
                      elif mismatches == min_mismatches:
                          best_locations.append((ref_id, start_pos, mismatches))

           # Store best mapping locations for this read
          mapping_results[read_id]=best_locations

          for ref_id, pos, dist in best_locations: 
          # Writing all best mapping locations to output file
                f_out.write(f"{read_id}, {ref_id}, {pos}, {dist}\n")
    
    return mapping_results

def compute_reference_coverage(mapping_results, ref_sequences, reads, output_dir):
    """Computes coverage (read abundance) per reference and writes to output2.txt."""
    # Initializing coverage counts for each reference sequence
    coverage_counts = {}            
    for ref_id in ref_sequences:             
        coverage_counts[ref_id] = 0

    for read_id, best_locs in mapping_results.items():
        read_len = len(reads[read_id])
        # A read with multiple best locations contributes its length to each 
        for ref_id, pos, dist in best_locs:
            coverage_counts[ref_id] += read_len

    out2_path = os.path.join(output_dir, "output2.txt")
    with open(out2_path, 'w') as f:
        for ref_id, total_bases in coverage_counts.items():
            ref_len = len(ref_sequences[ref_id])
            # Coverage = total nucleotides of mapped reads / sequence length 
            coverage = total_bases / ref_len
            f.write(f"{ref_id},{coverage}\n")

def plot_mapping_histogram(mapping_results, output_dir):
    """Plots a histogram showing how many best-mapping locations each read has."""
    # Get the count of best locations for every read
    counts = []
    for locs in mapping_results.values():
        if locs:
            counts.append(len(locs))
    
    plt.figure()
    # Using 50 bins as suggested
    plt.hist(counts, bins=50)
    plt.xlabel("Nr. best locations") 
    plt.ylabel("Count") 
    plt.title("Number of Best Mapping Locations per Read")
    
    
    # Saving as PDF with an informative filename 
    plt.savefig(os.path.join(output_dir, "mapping_histogram.pdf"))

def main():
    parser = argparse.ArgumentParser(description="Short-read mapping program: ")
    parser.add_argument("-k", type = int, required = True, help = "k-mer size")
    parser.add_argument("ref_file", help = "Path to reference FASTA file")
    parser.add_argument("reads_file", help = "Path to reads FASTA file")
    parser.add_argument("output_dir", help = "Directory for output files")
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    try:
        print("Reading: ", args.ref_file)
        ref_data = load_reference_fasta(args.ref_file)

        print("Reading: ", args.reads_file)
        read_data = load_reads_fasta(args.reads_file)

        if not ref_data or not read_data:
            print("One of the files is empty or not in correct FASTA format.")
            return 
    
    except FileNotFoundError:
        print("Could not find one of the files, please check the paths.")
        return 
    
    except Exception as e:
        print("An unexpected error occured", e)
        return 

    print("Indexing the reference and reads with k = ", args.k)
    ref_kmers = build_ref_kmer_index(ref_data, args.k)
    read_kmers = build_read_kmers(read_data, args.k)

    print("Mapping reads to reference ")
    out1_path = os.path.join(args.output_dir, "output1.txt")
    mapping_results = map_reads(read_kmers, ref_data, ref_kmers, read_data, args.k, out1_path)

    print("Calculating reference coverage ")
    compute_reference_coverage(mapping_results, ref_data, read_data, args.output_dir)

    print("Generating histogram")
    plot_mapping_histogram(mapping_results, args.output_dir)

    print("Processing completed. Results saved in ", args.output_dir)

if __name__ == "__main__":
    main()