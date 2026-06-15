# test_readmapper.py
#Submitted by Shrujana Kandadai
""" This test-readmapper is used to run a set of basic tests for the readmapper module.

    Tests included:
    1. Parsing reference genome and reads from FASTA files.
    2. Breaking sequences into k-mers and checking indexing.
    3. Calculating Hamming distance between sequences.
    4. Mapping reads to the reference genome and checking mismatches.
    5. Calculating coverage and verifying output files.

    This test uses a tiny synthetic reference and read set and writes outputs to a temporary directory which is automatically cleaned up.
"""
import os
import tempfile
from readmapper import (
    load_reference_fasta,
    build_ref_kmer_index,
    load_reads_fasta,
    build_read_kmers,
    hamming_distance,
    map_reads,
    compute_reference_coverage
)

def run_basic_tests():
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create tiny test reference and reads
        ref_file = os.path.join(tmpdir, "ref.fasta")
        reads_file = os.path.join(tmpdir, "reads.fasta")

        with open(ref_file, "w") as f:
            f.write(">ref1\nABCDEFABCDEF\n")
        with open(reads_file, "w") as f:
            f.write(">read1\nABCDEFABC\n")
            f.write(">read2\nCDEFABCD\n")

        # 1. Test parsing
        ref_data = load_reference_fasta(ref_file)
        read_data = load_reads_fasta(reads_file)
        assert "ref1" in ref_data
        assert "read1" in read_data

        # 2. Test k-mer indexing
        ref_kmers = build_ref_kmer_index(ref_data, 3)
        read_kmers = build_read_kmers(read_data, 3)
        assert "ABC" in ref_kmers
        assert len(read_kmers["read1"]) == len(read_data["read1"]) - 3 + 1

        # 3. Test Hamming distance
        assert hamming_distance("ABC", "ABC") == 0
        assert hamming_distance("ABC", "ABD") == 1

        # 4. Test mapping
        out1_path = os.path.join(tmpdir, "output1.txt")
        mapping_results = map_reads(read_kmers, ref_data, ref_kmers, read_data, 3, out1_path)
        assert "read1" in mapping_results
        assert mapping_results["read1"][0][2] == 0  # mismatches = 0

        # 5. Test coverage calculation
        compute_reference_coverage(mapping_results, ref_data, read_data, tmpdir)
        out2_path = os.path.join(tmpdir, "output2.txt")
        assert os.path.exists(out2_path)

        print("All basic tests passed.")

if __name__ == "__main__":
    run_basic_tests()
