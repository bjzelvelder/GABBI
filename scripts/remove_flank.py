#!/bin/python

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

def find_core_alignment(seq, min_block_size=10):
    """
    Identify the core region of a sequence based on blocks of real bases (> min_block_size)
    separated by gaps.
    """
    blocks = []
    for match in re.finditer(r'[^-]+', seq):
        block_start, block_end = match.span()
        block_length = block_end - block_start
        blocks.append((block_start, block_end, block_length))

    # Filter blocks based on min_block_size
    valid_blocks = [(start, end) for start, end, length in blocks if length >= min_block_size]
    
    if not valid_blocks:
        raise ValueError("No valid block found with sufficient length.")

    true_start = valid_blocks[0][0]
    true_end = valid_blocks[-1][1]
    
    return true_start, true_end

def trim_alignment(input_fasta, ref_name, output_fasta, min_block_size=10):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    
    ref_record = next((r for r in records if r.id == ref_name), None)
    if not ref_record:
        print(f"Reference '{ref_name}' not found.")
        sys.exit(1)

    ref_seq = str(ref_record.seq)
    
    # Find the core region
    try:
        start, end = find_core_alignment(ref_seq, min_block_size=min_block_size)
    except ValueError as e:
        print(e)
        sys.exit(1)

    # Trim all sequences based on the ref start/end
    trimmed_records = []
    for r in records:
        trimmed_seq = r.seq[start:end]
        trimmed_records.append(SeqRecord(trimmed_seq, id=r.id, description=""))

    # Write the trimmed alignment
    with open(output_fasta, "w") as out_f:
        SeqIO.write(trimmed_records, out_f, "fasta")

    print(f"Trimmed alignment written to {output_fasta} (from {start} to {end})")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <alignment.fasta> <reference_name> <output.fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    ref_name = sys.argv[2]
    output_fasta = sys.argv[3]

    trim_alignment(input_fasta, ref_name, output_fasta)

