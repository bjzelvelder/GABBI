#!/bin/python

import re
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_header(header):
    match = re.search(r'>(slice_[0-9]+)\s*\|(\S+):(\d+)(?:\.0)?-(\d+)(?:\.0)?', header)
    if not match:
        raise ValueError(f"Header not recognized: {header}")
    slice_name, contig, start, end = match.groups()
    return slice_name, contig, int(start), int(end)

def merge_sequences(seqs):
    if not seqs:
        return []

    merged = []
    seqs.sort(key=lambda x: x[1])  # Sort by start

    current_slice, current_start, current_end, current_seq = seqs[0]

    for slice_name, start, end, seq in seqs[1:]:
        if start <= current_end:  # overlap
            overlap = current_end - start
            if overlap < len(seq):
                current_seq += seq[overlap:]
            current_end = max(current_end, end)
        else:
            merged.append((current_slice, current_start, current_end, current_seq))
            current_slice, current_start, current_end, current_seq = slice_name, start, end, seq

    merged.append((current_slice, current_start, current_end, current_seq))
    return merged

def process_fasta(input_file, output_file):
    contig_dict = defaultdict(list)

    with open(input_file) as f:
        records = list(SeqIO.parse(f, "fasta"))

    for record in records:
        slice_name, contig, start, end = parse_header(">" + record.description)
        contig_dict[contig].append((slice_name, start, end, str(record.seq)))

    merged_records = []
    for contig, seqs in contig_dict.items():
        merged = merge_sequences(seqs)
        for slice_name, start, end, sequence in merged:
            header_id = slice_name
            header_desc = f"{contig}:{start}.0-{end}.0"
            merged_records.append(SeqRecord(Seq(sequence), id=header_id, description=header_desc))

    SeqIO.write(merged_records, output_file, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage : python merge_overlapping_seq_in_fasta.py input.fasta output.fasta")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]

    process_fasta(input_fasta, output_fasta)

