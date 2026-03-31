#!/bin/python3

import os
import sys

def parse_alignment_file(alignment_file, chromosomes):
    alignments = {}
    current_alignment = []

    with open(alignment_file, 'r') as f:
        for line in f:
            if line.startswith('a'):
                if current_alignment:
                    for chrom_name in chromosomes:
                        if any(chrom_name in seq_line for seq_line in current_alignment):
                            if chrom_name not in alignments:
                                alignments[chrom_name] = []
                            alignments[chrom_name].append(current_alignment)
                            break
                current_alignment = [line]
            else:
                current_alignment.append(line)
                
        # Ensure the last alignment is also processed
        if current_alignment:
            for chrom_name in chromosomes:
                if any(chrom_name in seq_line for seq_line in current_alignment):
                    if chrom_name not in alignments:
                        alignments[chrom_name] = []
                    alignments[chrom_name].append(current_alignment)
                    break
    
    return alignments

def write_chromosome_files(alignments, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    for chrom_name, alignment_blocks in alignments.items():
        with open(os.path.join(output_dir, f'{chrom_name}.maf'), 'w') as out_file:
            for block in alignment_blocks:
                for line in block:
                    out_file.write(line)

def main(alignment_file, bed_file, output_dir):
    with open(bed_file, 'r') as f:
        chromosomes = sorted(set(line.strip().split('\t')[0] for line in f))
    
    alignments = parse_alignment_file(alignment_file, chromosomes)
    write_chromosome_files(alignments, output_dir)
    print("Extraction finished.")

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} alignment_file bed_file output_dir")
        sys.exit(1)
    
    alignment_file = sys.argv[1]
    bed_file = sys.argv[2]
    output_dir = sys.argv[3]

    main(alignment_file, bed_file, output_dir)

