import sys
import os
from Bio import SeqIO

# Usage: python split_fasta.py input.fasta n_splits output_folder
input_fasta = sys.argv[1]
output_folder = sys.argv[2]
n_parts = int(sys.argv[3])

records = list(SeqIO.parse(input_fasta, "fasta"))
total = len(records)

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Compute base size and remainder
base = total // n_parts
remainder = total % n_parts

index = 0
for i in range(n_parts):
    size = base + (1 if i < remainder else 0)
    chunk = records[index:index+size]
    index += size
    with open(f"{output_folder}/split_{i+1:03}.fasta", "w") as out:
        SeqIO.write(chunk, out, "fasta")

