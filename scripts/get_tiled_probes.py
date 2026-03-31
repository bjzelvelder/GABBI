#!/bin/python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
from collections import defaultdict

def generate_probes(seq, probe_size, step):
    length = len(seq)
    total_needed = (length - probe_size) // step + 1
    span = probe_size + (total_needed - 1) * step

    if span > length or total_needed < 1:
        return []

    start = (length - span) // 2
    probes = [seq[start + i * step : start + i * step + probe_size] for i in range(total_needed)]
    return probes

def split_probes(input_fasta, output_fasta, probe_size, overlap_size):
    records = []
    step = probe_size - overlap_size

    if step <= 0:
        raise ValueError("Overlap must be smaller than probe size")

    probe_counters = defaultdict(int)

    for record in SeqIO.parse(input_fasta, "fasta"):
        header = record.id
        match = re.search(r"(uce-\d+)\|([^|]+)", header)
        if not match:
            continue

        uce_id = match.group(1)          # uce-XXX
        species = match.group(2)         # Species

        probes = generate_probes(str(record.seq).upper(), probe_size, step)

        for probe in probes:
            probe_num = probe_counters[uce_id]
            probe_counters[uce_id] += 1
            probe_id = f"{uce_id}_p{probe_num + 1} |{species}"
            records.append(SeqRecord(Seq(probe), id=probe_id, description=""))

    with open(output_fasta, "w") as out:
        SeqIO.write(records, out, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Split FASTA sequences into overlapping probes centered on the sequence.")
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("output", help="Output FASTA file")
    parser.add_argument("--probe-size", type=int, default=120, help="Length of each probe (default: 120)")
    parser.add_argument("--overlap-size", type=int, default=80, help="Overlap between probes (default: 80)")

    args = parser.parse_args()
    split_probes(args.input, args.output, args.probe_size, args.overlap_size)

if __name__ == "__main__":
    main()

