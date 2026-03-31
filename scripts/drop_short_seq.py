#!/usr/bin/python
import sys 
from Bio import SeqIO

# Usage : python drop_short_seq.py input output threshold
# The threshold corresponds to the minimum % of non-ambiguous characters required to keep a sequence, others are dropped

FastaFile = open(sys.argv[1], 'r')
FastaDroppedFile = open(sys.argv[2], 'w')
drop_cutoff = float(sys.argv[3])

if (drop_cutoff < 0) or (drop_cutoff > 100):
    print('\n Sequence drop cutoff must be in 0-100 range\n')
    sys.exit(1)

for seqs in SeqIO.parse(FastaFile, 'fasta'):
    name = seqs.id
    seq = seqs.seq
    seqLen = len(seqs)
    gap_count = sum(1 for base in seq if base in ['-', 'N', 'n', '?'])

    if seqLen == 0 or 100-(gap_count/float(seqLen))*100 <= drop_cutoff:
        print('dropped %s' % name)
    else:
        # Write as one-liner fasta
        FastaDroppedFile.write(f'>{name}\n{seq}\n')

FastaFile.close()
FastaDroppedFile.close()
