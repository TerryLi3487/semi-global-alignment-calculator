import numpy as np
from cell import Cell
import sys
from Bio import SeqIO
from assn1_part1 import fit_align

query_seq = []
library_seqs = []

# Parse the input files 
query_rec = list(SeqIO.parse(sys.argv[1], "fasta"))
library_rec = list(SeqIO.parse(sys.argv[2], "fasta"))

query_seq = [r.seq for r in query_rec]
library_seqs = [r.seq for r in library_rec]


# array for scores 
scores = []
# the aligned sequences for the sequence pairs 
seq_aligns = []

for seq in library_seqs:
  score, seq1, seq2 = fit_align(str(query_seq[0]), seq)
  scores.append(score)
  seq_aligns.append((seq1, seq2))

# Find the best alignment score 
max_score = max(scores)
arr = np.array(scores)
max_index = arr.argmax()

# Write to output file 
f = open(sys.argv[3], "w")
f.write(">")
f.write(str(library_rec[max_index].id))
f.write("\n")
f.write(str(library_rec[max_index].seq))
f.write("\n")
f.write(str(max_score))
f.write("\n")
f.write(''.join(map(str, seq_aligns[max_index][0])))
f.write("\n")
f.write(''.join(map(str, seq_aligns[max_index][1])))
f.write("\n")
