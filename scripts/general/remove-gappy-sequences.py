from Bio import SeqIO
import sys

gap_frac = float(sys.argv[3])

# Read alignment file
alignment_file = sys.argv[1]
remove_seq = []
for record in SeqIO.parse(alignment_file, "fasta"):
    if record.seq.count("-")/len(record.seq) > gap_frac:
        print(f"Remove {record.id}")
        remove_seq.append(record.id)

# Read sequence file
sequence_file = sys.argv[2]
keep = []
for record in SeqIO.parse(sequence_file, "fasta"):
    if record.id not in remove_seq:
        keep.append(record)

SeqIO.write(keep, sys.argv[3], "fasta")
print(f"Removed {len(remove_seq)} sequences")