
"""
Script to split HZS-BC gene into B and C subunit.
"""
from Bio import SeqIO
import copy
import sys

for record in SeqIO.parse(sys.argv[1], "fasta"):
    if record.id == "GAX62881.1":
        hzs_b = copy.deepcopy(record)
        hzs_c = copy.deepcopy(record)

# Cut HZS-BC so that it fits with
hzs_b.seq = hzs_b.seq[:375]
SeqIO.write(hzs_b, sys.argv[2], "fasta")
hzs_c.seq = hzs_c.seq[375:]
SeqIO.write(hzs_c, sys.argv[3], "fasta")
