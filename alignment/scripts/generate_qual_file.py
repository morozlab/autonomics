import argparse
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser(description="This script generates a default 454-style quality file for a FASTA file")
parser.add_argument("fastaFile", help="FASTA file for which we are creating default quality scores.")

args = parser.parse_args()

for seq_record in SeqIO.parse(args.fastaFile, "fasta"):
    sys.stdout.write(">" + seq_record.id + "\n")
    count = 0
    for c in seq_record.seq:
        sys.stdout.write("40 ")
        count += 1
        if((count % 20) == 0):
            sys.stdout.write("\n")
    
    sys.stdout.write("\n")

