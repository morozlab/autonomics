import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="This script changes the sequence ID on genomic scaffolds based on the set command-line flag")

parser.add_argument("-p1", "--first-period", dest="firstPeriod", default=False, const=True, action='store_const')
parser.add_argument("sequenceFile")
parser.add_argument("sequenceFormat")

args = parser.parse_args()

seqs = SeqIO.parse(args.sequenceFile, args.sequenceFormat)

for seq in seqs:
    seqID = seq.id
    seq = str(seq.seq)
    if(args.firstPeriod and seqID.find(".")):
        seqID = seqID.split(".")[0]
    if(args.sequenceFormat == "fasta"):
        print(">" + seqID)
        print(seq)
    
    
        
    
        



