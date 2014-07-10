import sys
from Bio import SeqIO

##############################################################################
#
#  Peter L. Williams
#  December 28, 2012
#
#  fastq2fasta.py
#
#    Call: fastq2fasta.py fastq.file fasta.file
#
#    converts fastq.file into fasta format
#
##############################################################################

argc = len(sys.argv)
if (argc != 3):
   print "Call: fastq2fasta.py fastq.file fasta.file"
   sys.exit()
infile = sys.argv[1]
outfile = sys.argv[2]

n = SeqIO.convert(infile, "fastq", outfile, "fasta")
print "Converted ", n, " seqs"


