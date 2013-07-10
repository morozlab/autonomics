import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to format sequences for use with the ABrowse genome browser.")
parser.add_argument('-f', '--input-file', dest='file', required=True)
parser.add_argument('-t', '--format', dest='fileFormat', default='fasta')
parser.add_argument('-d', '--output-directory', dest='outputDir', required=True)

args = parser.parse_args()

seqs = SeqIO.parse(args.file, args.fileFormat)

chromList = open('chromosome_length_list.txt', 'w')
index = 1;
for record in seqs:
    outFile = open(args.outputDir + "/chr" + str(index) + ".con", 'w')
    outFile.write(str(record.seq) + "\n")
    outFile.flush()
    outFile.close()
    chromList.write("chr" + str(index) + "=" + str(len(record.seq)) + "\n")
    index += 1

chromList.close()
    