import argparse
from file.SequenceUtils import SequenceUtils
 
parser = argparse.ArgumentParser(description='This script gets sequences from various file formats')
parser.add_argument('-f', '--seq-file', dest='seqFile', required=True)
parser.add_argument('-t', '--seq-file-type', dest='fileType', default='tab')
parser.add_argument('-c1', '--id-column', dest='column1', default=0)
parser.add_argument('-c2', '--seq-column', dest='column2', default=1)
parser.add_argument('-o', '--output-format', dest='outputFormat', default="fasta")
 
args = parser.parse_args()
 
if(args.fileType == 'MS'):
    SequenceUtils.printMassSpecSequences(args.seqFile, "tab", int(args.column2) - 1)
