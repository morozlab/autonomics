import argparse
from Bio import SeqIO
from Sequences.SequenceUtils import SequenceUtils



def main():
    parser = argparse.ArgumentParser(description="This script sorts sequence files!")
    parser.add_argument('-f', '--seq-file', dest='toSort', help='The sequence file you wish to sort.', required=True)
    parser.add_argument('-s', '--sort-type', dest='sortType', help='The type of sort you want to perform', required=True)
    parser.add_argument('-t', '--seq-type', dest='sequenceType', help='The format of the sequence file', requred=True)
    
    args = parser.parse_args()
    
    sortedSeqs =  SequenceUtils.sort(SeqIO.parse(args.toSort, args.sequenceType), args.sortType)
    
    for record in sortedSeqs:
        if(args.sequenceType == 'fastq'):
            print("@" + record.id + "\n")
            print(record.seq + "\n")
            print("+" + record.id + "\n")
            print(record.letter_annotations['phred_quality'])
        elif(args.sequenceType == 'fasta'):
            print(">" + record.id + "\n")
            print(record.seq + "\n")
            

if __name__ == '__main__':
    main()