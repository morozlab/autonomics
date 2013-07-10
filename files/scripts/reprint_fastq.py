import argparse
import files.Sequences

parser = argparse.ArgumentParser(description="This script reprints a FASTQ file, hopefully getting rid of weird errors while trying to use it")

parser.add_argument('inputFile')

args = parser.parse_args()

files.Sequences.SequenceUtils.reprintFastQFile(args.inputFile)

