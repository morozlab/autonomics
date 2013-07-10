#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Mathew Citarella
#
# Created:     01/10/2012
# Copyright:   (c) Mathew Citarella 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from Bio import SeqIO
import argparse
import subprocess


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-q", dest="queryFile", required=True)
    parser.add_argument("-t", dest="targetFile", required=True)
    parser.add_argument("--query-type", dest="queryType", default="dna")
    parser.add_argument("--target-type", dest="targetType", default="dna")
    parser.add_arugment("--query-format", dest="queryFormat", default="fasta")
    parser.add_argument("--cpu", dest="cpu", default="1")
    parser.add_argument("--identity", dest="minIdentity", default="90")
    parser.add_argument("--score", dest="minScore", default="30")
    parser.add_argument("--no-head", dest="noHead", default=False,const=True, action="store_const")

    #parse arguments
    args = parser.parse_args()

    #get a dictionary of the query sequences
    q = open(args.queryFile, 'r')
    qDict = SeqIO.to_dict(SeqIO.parse(q, args.queryFormat))

    #determine how many sequences go in each query file
    numPerFile = int(len(qDict))/int(args.cpu) + 1

    #split the query file and start the individual BLAT jobs
    for seqid, seq in qDict.iteritems():


if __name__ == '__main__':
    main()
