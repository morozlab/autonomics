import argparse
import os

#collect arguments for calculating the contig coverage
parser = argparse.ArgumentParser(description='Calculates coverage for contigs in a given sequence set')
parser.add_argument('-a', '--alignment-file', dest='alignmentFile', required=True)
parser.add_argument('-f', '--format', dest='format', required=True)
parser.add_argument('-c', '--coverage-cutoff', dest='cutoff', required=True)
args = parser.parse_args()

coverage = {}
coverageSeqs = {}


def storeCoverage(query, start, end):
        for i in range(int(start), int(end)+1):
            coverage[query][i] = 1

#open the file
alignments = open(args.alignmentFile, 'r')
#check the format of the alignment file
if(args.format == "blast"):
    line = alignments.readline()
    currentQuery = ""
    while line:
        columns = line.split("\t")
        thisQuery = columns[0]
        #if this is a new query sequence, initiate the coverage for it in the dictionary
        if(thisQuery != currentQuery):
            currentQuery = thisQuery
            #initialize the coverage for the new sequence
            coverage[currentQuery] = []
            coverageSeqs[currentQuery] = ""
            for i in range(0, int(columns[12]) + 1):
                coverage[currentQuery].append(0)
        
        #if the coverage is above our threshold for either sequence, add the coverage to the total coverage
        if((float(columns[3])/float(columns[12]) >= float(args.cutoff)) or (float(columns[3])/float(columns[13]) >= float(args.cutoff))):
            storeCoverage(currentQuery, columns[6], columns[7])
            coverageSeqs[currentQuery] = coverageSeqs[currentQuery] + columns[1] + ";"
        line = alignments.readline()
        
#iterate over the dictionary, printing the coverage for each query and the sequences contributing to the coverage
for key in coverage.keys():
    seqLen = len(coverage[key]) - 1
    total = 0
    for i in range(len(coverage[key])):
        #print(i)
        #print(coverage[key][i])
        if(coverage[key][i] == 1):
            total += 1
    print(key + "\t" + str(float(total)/float(seqLen)) + "\t" + coverageSeqs[key])