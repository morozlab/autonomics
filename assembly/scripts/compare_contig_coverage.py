import argparse

parser = argparse.ArgumentParser(description="Compares the coverage of two sets of sequences in a third set of sequences - basically a Venn Diagram of the coverage between the two aligned projects.")
parser.add_argument('-f1', '--first-coverage-file', dest='file1', required=True)
parser.add_argument('-f2', '--second-coverage-file', dest='file2', required=True)

args = parser.parse_args()

first = open(args.file1, 'r')
second = open(args.file2, 'r')

#generate a dictionary of the sequences that covered sequences from the first file
coveringSeqs = {}

line = first.readline()
firstTotal = 0
while line:
    columns = line.split("\t")
    ids = columns[2].split(";")
    for seqID in ids:
        if(seqID != ""):
            coveringSeqs[seqID] = 1
    line = first.readline();

first.close()
    
sharedMap = {}
secondTotal = 0


#determine which sequences also cover sequences from the second file
line = second.readline()
uniqueToSecond = {}
while line:
    columns = line.split("\t")
    ids = columns[2].split(";")
    for seqID in ids:
        if(seqID != ""):
            if(seqID in coveringSeqs):
                sharedMap[seqID] = 1
            else:
                uniqueToSecond[seqID] = 1

    line = second.readline()
    
second.close()

shared = len(sharedMap.keys())
firstTotal = len(coveringSeqs.keys())
secondTotal = len(uniqueToSecond.keys())

#print out the stats for the two ets of sequences
print(str(firstTotal - shared) + " sequences align to reads from only the first data set.")
print(str(secondTotal)  + " sequences align to reads from only the second data set.")
print(str(shared) + " sequences align to reads from both data sets.")

