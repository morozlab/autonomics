from Bio import SeqIO
from files.io import append_to_filename
from operator import methodcaller
import re
import os
import sys

class Utils:

    def __init__(self):
        '''
        Constructor
        '''

    @staticmethod
    def getDictFromFile(seqFile, seqFormat):
        returnDict = {}
        f = open(seqFile, 'Ur')
        for record in SeqIO.parse(f, seqFormat):
            returnDict[record.id] = record

        return returnDict

    @staticmethod
    def match(m,str1,str2): #Matches strings with m or less substitutions at the head
        for i,_ in enumerate(str1):
                if str1[i]!=str2[i]: m -= 1
                if m < 0: return False,i
        return True,(i+1)

    @staticmethod
    def printFastqFromDict(seqDict):
        for key, record in seqDict.iteritems():
            record = record.upper()
            sys.stdout.write(record.format('fastq-illumina'))

    @staticmethod
    def sort(sequences, sortType):
        if(sortType == 'length' or sortType == 'len'):
            return sorted(sequences, key=lambda record: len(str(record.seq)), reverse=True)

    @staticmethod
    def printMassSpecSequences(seqFile, fileFormat, index = None):
        if(fileFormat == "tab"):
            #assumes the sequence ID is in the first column, prints out forward and reverse sequences
            MSFile = open(seqFile, 'r')
            for line in MSFile:
                seqid = line.split("\t")[0]
                seq = line.split("\t")[index]
                seq = re.sub("\-{0,1}\(.+\)", '', seq)
                if(seq.count(".")):
                    seq = seq.split(".")[1]
                print(">" + seqid + " " + seq)
                print(seq)
                reverse = seq[::-1]
                print(">" + seqid + " " + reverse)
                print(reverse)

    @staticmethod
    def remove_duplicates(seq_file, seq_format):
        seq_ids = set()
        fh = open(seq_file, 'r')
        base, ext = os.path.splitext(seq_file)
        tmp = open(append_to_filename(seq_file, '_no_duplicates'), 'w')
        for record in SeqIO.parse(fh, seq_format):
            if(not record.id in seq_ids):
                seq_ids.add(record.id)
                tmp.write(record.format(seq_format))

        tmp.close()
        fh.close()

    @staticmethod
    def reprintFastQFile(seqFile):
        toRead = open(seqFile, 'r')
        state = "query"
        for line in toRead:
            #sys.stderr.write(line)
            if(line.startswith("\n")):
                continue
            line = line.strip()
            if(state == "query"):
                if(line.startswith("@")):
                    sys.stdout.write(line + "\n")
                    state = "seq"
            elif(state == "seq"):
                sys.stdout.write(line + "\n")
                state = "id"
            elif(state == "id"):
                if(line.startswith("+")):
                    sys.stdout.write(line + "\n")
                    state = "qual"
                else:
                    sys.stderr.write("Invalid FASTQ format, exiting...")
                    sys.exit()
            elif(state == "qual"):
                sys.stdout.write(line + "\n")
                state = "query"

    @staticmethod
    def splitSeqFile(queryFile, seqsPerFile):
        fh = open(queryFile, 'r')
        directory, file = os.path.split(queryFile)
        fileBase, fileExt = os.path.splitext(file)
        seqCounter = 0
        fileCounter = 1
        filenames = []
        if(directory == ''):
            curFile = fileBase + "_" + str(fileCounter) + fileExt
        else:
            curFile = directory + "/" + fileBase + "_" + str(fileCounter) + fileExt
        fw = open(curFile, 'w')
        filenames.append(curFile)
        for line in fh:
            line = line.strip()
            if(">" in line):
               seqCounter += 1
               if(seqCounter > int(seqsPerFile)):
                    #close the open partition file
                    fw.close()
                    #increment the file counter and open the next file
                    fileCounter += 1
                    if(directory == ''):
                        curFile = fileBase + "_" + str(fileCounter) + fileExt
                    else:
                        curFile = directory + "/" + fileBase + "_" + str(fileCounter) + fileExt
                    fw = open(curFile, 'w')
                    filenames.append(curFile)
                    #write this line to the new file
                    fw.write(line + "\n")
                    #reset the seqCounter to 1
                    seqCounter = 1

               else:
                    fw.write(line + "\n")

            else:
                fw.write(line + "\n")

        fw.close()
        return filenames

    @staticmethod
    def trimBarcode(seqs, barcode, m):
        for key, seq in seqs.iteritems():
            matched, index = Utils.match(1, seq, barcode)
            if(matched):
                ###FINISH IMPLEMENTING HERE, DAWG!
                pass


class Filter():

    seqs = {}
    verbose = False

    def __init__(self):
        '''
        Constructor!!
        '''

    def filterByLength(self, seqFile, minLen):
        pass

    def filterByQualityFasta(self, seqFile, qualFile, minQual, percent):
        #get a dictionary of sequences
        pass


    def filterFastqByQuality(self, seqFile, minQual, percent):
        returnDict = {}
        seqs = Utils.getDictFromFile(seqFile, "fastq")
        for key, seq in seqs.iteritems():
            belowThresh = 0
            for qual in seq.letter_annotations['phred_quality']:
                if(qual < minQual):
                    belowThresh += 1
            if(float(belowThresh)/float(len(seq.letter_annotations['phred_quality'])) > percent/100.0):
                returnDict[key] = seq
            elif(self.verbose):
                sys.stderr.write(seq.id + " fell below quality requirement of " + str(minQual) + "(" + str(belowThresh/len(seq.letter_annotations['phred_quality'])) + ") quality in at least " + str(percent) + " percent of bases.\n")

        return returnDict
