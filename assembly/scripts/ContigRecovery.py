import argparse
import os
import re
import subprocess
INSTALL_DIR = "/home/girardo/tuber/moroz-python/"
SCRIPT_ASSEM = "assembly/scripts/"
SCRIPTPATH = INSTALL_DIR + SCRIPT_ASSEM #Permission problems... cloned everything. Put this in settings!


parser = argparse.ArgumentParser(description="Helper script for recovering contig sequences from a finished assembly")
parser.add_argument("singletList", help="File containing a list of identifiers considered singlets from the final assembly")
parser.add_argument("singletFasta", help="FASTA file of singlets from the second assembly")
args = parser.parse_args()

def findContigs(file, outfile, pattern):
    line = file.readline()
    while line:
        if(re.search(pattern, line)):
            line = line.replace(">", "", 1)
            outfile.write(line)
        line = file.readline()

splitPath = os.path.split(args.singletFasta)
outfile = open(splitPath[0] + "/recoveredContigs.list", 'w')
singletList = open(args.singletList, "r")
findContigs(singletList, outfile, "contig")
singletList.close()
singletList = open(args.singletList, "r")
findContigs(singletList, outfile, "c\d+")
singletList.close()
outfile.close()

process = subprocess.Popen("perl " + SCRIPTPATH + "getSequencesByID.pl " + args.singletFasta + " " + splitPath[0] + "/recoveredContigs.list" + " fasta 50 > " + splitPath[0] + "/recoveredContigs.fa", shell=True)
process.wait()
