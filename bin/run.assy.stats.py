#!/usr/local/bin/python

'''

Author: Peter L. Williams
Version 1.1
2/27/2013

'''

from files.sequences import Utils
from multiprocessing import Process
import datetime
import os
import shutil
import ssl
import sys
import time
from Bio import SeqIO

def main():

    num_args = len(sys.argv)
    arg_list = sys.argv
    if num_args != 2:
       print "call is: ", arg_list[0], " project.fasta file name"
       sys.exit()

    data_file = arg_list[1]

    if os.path.exists(data_file):
        cmd = "python /home/pwilliams/autonomics/scripts/assemstats.py 0 " + data_file
        for x in os.popen(cmd).readlines():
            res = x[:-1]
        res = res.split()
#        print "num_assembled_bases = ", res[1]
#        print "num_transcripts = ", res[2]
        # res[3] = trim_n
#        print "min_transcript_len = ", res[4]
        # res[5] = med
#        print "mean_transcript_len = ", res[6]
#        print "max_transcript_len = ", res[7]
        print "n50: ", res[8]
        print "n50_len: ", res[9]
        print "n90: ", res[10]
        print "n90_len: ", res[11]
    else:
       sys.stderr.write(data_file + " does not exist\n")


if __name__ == '__main__':
    main()
