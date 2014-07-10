'''
Peter L. Williams
November 14, 2012 .. April 27, 2013

This script implements the steps needed to complete read normalization.
(read normalization is normally done after running adapter_trim and
quality_trim, which leave the data in fastq format.)

If the data is paired end, the data from the two ends is interleaved.
Then the input data file is converted from fastq to fasta format.  If the
number of sequences in the fasta file >= settings.KHMER_MIN_SEQS, khmer
is run. After khmer completes, if the data is paired-end, the data is
de-interleaved.

args:
  -in_file : fastq file with path to be operated on;
             expects file name is:  <sampleName>.fastq

  -paired_end : use only if data is paired end; 
                expects second file name is: <sampleName>.fastq.end2

     output: if single input file: <sampleName>.khmer.out
             if paired end data: 
                        <sampleName>.khmer.out.1 <sampleName>.khmer.out.2

     runs:    -- if paired end data, interleaves pair of cutadapted files.
              -- converts fastq files to fasta format
              -- counts seqs in fasta file, and if number seqs > 
                                                       $KHMER_MIN_SEQS, runs:
                       khmer to normalize fasta file
                       (normalize-by-median.py -C 30 -k 20 -N 1 -x 4.0e9)
              -- if paired end data, deinterleaves above result
                 (SequenceUtils.py -f outfile -ft fasta  --split-paired-end)
'''

import argparse
import os
import sys
import re
from autonomics import netutils, settings

def main():

    ##########################################################################
    #
    #    Get & check args.  Set paired_end = 0/1; get work_dir & script paths
    #
    ##########################################################################

    parser = argparse.ArgumentParser(description= \
                                            'Submission script for running .')
    parser.add_argument('-in_file', dest='in_file', help = 
                                        'fastq file with path', required=True)
    parser.add_argument('-paired_end', dest='paired_end', default = 0)
    args = parser.parse_args()

    paired_end = 1
    if args.paired_end == 0:
        paired_end = 0
    in_file = args.in_file

    work_dir = os.path.dirname(in_file) + "/"
    in_file_basename = os.path.basename(in_file)

    khmer_path = settings.khmer_path
    python_scripts_path = settings.SCRIPTPATH

    KHMER_MIN_SEQS = settings.KHMER_MIN_SEQS

    ##########################################################################
    #
    #    Deal with paired end case and set single_file
    #
    ##########################################################################

    in_file1 = in_file
    in_file2 = ""
    single_file = 1
    if paired_end:
        in_file2 = in_file + ".end2"
        single_file = 0
    print 'starting pre_assemble'

    ##########################################################################
    #
    #    Get sample name
    #
    ##########################################################################

    match_obj = re.search("(.*?)\.fastq",in_file_basename)
    sample_name = ""
    if match_obj:
        sample_name =  match_obj.group(1)
    else:
        print "pre_assemble.py unable to parse sample_name from input: " + \
              in_file_basename      
        sys.exit()

    ##########################################################################
    #
    #    confirm all files/dirs actually exist before starting
    #
    ##########################################################################
    
    if not os.path.isfile(in_file1):
        print in_file1 + ' does not exist'
        sys.exit()
    if not single_file and not os.path.isfile(in_file2):
        print in_file2 + ' does not exist'
        sys.exit()

    ##########################################################################
    #
    #    Change to work dir
    #
    ##########################################################################

    os.chdir(work_dir)    

    ##########################################################################
    #
    #     If paired end data, interleave the pair of cutadapted fastq files
    #        interleaved_file = <sample_name>.fastq_interleaved.fastq
    #     else
    #        interleaved_file = <sample_name>.fastq
    #
    ##########################################################################

    interleaved_file = ""
    if not single_file:
	interleaved_file = in_file1 + "_interleaved.fastq"
	f1 = in_file1
	f2 = in_file2
	icmd = python_scripts_path + "interleave_fastq.py"
	cmd = 'python ' + icmd + ' -f1 ' + f1 + ' -f2 ' + f2 + ' -o ' + \
                                          interleaved_file + '  > /dev/null'
	print cmd
        sys.stdout.flush()
        res =os.system(cmd)
        if res:
            print 'os.system(' + cmd + ') returns: ' + str(res)
            print cmd + " FAILED"
            sys.exit()
    else:
        interleaved_file = in_file1

    ##########################################################################
    #
    #    Convert $interleaved_file to fasta format
    #       ==> $fasta_file = <sample_name>.fasta
    #         OR
    #       ==> $fasta_file = <sample_name>.fastq_interleaved.fasta
    #
    ##########################################################################

    py_cmd = python_scripts_path + "fastq2fasta.py"
    match_obj = re.search("(.*?)\.fastq",interleaved_file)
    fasta_file = ""
    if match_obj:
        fasta_file =  match_obj.group(1) + ".fasta"
    else:
        print "pre_assemble.py unable to parse fasta_file from \
               interleaved_file: " + interleaved_file      
        sys.exit()
    cmd = "python " + py_cmd + ' ' + interleaved_file + ' ' + fasta_file
    print cmd
    sys.stdout.flush()
    res = os.system(cmd)
    if res:
        print 'os.system(' + cmd + ') returns: ' + str(res)
        print cmd + " FAILED"
        sys.exit()

    ##########################################################################
    #
    #    Count seqs in fasta_file; if number seqs > settings.KHMER_MIN_SEQS, run khmer.
    #       ==> $fasta_file = <sample_name>.fasta
    #         OR
    #       ==> $fasta_file = <sample_name>.fastq_interleaved.fasta
    #     ==> <sample_name>.khmer.out
    #
    ##########################################################################

    seq_count = 0
    if not os.path.isfile(fasta_file):
        print fasta_file + ' does not exist'
        sys.exit()
    ff = open(fasta_file, 'r')
    for line in ff:
        match_obj = re.search("^>.*",line)
        if match_obj: seq_count = seq_count + 1
    print "Found " + str(seq_count) + ' seqs in ' + fasta_file

    khmer_out = work_dir + sample_name + ".khmer.out"
    if seq_count >= settings.KHMER_MIN_SEQS:
        khmer_out2 = fasta_file + ".keep"  # what khmer calls its output.
	khmer_cmd = khmer_path + "/normalize-by-median.py"

	cmd = "python " + khmer_cmd + ' -C 30 -k 20 -N 5 -x 4.0e9 ' + \
						    fasta_file + ' > /dev/null'
	print cmd
        sys.stdout.flush()
	res = os.system(cmd)
        if res:
            print 'os.system(' + cmd + ') returns: ' + str(res)
            print cmd + " FAILED"
            sys.exit()

	cmd = "mv " + khmer_out2 + ' ' + khmer_out
		  # sample_name.fasta.keep ==> sample_name.khmer.out
	print cmd
        sys.stdout.flush()
	res = os.system(cmd)
        if res:
            print 'os.system(' + cmd + ') returns: ' + str(res)
            print cmd + " FAILED"
            sys.exit()

    else:
        cmd = "mv " + fasta_file + ' ' + khmer_out
                       # sample_name.fasta ==> sample_name.khmer.out
        print cmd
        sys.stdout.flush()
	res = os.system(cmd)
        if res:
            print 'os.system(' + cmd + ') returns: ' + str(res)
            print cmd + " FAILED"
            sys.exit()

        print "Only " + str(seq_count) + " seqs so skipping khmer"
        sys.stdout.flush()

    ##########################################################################
    #
    #    Deinterleave data
    #
    #     expects input in khmer_out = <sample_name>.khmer.out
    #
    #     output = <sample_name.fasta> & <sample_name>.end2.fasta
    #
    #     else: = <sample_name.fasta>
    #
    ##########################################################################

    deinterleaved1 = khmer_out + ".1"
    deinterleaved2 = khmer_out + ".2"

    if not single_file:
	dcmd = python_scripts_path + "seqtools.py"
	cmd = "python " + dcmd + ' -f ' + khmer_out + \
                                  ' -ft fasta --split-paired-end > /dev/null'
	print cmd
        sys.stdout.flush()
        res = os.system(cmd)
        if res:
            print 'os.system(' + cmd + ') returns: ' + str(res)
            print cmd + " FAILED"
            sys.exit()

	out1 = sample_name + ".fasta"
	out2 = sample_name + ".fasta.end2"
	cmd = "mv " + deinterleaved1 + ' ' + out1
	print cmd
        sys.stdout.flush()
        res = os.system(cmd)
        if res:
            print 'os.system(' + cmd + ') returns: ' + str(res)
            print cmd + " FAILED"
            sys.exit()

	cmd = "mv " + deinterleaved2 + ' ' +out2
	print cmd
        sys.stdout.flush()
        res = os.system(cmd)
        if res:
            print 'os.system(' + cmd + ') returns: ' + str(res)
            print cmd + " FAILED"
            sys.exit()

    else:
	# dont need to de-interleave, rename khmer output to
        # sample_name . ".fasta"
	out = work_dir + sample_name + ".fasta"
	cmd = "mv " + khmer_out + ' ' + out
	print cmd
        sys.stdout.flush()
        res = os.system(cmd)
        if res:
            print 'os.system(' + cmd + ') returns: ' + str(res)
            print cmd + " FAILED"
            sys.exit()
    print "pre_assemble complete"
   
if __name__ == '__main__':
    main()
