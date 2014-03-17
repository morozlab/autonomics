'''

  Peter L. Williams
  November 14, 2012,  December 30, 2012,    January 10, 2013

  run_trinity.py

     input: -- 1 or 2 fast[a|q] files, WITH PATH .. if 1 file, assumes name
               is <sample_name>.fast[a|q], if 2 files: assumes names are
               <sample_name>.fast[a|q] & <sample_name>.end2.fast[a|q]
            -- number cpus to be used by trinity

     Example call: 
          python run_trinity.py --cpus NN  --in_file1
                      /srv/data2/pipeline/Pleurobrachia/Pleurobrachia.fasta
              OR
          python run_trinity.py  --cpus NN  --in_file1
                      /srv/data2/pipeline/Pleurobrachia/Pleurobrachia.fasta
                      -in_file2  
                   /srv/data2/pipeline/Pleurobrachia/Pleurobrachia.end2.fasta

     output: assembled data named: <sample name>_project.fasta

     runs:
        -- fastq2fasta if suffix is fastq
        -- trinity assembler:
                 Trinity.pl  --seqType fa --left xxx.khmer.out.1
                         --right xxx.khmer.out.2 --JM 100G
                         --inchworm_cpu <cpus> --CPU <cpus>
               -- removes intermediate files
'''

import argparse
import os
import re
from autonomics import netutils, settings

def main():
    parser = argparse.ArgumentParser(description='Submission script for running Trinity.')
    parser.add_argument('-in_file1', dest='in_file1', required=True)
    parser.add_argument('-in_file2', dest='in_file2', default = "")
    parser.add_argument('-cpus', dest='cpus', required=True)
    args = parser.parse_args()

    paired_end = 1
    if args.in_file2 == "":
        paired_end = 0
    in_file1 = args.in_file1
    cpus = args.cpus
    work_dir = os.path.dirname(in_file1) + "/"

    ###########################################################################
    #
    #    Check extension of in_files
    #    Get sample name & path to trinity
    #
    ###########################################################################

    in_file1_basename = os.path.basename(in_file1)
    match_obj = re.search("(.*?)\.(fast.*)",in_file1_basename)
    ext = ""
    sample_name = ""
    if match_obj:
        sample_name =  match_obj.group(1)
        if match_obj.group(2) == 'fastq': ext = 'fastq'
        elif match_obj.group(2) == 'fasta': ext = 'fasta'
        else:
            print "ERROR: bad extension: " + ext + " on " + in_file + " to run_trinity"
            sys.exit()
    else:
        print "run_trinity.py unable to parse sample_name from input: " + \
              in_file1_basename      
        sys.exit()
 
    trinity_path = settings.trinity_path

    ###########################################################################
    #
    #    change to work dir
    # 
    ###########################################################################

    os.chdir(work_dir)    

    ###########################################################################
    #
    #    Convert file to fasta if necessary
    #
    ###########################################################################

    cmd = settings.SCRIPTPATH + 'fastq2fasta.py'
    in_file1_fasta = work_dir + sample_name + ".fasta"
    in_file2_fasta = work_dir + sample_name + ".fasta.end2"

    if ext == 'fastq':
       cmd = "python " + cmd + ' ' + in_file1 + ' ' + in_file1_fasta
       print cmd
       os.system(cmd)
       if paired_end:
           cmd = "python " + cmd + ' ' + in_file2 + ' ' + in_file2_fasta
           print cmd
           os.system(cmd)

    ###########################################################################
    #
    #    Run trinity assembler
    #
    ###########################################################################

    trinity_out_dir = work_dir + 'trinity_out'
    trinity_out = trinity_out_dir + "/Trinity.fasta"
    tcmd = trinity_path + "/Trinity.pl"

#    mem = 40G
    mem = '256G'

    if paired_end:
        cmd = "perl " + tcmd + " --seqType fa --left " + in_file1_fasta + \
              ' --right ' + in_file2_fasta + ' --JM ' + mem + ' --inchworm_cpu ' + \
              cpus + " --CPU " + cpus + " --output " + trinity_out_dir + \
              " > /dev/null"
        print cmd
        os.system(cmd)
    else:
        cmd = "perl " + tcmd + " --seqType fa --single " + in_file1_fasta + \
              ' --JM ' + mem + ' --inchworm_cpu ' + cpus + " --CPU " + cpus + \
              " --output " + trinity_out_dir + " > /dev/null"
        print cmd
        os.system(cmd)

    ###########################################################################
    #
    #    Move assembled data to dir above
    #
    ###########################################################################

#   work_dir = <pipeline_path>/<proj_dir>
#    trinity_out_dir = work_dir + 'trinity_out'  == /srv/data2/pipeline/<proj>/trinity_out
#    trinity_out = trinity_out_dir + "/Trinity.fasta" == /srv/data2/pipeline/<proj>/trinity_out/Trinity.fasta

#    tgt = trinity_out_dir + '/' + sample_name + "_project.fasta"
    tgt = work_dir + sample_name + "_project.fasta"

    cmd = "mv " + trinity_out + " " + tgt
    print cmd
    os.system(cmd)

    src = trinity_out_dir + '/Trinity.timing'
    cmd = "mv " + src + " " + work_dir
    print cmd
    os.system(cmd)

    ##########################################################################
    #
    #    Remove intermediate files
    #
    ##########################################################################

    cmd = "rm -fr " + trinity_out_dir
    print cmd
    os.system(cmd)

if __name__ == '__main__':
    main()
