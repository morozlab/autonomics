'''

Author: Peter L. Williams
Version 1.1
2/27/2013

Module computes statistics for various assembly and
annotation results and stores the results in the
run_stats mysql table.  See the individual funtion
documentation below for details of specific stats
computed.

'''

from files.sequences import Utils
from multiprocessing import Process
from sqlalchemy.sql import select, and_, or_
import datetime
import hashlib
import imaplib
import os
import pickle
import redis
import shutil
import ssl
import subprocess
import sys
import threading
import time
from Bio import SeqIO
from autonomics import settings, utility, netutils

err_written = 0

def do_stats_update( rdict, proj_id, special_run, dir):
    '''
    rdict: a dictionary {statistics_name : statistic, ...} returned by the calculate
           functions below
    proj_id:  project_id of the project
    special_run: 0 if default; 1 if custom project configuration
    dir: project_name directory to write stats file to if special_run = 1

    This function takes a dictionary of statistics and enters them into the run_stats
    table if special_run = 0; otherwise writes the stats to a file named stats in dir.
    '''

    session = netutils.make_db_session()
    stats = netutils.get_table_object('run_stats', session)
    for key, value in rdict.items():
        if value == 'NULL':
            continue
        pair = ""
        pair = {str(key): utility.convert_if_int(value)}
        print "STATS: ", str(key),"  ", utility.convert_if_int(value)
        if not special_run:
            print "updating run_stats for project_id: ", proj_id
            u = stats.update().where(
                             stats.c.project_id==proj_id).values(pair).execute()
        else:
            file = dir + 'stats'
            cmd = 'echo ' + str(key) + ' ' + utility.convert_if_int(value) + ' >> ' + file
            os.system(cmd)
    session = None


def run_cmd(file, cmd, txt ):
    '''
    file: name of file to run write cmd on.
    cmd:  command to run (grep, cat, cut, paste, etc), chosen by caller.
    txt: name of statistic being worked on.

    This function takes a pipe cmd (like grep, cut, cat, etc) and reads the 
    output (res), returning a single dictionary {txt : res}
    '''
    global err_written
    if os.path.exists(file):
        for x in os.popen(cmd).readlines():
            res = x[:-1]
        return {txt:res}
    else:
        if not err_written:
            sys.stderr.write(file + " does not exist\n")
            err_written = 1
        return {txt:'NULL'}

def CalculatePfamStatistics(data_file):
    '''
    data_file: file on which stats are to be calculated

    This function counts the number of pfam hits and the number of 
    transcripts with pfam hits in the input data_file, returning a 
    dictionary of this information
    '''
    global err_written
    err_written = 0
    cmd = "grep PF " + data_file + " | wc -l"
    txt = "num_pfam_hits"
    dict1 = run_cmd(data_file, cmd, txt)
    cmd = 'grep PF ' + data_file + ' | cut -f1,1 -d " " | sort -u | wc -l';
    txt = "num_transcripts_with_pfam_hits"
    dict2 = run_cmd(data_file,cmd, txt)
    return dict(dict1.items() + dict2.items())

def CalculateBlastNRStatistics(data_file):
    '''
    data_file: file on which stats are to be calculated

    This function counts the number of blast_nr hits and the number of 
    transcripts with blast_nr hits in the input data_file, returning a 
    dictionary of this information
    '''
    global err_written
    err_written = 0
    cmd = 'grep "^>" ' + data_file + ' | wc -l'
    txt = "num_blast_nr_hits"
    dict1 = run_cmd(data_file, cmd, txt)
    cmd = 'grep "Sequences producing significant alignments:" ' + \
                                             data_file + ' | wc -l';
    txt = "num_transcripts_with_blast_nr_hits"
    dict2 = run_cmd(data_file,cmd, txt)
    return dict(dict1.items() + dict2.items())


def CalculateGOStatistics(data_file):
    '''
    data_file: file on which stats are to be calculated

    This function counts the number of GO hits and the number of 
    transcripts with GO hits in the input data_file, returning a 
    dictionary of this information
    '''
    global err_written
    err_written = 0
    cmd = 'cat ' + data_file + ' | wc -l'
    txt = "num_go_hits"
    dict1 = run_cmd(data_file,cmd, txt)
    cmd = 'cut -f1,1 -d " " ' + data_file + ' | sort -u | wc -l'
    txt = "num_transcripts_with_go_hits"
    dict2 = run_cmd(data_file,cmd, txt)
    return dict(dict1.items() + dict2.items())
    

def CalculateKEGGStatistics(data_file):
    '''
    data_file: file on which stats are to be calculated

    This function counts the number of KEGG hits and the number of 
    transcripts with KEGG hits in the input data_file, returning a
    dictionary of this information
    '''
    global err_written
    err_written = 0
    cmd = 'cat ' + data_file + ' | wc -l'
    txt = "num_kegg_hits"
    dict1 = run_cmd(data_file,cmd, txt)
    cmd = 'cut -f1,1 -d " " ' + data_file + ' | sort -u | wc -l'
    txt = "num_transcripts_with_kegg_hits"
    dict2 = run_cmd(data_file,cmd, txt)
    return dict(dict1.items() + dict2.items())


def CalculateBlastSwissprotStatistics(data_file):
    '''
    data_file: file on which stats are to be calculated

    This function counts the number of blast_swissprot hits and the number of 
    transcripts with blast_swissprot hits in the input data_file, returning a
    dictionary of this information
    '''
    global err_written
    err_written = 0
    cmd = 'grep "^>" ' + data_file + ' | wc -l';
    txt = "num_blast_swissprot_hits"
    dict1 = run_cmd(data_file,cmd, txt)
    cmd = 'grep "Sequences producing significant alignments:" ' + \
                                                      data_file + '| wc -l'
    txt = "num_transcripts_with_blast_swissprot_hits"
    dict2 = run_cmd(data_file,cmd, txt)
    return dict(dict1.items() + dict2.items())


def CalculateAssemblyStatistics(data_file):
    '''
    data_file: file on which stats are to be calculated

    This function counts the number of assembled bases, num_transcripts created,
    minimum, mean, maximum transcript lengths, and n50, n50_len, n90 and n90_len
    measures on the input data_file.  Returns a dictionary of this information.
    '''
    adict = dict()
    if os.path.exists(data_file):
        res = ""
        cmd = "python " + settings.SCRIPTPATH + "assemstats.py 0 " + \
                                                                  data_file
        for x in os.popen(cmd).readlines():
            res = x[:-1]
        res = res.split()
        adict["num_assembled_bases"] = res[1]
        adict["num_transcripts"] = res[2]
        # res[3] = trim_n
        adict["min_transcript_len"] = res[4]
        # res[5] = med
        adict["mean_transcript_len"] = res[6]
        adict["max_transcript_len"] = res[7]
        adict["n50"] = res[8]
        adict["n50_len"] = res[9]
        adict["n90"] = res[10]
        adict["n90_len"] = res[11]
        return adict
    else:
       sys.stderr.write(data_file + " does not exist\n")
       return {'assembly':'NULL'}


def CalculateReadsStatistics(data_file):
    '''
    data_file: file on which stats are to be calculated

    This function counts the number of raw reads and number of raw bases in
    the input data_file.  Returns a dictionary of this information.
    '''
    num = 0
    bases = 0
    if os.path.exists(data_file):
        for record in SeqIO.parse(data_file, 'fastq'):
            num+= 1
            bases += len(str(record.seq))
        rdict = dict()
        rdict['num_raw_reads'] = num
        rdict['num_raw_bases'] = bases
        return rdict
    else:
        sys.stderr.write(data_file  + " does not exist\n")
        return {'reads':'NULL'}


def CalculateQuantificationStatistics(data_file):
    '''
    data_file: file on which stats are to be calculated

    This function returns a dictionary with a flag
    indicating whether 'quantification_done' is 'Y' or 'N'.
    '''
    if os.path.exists(data_file):
        return {"quantification_done":'Y'}
    else:
        return {"quantification_done":'N'}
