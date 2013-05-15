'''
Peter L. Williams
January 3, 2013 .. April 25, 2013

    args: --job_type <job_type> (required, e.g. adapter_trim, quality_trim,
                                 read_normalization, assemble, blast_nr,
                                 blast_swissprot, panther, go, kegg, 
                                 pfam, quantification) 

          --job_name <job_name> (required; unique name for job)
          --input_file <input_file_with_path> (required for all job types except
                                                                       quantification)
          --go_quant_file', (required only for GO jobs, use either the quantification
                              to use or 'OMIT' to not use any file.
          --args <nondefault arg(s) to use. (optional; args to use instead of those 
                                                                in default_args table)
            Format: --args "<arg_name>;<arg value(s)>" "<arg_name>;<arg value(s)>" ..]
            Examples:
              --args "pipeline_args; --aligner bowtie --db-type NT "
                     "process_args;-query <input> -evalue 10" 
                     "resources;cpu:100,mem:250mb"


This script inserts a single custom user job (e.g. blast_nr, assemble, etc.) into
the zeroclick pipeline.

This code: 
   reads, parses and verifies the args
   verifies input_file exists
   verifies job_type is valid
   checks that job_name is unique
   enters job in pn_mapping and jn_mapping tables.
   reads in default_args for this job type
   uses default_args except where --args has a value
   inserts args in the args table for this job_id
   inserts job into special_queue
'''

from zeroclick import settings, netutils
from zeroclick.file_io import FileExtensions
import os
from sqlalchemy import *
import argparse
import datetime
import shutil
import time
import subprocess
import sys
from os.path import isfile, join, isdir, exists
from ftplib import FTP
import argparse
from files.io import TabFile, convert_if_int, convert_if_number
import operator

def main():

    ###########################################################################
    #
    #   deal with sql database
    #
    ###########################################################################

    session = netutils.DBSession("localhost", settings.ZC_DB_NAME, 
                                settings.db_cred.user, settings.db_cred.passwd)
    default_args_table = netutils.getTableObject("default_args", session)
    args_table = netutils.getTableObject("args", session)
    queue_table = netutils.getTableObject(settings.special_queue, session)
    pn_mapping_table = netutils.getTableObject("pn_mapping", session)
    jn_mapping_table = netutils.getTableObject("jn_mapping", session)

    ###########################################################################
    #
    #   get all valid job_types from default_args table
    #
    # mysql> select job_type from default_args;
    # +--------------------+
    # | job_type           |
    # +--------------------+
    # | adapter_trim       |
    # | assemble           |
    # | blast_nr           |
    # | blast_swissprot    |
    # | go                 |
    # | kegg               |
    # | pfam               |
    # | quality_trim       |
    # | quantification     |
    # | read_normalization |
    # | upload             |
    # +--------------------+
    #
    ###########################################################################

    all_job_types = []
    job_type_str = ""
    results = session.conn.execute("select job_type from default_args")
    rows = results.fetchall()
    for row in rows:
        jt = row.job_type
        if jt == 'upload':
            continue
        all_job_types.append(jt)
        job_type_str = job_type_str + " " + jt
    all_job_types.append('panther')
    job_type_str = job_type_str + " panther"

    ###########################################################################
    #
    #   parse args
    #
    ###########################################################################

    print " "
    parser = argparse.ArgumentParser() 
    parser.add_argument('--job_type', help=job_type_str,required=True)
    parser.add_argument('--job_name', help='unique name for job',required=True)
    parser.add_argument('--input_file', help='file with full path to use as \
          input to job_type (not used for quantification jobs)',required=False)
    parser.add_argument('--go_quant_file', help='ONLY required for \
      job_type=go; quantification file with full path or enter the word OMIT',\
      required=False)
    parser.add_argument('--args', nargs='+', help='--args "pipeline_args; \
             --aligner bowtie --db-type NT --db /path ../XXX_project.fasta \
             --query /path ../YYY.fastq"')
    args = parser.parse_args()
    input_file = args.input_file
    go_quant_file = args.go_quant_file
    job_type = args.job_type
    job_name = args.job_name
    all_args = args.args
    print "job_type: ", job_type
    print "job_name: ", job_name
    print "input_file: ", input_file
    if ((job_type=='go') and (go_quant_file !=None)):  
        print "go_quant_file: ", go_quant_file
    print "args:", args.args

    len_all_args = 0
    if (all_args):
        len_all_args = len(all_args)
    running_mira = 0
    bowtie = 1
    if (len_all_args) :
        for arg in all_args:
            if 'mira' in arg: running_mira = 1

    ###########################################################################
    #
    #   get all valid column names in default_args table
    #
    ###########################################################################
    all_arg_names=[]
    results = session.conn.execute('select COLUMN_NAME from \
                   information_schema.columns WHERE TABLE_NAME="default_args"')
    rows = results.fetchall()
    for row in rows:
        cn = row.COLUMN_NAME
        if cn not in all_arg_names and cn != 'job_type':
            all_arg_names.append(cn)

    ###########################################################################
    #
    #   get valid input file extensions
    #
    ###########################################################################

    inputExtensions = {"panther": ".fasta", "quality_trim": ".fastq", 
        "adapter_trim": ".fastq", "read_normalization": ".fastq", 
        "blast": ".fasta", "blast_swissprot": ".fasta",
        "blast_nr": ".fasta", "blat" : ".fasta", "pfam": ".fasta",
        "go": "_blast_swissprot.txt", "kegg": "_blast_swissprot.txt"}

    outputExtensions = FileExtensions.outputExtensions

    ###########################################################################
    #
    #   initialize dictionaries for args
    #
    ###########################################################################

    cmd_line_args = {}
    default_args = {}
    args_to_use = {}

    ###########################################################################
    #
    #   enter cmd line args into cmd_line_args dictionary & find quant_db
    #   file, quant_query_file & quant_aligner quant_db_type for quantification
    #
    ###########################################################################

    found_pipeline_args = 0
    found_process_args = 0
    found_quant_db =0
    found_quant_aligner =0
    found_quant_db_type =0
    found_quant_query = 0
    quant_db_file = ""
    quant_query_file = ""
    quant_aligner = ""
    quant_db_type = ""

    if (len_all_args) :
        for arg in all_args:
            arg_list = arg.split(';')
            len_arg_list = len(arg_list)
            if len_arg_list < 2:
                 print "\nERROR: incomplete arg list to ", arg_list[0]
                 sys.exit()
            arg_name = arg_list[0]
            if arg_name not in all_arg_names:
                print "\nERROR: invalid arg_name: ", arg_name
                sys.exit()
            arg_value = arg_list[1]
            cmd_line_args[arg_name]=arg_value
            if arg_list[0] == 'pipeline_args':
                found_pipeline_args = 1
                if 'blast' in arg_list[1]:
                    bowtie = 0
                words = arg_list[1].split()
                len_words = len(words)
                if '--db' in words:
                  found_quant_db = 1
                  index = words.index('--db')
                  if index+1 < len_words:
                      quant_db_file = words[index+1]
                      if (not isfile(quant_db_file)):
                          print "\nERROR: --db <file> does not exist: ",quant_db_file
                          sys.exit()
                  else:
                      print "\nERROR: <file> in --db <file> is missing in pipeline_args"
                      sys.exit()
                if '--query' in words:
                    found_quant_query = 1 
                    index = words.index('--query')
                    if index+1 < len_words:
                        quant_query_file = words[index+1]
                        if (not isfile(quant_query_file)):
                            print "\nERROR: --query <file> does not exist: ",quant_query_file
                            sys.exit()
                    else:
                        print "\nERROR: <file> in --query <file> is missing in pipeline_args"
                        sys.exit()
                if '--aligner' in words:
                    found_quant_aligner = 1 
                    index = words.index('--aligner')
                    if index+1 < len_words:
                        quant_aligner = words[index+1]
                    else:
                        print "\nERROR: <arg> in --aligner <arg> is missing in pipeline_args"
                        sys.exit()
                if '--db-type' in words:
                    found_quant_db_type = 1 
                    index = words.index('--db-type')
                    if index+1 < len_words:
                        quant_db_type = words[index+1]
                    else:
                        print "\nERROR: <type> in --db-type <type> is missing in pipeline_args"
                        sys.exit()

    if job_type == 'quantification' and (not found_quant_query or not found_quant_db):
      print "\nERROR: pipeline_args with --db and --query must be specified for quantification jobs"
      sys.exit()

  #############################################################################
  #
  #   verify veracity of args
  #
  #############################################################################

    quant_db_file_ext = ""
    quant_query_file_ext = ""
    if job_type not in all_job_types:
      print "\nERROR: invalid job_type: ", job_type
      print "Allowable job types:"
      for jt in all_job_types:
         print jt
      sys.exit()

    if (input_file == None) and job_type != 'quantification':
      print "\nERROR: input_file must be given for job_type: ", job_type
      sys.exit()

    if (not isfile(input_file)) and (job_type != 'quantification'):
      print "\nERROR: input_file does not exist"
      sys.exit()

    go_qcode = ""
    if job_type=='go':
      if go_quant_file!=None:
        go_qcode = go_quant_file.upper()
      else:
        print "for go jobs, --go_quant_file must be specified and be a valid \
                        quantification file with path, or --go_quant_file OMIT"
        sys.exit()
    if ((job_type=='go') and (go_qcode != 'OMIT')):
      if (not isfile(go_quant_file)):
        print "\nERROR: go_quant_file does not exist"
        sys.exit()

    if(go_qcode == 'OMIT'): 
        print "Running go job without a go_quantification file"

    s = jn_mapping_table.select(jn_mapping_table.c.job_name==job_name)
    results = s.execute()
    row = results.fetchone()
    if row:
      print "\nError job_name: ", job_name, " already used & is not unique"
      sys.exit()

    input_file_ext = ""
    output_file_ext = outputExtensions[job_type]

    if job_type != 'assemble':
      if job_type != 'quantification':
        input_file_ext = inputExtensions[job_type]
        if (input_file.endswith(input_file_ext)):
          pass  
        else:
          print "This job type needs input file ending with ", input_file_ext
          sys.exit()

      else: # running quantification
        if bowtie:   #  db = _project.fasta  query = raw reads (fastq)
          quant_db_file_ext = "_project.fasta"
          quant_query_file_ext = ".fastq"
          if quant_db_file.endswith(quant_db_file_ext):
            pass;
          else:
            print "quant_db_file needs to end with: ", quant_db_file_ext
            sys.exit()
          if quant_query_file.endswith(quant_query_file_ext):
            pass;
          else:
            print "quant_query_file needs to end with: ", quant_query_file_ext
            sys.exit()
        else: # blast : both fasta (user_job.py must require this)
          quant_db_file_ext = ".fasta"
          quant_query_file_ext = ".fasta"
          if quant_db_file.endswith(quant_db_file_ext):
            pass;
          else:
            print "quant_db_file needs to end with ", quant_db_file_ext
            sys.exit()
          if quant_query_file.endswith(quant_query_file_ext):
            pass;
          else:
            print "quant_query_file needs to end with: ", quant_query_file_ext
            sys.exit()
    else:  # running assembler
      if running_mira:  input_file_ext = '.fastq'
      else: 
        if input_file.endswith('fasta'):
          input_file_ext = '.fasta'
        elif input_file.endswith('fastq'):
          input_file_ext = '.fastq'
        else:
          print "for assembler, input file must end with fasta or fastq"
          sys.exit()

    ###########################################################################
    #
    #   enter "user_job" in pn_mapping table and jn_mapping table.
    #
    ###########################################################################

    pname = 'user_job'
    project_id = 0
    s = pn_mapping_table.select(pn_mapping_table.c.project_name==pname)
    results = s.execute()
    row = results.fetchone()
    if row: 
      project_id = row.project_id
    else:
      i = pn_mapping_table.insert()
      i.execute(project_name='user_job')
      s = pn_mapping_table.select(pn_mapping_table.c.project_name==pname)
      results = s.execute()
      row = results.fetchone()
      project_id = row.project_id

    ###########################################################################
    #
    #   get default args
    #
    ###########################################################################

    if (job_type != 'panther'):
      s = default_args_table.select(default_args_table.c.job_type==job_type)
      results = s.execute()
      row = results.fetchone()

      for arg_name in all_arg_names:
        default_args[arg_name] = row[arg_name]

    ###########################################################################
    #
    #  set args_to_use to cmd_line_arg if given else to default value
    #
    ###########################################################################

    args_to_use['pipeline_args'] = ""

    for arg_name in all_arg_names:
      if arg_name in cmd_line_args: 
        if arg_name == 'pipeline_args' and job_type == 'quantification':
          args_to_use[arg_name] = ' --db-type ' + quant_db_type + \
                                                ' --aligner ' +  quant_aligner
        else:
          args_to_use[arg_name] = cmd_line_args[arg_name]
      else:
        if job_type == 'panther':
          args_to_use['loc'] = 'local'
          args_to_use['resources'] = 'cpu:1'
          args_to_use['priority'] = 1
          args_to_use['pipeline_args'] = ""
          args_to_use['executable'] = ""
          args_to_use['process_args'] = ""
        else:
          args_to_use[arg_name] = default_args[arg_name]

    ###########################################################################
    #
    #  create name for output file and append -output_file <output_file> 
    #  to args_to_use['pipeline_args']
    #
    ###########################################################################

    dest_dir = settings.special_runs_dir + job_name
    output_file = dest_dir + "/" + job_name + output_file_ext
    if os.path.exists(dest_dir):
      print "\nError job_name: ", job_name, " already used and is not unique"
      sys.exit()
    os.mkdir(dest_dir)
    args_to_use['pipeline_args'] = args_to_use['pipeline_args'] + \
                                                 " -output_file " + output_file

    ###########################################################################
    #
    #  append -input <input_file> to args_to_use['pipeline_args'] AND copy 
    #  input file to dest_dir for quantification jobs, copy query and db files
    #  to dest_dir and add them to process_args
    #
    ###########################################################################

    quant_query_file_new = ""
    quant_db_file_new = ""

    if job_type == 'go':
      if go_quant_file != 'OMIT':
        go_quant_file_new = output_file + '_quantification.txt'
        shutil.copyfile(go_quant_file,go_quant_file_new)
    if (job_type == 'quantification'):
      quant_db_file_new = dest_dir + '/' + job_name + quant_db_file_ext
      shutil.copyfile(quant_db_file,quant_db_file_new)
      quant_query_file_new = dest_dir + '/' + job_name + quant_query_file_ext
      shutil.copyfile(quant_query_file,quant_query_file_new)
      args_to_use['pipeline_args'] = args_to_use['pipeline_args'] + \
                                           " --query " + quant_query_file_new
      args_to_use['pipeline_args'] = args_to_use['pipeline_args'] + \
                                                 " --db " + quant_db_file_new
    else:
      path,file=os.path.split(input_file)
      file_new = job_name + input_file_ext
      input_file_new = dest_dir + '/' + file_new
      shutil.copyfile(input_file,input_file_new)
      args_to_use['pipeline_args'] = args_to_use['pipeline_args'] + \
                                                    " -input " + input_file_new

    ###########################################################################
    #
    #  append -job_name  <job_name> to args_to_use['pipeline_args']
    #
    ###########################################################################

    args_to_use['pipeline_args'] = args_to_use['pipeline_args'] + \
                                                       " -job_name " + job_name

    ###########################################################################
    #
    #  enter row for job into settings.special_queue
    #
    ###########################################################################

    i = jn_mapping_table.insert()
    i.execute(project_id=project_id,job_type=job_type,job_name=job_name,
                                                      started='N',finished='N')
    results = session.conn.execute(
        "select job_id from jn_mapping where ((job_name = ?) and \
                                        (job_type = ?))", (job_name,job_type,))
    row = results.fetchone()
    job_id = ""
    if not row:
       print "ERROR ", job_type, job_name, " not found in jn_mapping table"
       sys.exit()
    else:
       job_id = row.job_id
    i = queue_table.insert()
    i.execute(project_id=project_id,job_id=job_id,job_type=job_type,
                                              priority=args_to_use['priority'])
    i = args_table.insert()
    i.execute(job_id=job_id,executable=args_to_use['executable'],
              loc=args_to_use['loc'], process_args=args_to_use['process_args'],
              pipeline_args=args_to_use['pipeline_args'],
          resources=args_to_use['resources'],priority=args_to_use['priority'])
    print job_name, ''' successfully inserted into args_table & special_queue
                      & jn_mapping'''

if __name__ == '__main__':
    main()
