#! /usr/bin/env python
'''
Peter L. Williams
May 16, 2013
'''

from autonomics import settings, netutils
from autonomics.file_io import FileExtensions
from autonomics.queue import queue_all
import os
from sqlalchemy import *
from sqlalchemy.sql import functions, select, and_
import argparse
import datetime
import time
import subprocess
import sys
from os.path import isfile, join, isdir, exists
from ftplib import FTP

def main():
    all_job_types = ['adapter_trim', 'quality_trim', 'read_normalization','assemble', 'blast_swissprot', 'go', 'kegg', 'pfam', 'quantification', 'blast_nr','all']

    num_args = len(sys.argv)

    arg_list = sys.argv
    if ((num_args != 3) and (num_args !=4)):
       print "call is: ", arg_list[0], " <project_name> <job_type> <prioity[OPTIONAL]>"
       for job_type in all_job_types: print job_type
       sys.exit()
    project_name = arg_list[1]
    job_type = arg_list[2]
    if (num_args != 4): priority = 1
    else: priority = arg_list[3]
    print "project_name: ", project_name, " job_type: ", job_type, " priority: ",priority

    if job_type not in all_job_types:
       print "invalid job_type: ", job_type
       print "FAILED\nFAILED\nFAILED\nFAILED\nFAILED\nFAILED\nFAILED\n"
       sys.exit()


    has_fastq = 0
    fname = settings.PROJECT_DIR + "/" + project_name + "/" + project_name +'.fastq'
    if os.path.exists(fname): has_fastq = 1

    if ((job_type == 'quantification') and (not has_fastq)):
      print "project has no fastq fastq file and one is needed for quantification"
      sys.exit()

    # setup database connections    
    session = netutils.DBSession("localhost", settings.ZC_DB_NAME,
                                settings.db_cred.user, settings.db_cred.passwd)
    pn_mapping_table = netutils.get_table_object("pn_mapping", session)
    jn_mapping_table = netutils.get_table_object("jn_mapping", session)
    quenew_table = netutils.get_table_object("quenew", session)

    # find project_id
    s = pn_mapping_table.select(
         pn_mapping_table.c.project_id).where(
              pn_mapping_table.c.project_name==project_name)
    results = s.execute()
    row = results.fetchone()
    project_id = 0
    if row:
       project_id = row.project_id            
       print "project_id: ", project_id
    else:
       print "project_name not in pn_mapping_table"
       print "FAILED\nFAILED\nFAILED\nFAILED\nFAILED\nFAILED\nFAILED\n"
       sys.exit()

    # for individual job, set jn_mapping: started='N',finished='N',queued='N'
    if job_type != 'all':
       res = jn_mapping_table.select(jn_mapping_table.c.project_id).where(and_(jn_mapping_table.c.project_id==project_id,
              jn_mapping_table.c.job_type==job_type)).execute()
       job_name = project_name + '_' + job_type
       if not res.fetchone():
         i = jn_mapping_table.insert()
         i.execute(project_id = project_id,job_type=job_type,job_name=job_name,started='N',finished='N',queued='Y',priority=priority)
         print "successfully restarted (by INSERT) ", job_type, " for ", project_name
       else:
         u = jn_mapping_table.update().where(and_(jn_mapping_table.c.job_type==job_type, jn_mapping_table.c.project_id==project_id)).values(priority=priority,started = 'N', finished = 'N', queued = 'Y').execute()
         print "successfully restarted (by UPDATE) ", job_type, " for ", project_name

       # put job in quenew
       res = jn_mapping_table.select().where(and_(jn_mapping_table.c.project_id==project_id,
              jn_mapping_table.c.job_type==job_type)).execute()
       row = res.fetchone()
       jid = row.job_id
       i = quenew_table.insert()
       i.execute(project_id = project_id,job_id=jid,job_type=job_type,priority=priority)
       print "added job to quenew"

       print "project_id = ",project_id," job_id= ", jid, " job_type= ", job_type, " priority= ",priority

    # for all job_types, set jn_mapping: started='N',finished='N',queued='N'
    if job_type == 'all':
      for job_type in all_job_types:
          u = jn_mapping_table.update().where(and_(jn_mapping_table.c.job_type==job_type, jn_mapping_table.c.project_id==project_id)).values(priority=priority,started = 'N', finished = 'N', queued = 'N').execute()
          print "successfully restarted ", job_type, " for ", project_name

if __name__ == '__main__':
    main()
