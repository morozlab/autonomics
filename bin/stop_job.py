#!/usr/bin/env python
'''
Peter L. Williams
May 16, 2013
'''

from autonomics import settings, netutils
from autonomics.file_io import FileExtensions
# from autonomics.queue import queue_all
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

    num_args = len(sys.argv)

    arg_list = sys.argv
    if num_args != 3:
       print "call is: ", arg_list[0], " project_name", "job_type"
       sys.exit()

    project_name = arg_list[1]
    job_type = arg_list[2]
    print "project_name: ", project_name
    print "job_type: ", job_type

    session = netutils.DBSession("localhost", settings.ZC_DB_NAME,
                                settings.db_cred.user, settings.db_cred.passwd)
    pn_mapping_table = netutils.get_table_object("pn_mapping", session)
    jn_mapping_table = netutils.get_table_object("jn_mapping", session)
    quenew_table = netutils.get_table_object("quenew", session)
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
       sys.exit()
    
    res = jn_mapping_table.select(jn_mapping_table.c.project_id).where(and_(jn_mapping_table.c.project_id==project_id,
           jn_mapping_table.c.job_type==job_type)).execute()
    if not res.fetchone():
      job_name = project_name + '_' + job_type
      i = jn_mapping_table.insert()
      i.execute(project_id = project_id,job_type=job_type,job_name=job_name,started='Y',finished='Y',queued='Y')
      print "successfully stopped (by INSERT) job: ", job_type ," for project: ", project_name
    else:
      u = jn_mapping_table.update().where(and_(jn_mapping_table.c.job_type==job_type, jn_mapping_table.c.project_id==project_id)).values(started = 'Y', finished = 'Y', queued = 'Y').execute()
      print "successfully stopped (by UPDATE) job: ", job_type ," for project: ", project_name
      quenew_table.delete().where(and_(quenew_table.c.project_id==project_id, quenew_table.c.job_type==job_type)).execute()


if __name__ == '__main__':
    main()
