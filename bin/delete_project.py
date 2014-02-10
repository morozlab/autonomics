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
    if num_args != 2:
       print "call is: ", arg_list[0], " project_name"
       sys.exit()

    project_name = arg_list[1]
    print "project_name: ", project_name
    
    session = netutils.DBSession("localhost", settings.ZC_DB_NAME,
                                settings.db_cred.user, settings.db_cred.passwd)
    pn_mapping_table = netutils.get_table_object("pn_mapping", session)
    jn_mapping_table = netutils.get_table_object("jn_mapping", session)
    run_stats_table = netutils.get_table_object("run_stats", session)
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
    u = jn_mapping_table.delete().where(jn_mapping_table.c.project_id==project_id).execute()
    u = pn_mapping_table.delete().where(pn_mapping_table.c.project_id==project_id).execute()
    u = run_stats_table.delete().where(run_stats_table.c.project_id==project_id).execute()
    print "successfully deleted: ", project_name

if __name__ == '__main__':
    main()
