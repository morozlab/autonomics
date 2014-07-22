#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      mat
#
# Created:     08/08/2012
# Copyright:   (c) mat 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from alignment.io import Reader
# from Bio.Blast import NCBIStandalone
from Bio import SeqIO
from sqlalchemy import *
from sqlalchemy.sql import functions
from sqlalchemy.sql import *
from autonomics.file_io import Record
from autonomics import settings, netutils, utility
import argparse
import datetime
import os
import MySQLdb
import re
import sqlalchemy
import subprocess
import sys
import time
import utils
from datetime import date
import warnings
import datetime

def main():
    parser = argparse.ArgumentParser(description = "This script contains various methods of loading data into NeuroBase.")
    parser.add_argument('--user', dest='user', required=True, help='Username for mysql database containing NeuroBase data')
    parser.add_argument('--password', dest='passwd', required=True, help='Password for mysql database containing NeuroBase data')
    parser.add_argument('--host', dest='host', required=True, help='Password for mysql database containing NeuroBase data')
    parser.add_argument('--db', dest='db', required=True, help='Password for mysql database containing NeuroBase data')
    parser.add_argument('--outfile', dest='outfile', required=True, help='file with data')
    args = parser.parse_args()

    session = utils.DBSession(args.user, args.passwd, args.host, args.db)
    session.activateConnection()
    query = ""
    try: 
        query = "LOAD DATA LOCAL INFILE '" + args.outfile + "' REPLACE INTO TABLE sorted_homology"        
#        query = "LOAD DATA INFILE '" + args.outfile + "' REPLACE INTO TABLE sorted_homology"
#        print "try-ing QUERY: ", query
        with warnings.catch_warnings():
           warnings.simplefilter("ignore")
           session.conn.execute(query)
#        print "exited try with success"
    except:
#        print "FAILED: ", query
        query = "LOAD DATA INFILE '" + args.outfile + "' REPLACE INTO TABLE sorted_homology"
#        query = "LOAD DATA LOCAL INFILE '" + args.outfile + "' REPLACE INTO TABLE sorted_homology"        
#        print "running except: ", query
        with warnings.catch_warnings():
          warnings.simplefilter("ignore")
          session.conn.execute(query)
#        print "exited except with success"
if __name__ == '__main__':
    main()
