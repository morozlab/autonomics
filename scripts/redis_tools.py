'''
Created: 8/2012

Author: Mathew Citarella

This module is designed to be called as a script from the command line, and includes functionality for loading kegg and GO associations into a redis keystore.

usage: redis_tools.py [-h] [--koa KOA_FILE] [--goa GOA_FILE] [-b BATCH_SIZE]

This script will load the koa and goa files into a redis keystore. Assumes a
non-password protected redis instance running on 127.0.0.1 port 6379

optional arguments:
  -h, --help            show this help message and exit
  --koa KOA_FILE        File of swissprot acc -> ko associations.
  --goa GOA_FILE        File of swissprot acc -> go term associations.
  -b BATCH_SIZE, --batch-size BATCH_SIZE
                        Number of set commands in each transaction with the
                        db. Default: 10000
'''

import argparse
import redis

def batch_submit_tab(file_path, pipe, pipe_method, batch_size, key_prefix, key_columns):
    num_in_batch = 0
    batch_no = 1
    fh = open(file_path, 'r')
    for line in fh:
        line = line.strip()
        el = line.split("\t")
        key = ":".join([el[c] for c in key_columns])
        pipe_method(key_prefix + key, line)   
        num_in_batch += 1 
        if(num_in_batch == batch_size):
            pipe.execute()
            num_in_batch = 0
            batch_no += 1
            
    fh.close()

def main():
    parser = argparse.ArgumentParser(description="This script will load the koa and goa files into a redis keystore. Assumes a non-password protected redis instance running on 127.0.0.1 port 6379")

    parser.add_argument('--koa', dest='koa_file', default=None, help='File of swissprot acc -> ko associations.')
    parser.add_argument('--goa', dest='goa_file', default=None, help='File of swissprot acc -> go term associations.')
    parser.add_argument('-b', '--batch-size', dest='batch_size', default=10000, help='Number of set commands in each transaction with the db. Default: 10000')
    
    args = parser.parse_args()
    
    r = redis.Redis(host='localhost', port=6379, db=0)
    pipe = r.pipeline()
    batch_size = args.batch_size
    
    if(args.koa_file):
        batch_submit_tab(args.koa_file, pipe, pipe.lpush, batch_size, "zero_click:koa:", [0])    
    if(args.goa_file):
        batch_submit_tab(args.goa_file, pipe, pipe.lpush, batch_size, "zero_click:goa:", [0])
    
        

if __name__ == '__main__':
    main()