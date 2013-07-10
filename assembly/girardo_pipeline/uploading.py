import argparse
import os
from subprocess import call

parser = argparse.ArgumentParser(description='Uploading script for the zero-click pipeline')
parser.add_argument('-c', '--cleanup-project', dest='cleanUpFlag', action='store_true', help='Flag to tell this script to cleanup the project')
parser.add_argument('-s', '--send-project', dest='sendFlag', action='store_true', help='Flag to tell this script to send the project for loading to the database')
parser.add_argument('projName',  help='Unique identifier (name) for this project. Should be consistently used across the pipeline.')

args = parser.parse_args()

os.chdir("/srv/data2/pipeline/" + args.projName + "/")

fh = open("listFile.txt", 'w')

if(args.sendFlag):
    files=[filename for filename in os.listdir(".") if filename[0] != '.']
    fh = open("listFile.txt", 'w')
    for filename in files:
        fh.write(filename + "\n")

    
fh.close()
        
call("curl -T listFile.txt ftp://ftpguest:7002ram@74.252.103.104", shell=True)

    
    
    


