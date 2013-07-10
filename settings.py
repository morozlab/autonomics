# -*- coding: utf-8 -*-
'''
Created on Oct 23, 2012

@author: Mathew Citarella
'''

import platform
import os
import re

hname = os.environ['HOST'] 
print hname

machine = ""
match_obj = re.search('oem',hname)
if match_obj:
  machine = 'oem'

match_obj = re.search('acis',hname)
if match_obj:
  machine = 'acis'

match_obj = re.search('hpc',hname)
if match_obj:
  machine = 'hpc'

if machine == "":
  print "settings.py unable to find machine it is running on"
  sys.exit()

class Credentials:
    '''
        This class implements functionality to store credentials for accessing
        remote resources, whether those resources be mail, ftp, or compute
        servers.
        
        Attributes:
            host (str):
                The string representation of the IP address for the remote
                server
            
            user (str):
                The username to use when connecting to the remote server
            
            password (str):
                The password to use when connecting to the remote server
            
            from_file (str)=None:
                The text file (optional) to use when loading credentials into
                this object
    '''

    def __init__(self, host=None, user=None, passwd=None, from_file=None):
        self.host = host
        self.user = user
        self.passwd = passwd
        if(not from_file is None):
            self.load_from_file(from_file)


    def __getitem__(self, item_name):
        if(hasattr(self, item_name)):
            return getattr(self, item_name)
        else:
            return None

    def __setitem__(self, item_name, item_value):
        setattr(self, item_name, item_value)

    def load_from_file(self, path):
        '''
            Loads credentials from the file given by 'path'. 
            
            Credential files should be in the form:
                attribute: value
            
            With each line containing only a single attribute.
        '''
        try:
            fh = open(path, 'r')
            for l in fh:
                el = l.strip().split(":")
                if(len(el) == 2):
                    self[el[0]] = el[1]
            fh.close()

        except Exception as e:
            print("Error reading credentials from file: " + path)
            print(e.message)        

    def update(self, u, p):
        self.user = u
        self.passwd = p

#detect the operating system the software is currently running on
THIS_OS = platform.system()

#The default data directory
#
#All project folders will be created under this directory
#
#Command-line and web-interface submitted jobs can be found under:
#    home_dir/special_jobs/
home_dir = "/srv/data2/pipeline/"
special_runs_dir = "/srv/data2/pipeline/special_runs/"


#Important system paths for finding various executables and data directories.
#
pfam_exec_path = "/srv/data/pfam/PfamScan/"
pfam_data_path = "/scratch/hpc/mcitar/pfam/"
panther_data_path = "/srv/data2/pwilliams/PANTHER8.0/"
khmer_path = "/srv/data2/software/Khmer/scripts/"
trinity_path = "/srv/data2/software/Trinity/"
cutadapt_path = "/usr/local/bin/"
mira_path = "/usr/bin/"

#Base directory to use for temporary files and folders on the remote HPC cluster
#
remote_dir = "/scratch/hpc/mcitar/"

#Main installation directories. INSTALL_DIR is the full path to where the 
#system was cloned from github. 
if machine == 'oem':
    INSTALL_DIR = "/home/oem/python_software/autonomics/"
    PERLPATH = "/home/oem/PerlScripts/Database/"
if machine == 'acis':
    INSTALL_DIR = "/home/pwilliams/python_software/autonomics/"
if machine == 'hpc':
    INSTALL_DIR = "/home/pwilliams/python_software/autonomics/"
    PERLPATH = "/home/pwilliams/PerlScripts/Database/"

#Directory containing all python scripts used by the system. This is the full 
#path to the 'scripts' subdirectory of the Autonomics repository.

python_scripts_path = INSTALL_DIR + "scripts/"

#SCRIPT_BASE describes where to find scripts, relative to the main Autonomics
#installation.
SCRIPT_BASE = "scripts/"
#SCRIPT_PATH is the full path to the scripts directory.
SCRIPTPATH = INSTALL_DIR + SCRIPT_BASE

#CYGWIN_INSTALL_DIR refers to the installation path for an instance of 
#Autonomics running on a cygwin machine.
CYGWIN_INSTALL_DIR = "/cygdrive/c/Documents and Settings/Administrator/My \
Documents/Dropbox/Work/moroz-python/zeroclick/"
#WINDOWS_INSTALL_DIR is the full path for the installation of Autonomics on
#a machine running a flavor of Windows.
WINDOWS_INSTALL_DIR = "C:/Documents and Settings/Administrator/My Documents\
/Dropbox/Work/moroz-python/zeroclick/"
WINDOWS_SCRIPT_PATH = WINDOWS_INSTALL_DIR + SCRIPT_BASE

#Change the installation paths depending on the detected platform.
if("Windows" in THIS_OS):
    INSTALL_DIR = WINDOWS_INSTALL_DIR
    SCRIPTPATH = WINDOWS_SCRIPT_PATH
if("CYGWIN" in THIS_OS):
    INSTALL_DIR = CYGWIN_INSTALL_DIR

#Path for credentials files
CRED_PATH = INSTALL_DIR + "credentials/"
#Path for configuration files
CONFIG_BASE = "proj_config/"
CONFIG_PATH = INSTALL_DIR + CONFIG_BASE

#Names for system queues
#These names correspond to the table names in the Autonomics database storing 
#the respective queues. 
normal_queue = 'quenew'
special_queue = 'quenew_special'

#System storage paths
CYGWIN_NEUROBASE_STORAGE = "//Moroz-500g/disk 2/zeroclick_load_data/"
# NEUROBASE_STORAGE_PATH = "G:/zeroclick_storage/"

if machine == 'oem':
   NEUROBASE_FASTADB_PATH = "/var/www/seq_view/database/"
   NEUROBASE_DATA_PATH = "/data/neurobase_load_data/"
if machine == 'hpc':
   NEUROBASE_FASTADB_PATH = "/var/www/html/neurobase/seq_view/database/"
   NEUROBASE_DATA_PATH = "/data/autonomics_data/"

#System credentials 
mail_cred = Credentials(from_file=CRED_PATH + "pipeline_mail_account")
db_cred = Credentials(from_file=CRED_PATH + "web_config_db")
hpc_cred = Credentials(from_file=CRED_PATH + "submit.hpc.ufl.edu")
zc_cred = Credentials(from_file=CRED_PATH + "128.227.70.246")
creds = {}
creds['128.227.123.35'] = Credentials(from_file=CRED_PATH + "128.227.123.35")
creds['128.227.70.246'] = zc_cred


#time the manager waits between checking job status
mainLoopSleepInterval = 30

#total seconds job must not be seen in HPC queue to be marked as complete
waitBeforeMarkComplete = 120
#do you want to restart jobs that run out of memory?!
restart_jobs = True

#which sequence file extenions does the pipeline support
SUPPORTED_FILETYPES = set([".fasta", ".fa", ".fastq", ".bam"])
PAIRED_EXTENSION = ".end2"

#Default resource allocations for the pipeline
DEFAULT_HPC_CPU = 100
DEFAULT_LOCAL_CPU = 1

#redis db connection details
REDIS_HOST = 'localhost'
REDIS_PORT = 6379

#key for unlocking credentials
CRED_KEY = "tg3wutpP6j7j"

#connection details for the ZC server
ZC_HOST = "128.227.70.246"
ZC_USER = "morozgroup"
ZC_PASSWD = "Whitney2011"
ZC_DB_NAME = "zero_click"

DISPATCHER_SLEEP_INTERVAL = 60

#database stuff
DEFAULT_MYSQL_DRIVER = "mysql+oursql"
MYSQLDB_DRIVER = "mysql+mysqldb"

#email configuration
MAIL_PROVIDER = "gmail.com"
MAIL_ACCOUNT = "morozhpc"
MAIL_ACCOUNT2 = "plw1080"

#should we push configuration data from the web server to zc server
PUSH_CONFIGURATION_DATA = True


####################### process-specific settings #########################
PAIRED_END_FLAG = "--paired-end"
PROCESS_MONITOR_INTERVAL = 60
SSH_RETRIES = 10
KHMER_MIN_SEQS = 10000000

#** resource_server.py
RESOURCE_SERVER_PORT = 8000

#** jobs.py
#mailbox in which new job submissions should be located
checkedMailbox = "INBOX"
#retry interval for reconnecting via SSH or pop3
retryInterval = 10

#which assemblers produce their own quantification data
QUANTIFICATION_ASSEMBLERS = ['mira', 'miranewbler']

MAX_HPC_CORES = 16
#** data_gremlin.py
#Ion Torrent chip-specific assemblers
ION_ASSEMBLERS = {"314": "mira", "316": "mira", "318": "mira", "318C": "mira", 
                  "316D": "mira", "900": "trinity"}
IT_BASE_DIR = "/var/www/"

#MiSeq-specific settings
MISEQ_SAMPLE_CSV = "SampleSheet.csv"

#How do you want to transfer files? Acceptable values are 'rsync' 
#(fast, requires public-key infrastructure) and 'sftp' (slow, uses password, 
#no keys required)
MISEQ_XFER_METHOD = 'sftp'
FTP_XFER_METHOD = 'sftp'
IT_XFER_METHOD = 'sftp'
ZC_XFER_METHOD = 'sftp'
