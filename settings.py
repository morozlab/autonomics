'''

@authors: Mathew Citarella & Peter Williams


'''

#-----------------------------------------------------------------------------------
# SETTINGS YOU NEED TO SET ARE BELOW AT  ======= START OF SETTINGS TO MODIFY =====
#-----------------------------------------------------------------------------------

import platform
import os
import re

debug_pipe = 0

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

# INSTALL_DIR is the full path to where the autonomics dir exists (the one cloned
# from github)  git clone https://github.com/morozlab/autonomics

INSTALL_DIR = os.environ['AUTONOMICS_PATH'] + "/"

# PROJECT_DIR is where your project folders will be created on the local machine
# that runs the autonomics pipeline (see README)

PROJECT_DIR = os.environ['PIPELINE_PATH'] + "/"


# ======= START OF SETTINGS TO MODIFY =====================


# see QUICK START section of README on setting environment vars.

# PATHS FOR FINDING VARIOUS EXECUTABLES AND DATA DIRECTORIES.
# this group of paths are on the local machine running the autonomics pipeline

pfam_exec_path = "/srv/data/pfam/PfamScan/"

khmer_path = "/srv/data2/software/Khmer/scripts/"
# https://github.com/ctb/khmer/blob/master/scripts/normalize-by-median.py
# we use normalize-by-median.py

trinity_path = "/srv/data2/pwilliams/trinityrnaseq_r20131110/"
cutadapt_path = "/usr/local/bin/"
mira_path = "/usr/bin/"
panther_data_path = "/srv/data2/pwilliams/PANTHER8.0/"  # panther is optional

# these paths are on the remote HPC cluster
pfam_data_path = "/scratch/lfs/moroz/pfam/"  # dir where the pfam_scan data files located
remote_dir = "/scratch/lfs/moroz/" #directory where a remote job can create a sub dir

# user name on remote HPC cluster
hpc_user = "plw1080"

# QSUB_OK is used by jobs.py to check for the success of the HPC cluster sub job submission command
# i.e. what "qsub xxx.qsub" returns when it successfully submits the qsub script on the cluster.
# on our HPC cluster qsub returns: <process_id>.moab.ufhpc, e.g. 1303263.moab.ufhpc
# our code checks for success by checking for the presence of 'moab' in what is returned.
# we check for success because there are many different reasons for failure making
# it hard to check for failure.
QSUB_OK = 'moab'

# MYSQL SERVER INFORMATION (for autonomics pipeline database)
ZC_DB_NAME = "zero_click"
ZC_DB_HOST = "localhost"

# next 3 only needed if using data_gremlin
ZC_HOST = "128.227.70.246"
ZC_USER = "zeroclick"
ZC_PASSWD = "XXXXXXXX"
IT_BASE_DIR = "/var/www/" # where data lives on sequencer server

# EMAIL CONFIGURATION (where HPC cluster will send error alerts)
MAIL_PROVIDER = "gmail.com"
MAIL_ACCOUNT = "morozhpc"
MAIL_ACCOUNT2 = "plw1080"

# PATH FOR CREDENTIALS FILES (see .../autonomics/credentials/ for examples)
# do not modify this
CRED_PATH = INSTALL_DIR + "/credentials/"

# SYSTEM CREDENTIALS: 
# set next to name of file in credentials dir with login info for pipeline mysql database

db_cred = Credentials(from_file=CRED_PATH + "autonomics_pipe_db");

# set to name of file in credentials dir with your email info (for error reports from the HPC)
mail_cred = Credentials(from_file=CRED_PATH + "pipeline_mail_account")

# set this to file with login info for your HPC cluster
hpc_cred = Credentials(from_file=CRED_PATH + "hpc")

# only needed if you use the data gremlin
zc_cred = Credentials(from_file=CRED_PATH + "data_gremlin")

# SOME EXAMPLES OF CREDS FILES:
# cd CRED_PATH (see above) (note: files end with no newline)
# cat pipeline_mail_account 
# host:pop.gmail.com
# user:morozhpc
# passwd:XXXXXXXXacis.ufl.edu{pwilliams}1222: cat autonomics_pipe_db
# host:localhost
# user:zeroclick
# passwd:XXXXXXXX%: cat hpc
# host:hipergator.hpc.ufl.edu
# user:plw1080
# passwd:XXXXXXXacis.ufl.edu{pwilliams}1226:

# MAX NUM LOCAL/REMOTE CPUS WE CAN USE AT ONE TIME (depends on your local/hpc_cluster policy/allocation)
MAX_NUM_LOCAL_CPUS = 61
MAX_NUM_HPC_CPUS = 1024

# MAX NUM LOCAL JOBS / REMOTE BLAST_NR JOBS WE CAN RUN AT ONE TIME
# (each job may use multiple cpus)  We limit number blast_nr jobs to 2 as they can take
# a long time and prevent other shorter jobs such as pfam or swissprot from running.
MAX_NUM_LOCAL_JOBS = 5
MAX_NUM_BLAST_NR_JOBS = 4

# you may need to adjust this depending on the efficiency of your hpc cluster
BLAST_NR_MAX_WALL_TIME = "99:00:00"
BLAST_NR_MAX_MEM = "12000mb"

BLAST_SWISSPROT_MAX_WALL_TIME = "48:00:00"
BLAST_SWISSPROT_MAX_MEM = "2000mb"

QUANTIFICATION_MAX_WALL_TIME = "24:00:00"
PFAM_MAX_WALL_TIME = "48:00:00"

# REDIS DB CONNECTION DETAILS (http://redis.io/topics/quickstart)
# these are probably correct for you as they are the defaults
REDIS_HOST = 'localhost'
REDIS_PORT = 6379

#  set to 0 if not using neurobase browser
USING_NB = 1

if USING_NB:
  NEUROBASE_FASTADB_PATH = os.environ['NEUROBASE_FASTADB_PATH'] + "/"
  NEUROBASE_DATA_PATH = os.environ['NEUROBASE_DATA_PATH'] + "/"

# ======= END OF SETTINGS TO MODIFY =====================

# special_runs_dir = "/srv/data2/pipeline/special_runs/"

SCRIPTPATH = INSTALL_DIR + "/scripts/"

#Path for configuration files
CONFIG_BASE = "/proj_config/"
CONFIG_PATH = INSTALL_DIR + CONFIG_BASE

#These names correspond to the table names in the Autonomics database storing the respective queues. 
normal_queue = 'quenew'
special_queue = 'quenew_special'

#time the manager waits between checking job status
mainLoopSleepInterval = 120

#total seconds job must not be seen in HPC queue to be marked as complete
waitBeforeMarkComplete = 120
#do you want to restart jobs that run out of memory?!
restart_jobs = True

#which sequence file extenions does the pipeline support
SUPPORTED_FILETYPES = set([".fasta", ".fa", ".fastq", ".bam"])
PAIRED_EXTENSION = ".end2"

DISPATCHER_SLEEP_INTERVAL = 30

#database stuff
DEFAULT_MYSQL_DRIVER = "mysql+oursql"
MYSQLDB_DRIVER = "mysql+mysqldb"
MYSQLDBNAME = 'moroz_lab'

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

#How do you want to transfer files? Acceptable values are 'rsync' 
#(fast, requires public-key infrastructure) and 'sftp' (slow, uses password, 
#no keys required)
MISEQ_XFER_METHOD = 'sftp'
FTP_XFER_METHOD = 'sftp'
IT_XFER_METHOD = 'sftp'
ZC_XFER_METHOD = 'sftp'
