'''

@author: Mathew Citarella

This module loads data into the comparative Neurogenomics database, as such data arrives, if all dependencies are satisfied.

Neurobase requires that the assembly sequences are loaded first, followed by quantification for the assembled transcripts. Any data file can then be loaded.

usage: neurobase_uploader.py [-h] [--sleep-interval SLEEP_INTERVAL]
                             [--zc-dbhost ZC_DB_HOST] [--zc-dbuser ZC_DB_USER]
                             [--zc-dbpass ZC_DB_PASSWD]

optional arguments:
  -h, --help            show this help message and exit
  --sleep-interval SLEEP_INTERVAL
  --zc-dbhost ZC_DB_HOST    host address of the MySQL database hosting neurobase, defaults in system settings.
  --zc-dbuser ZC_DB_USER    user to use while connecting to MySQL database, defaults in system settings.
  --zc-dbpass ZC_DB_PASSWD    password to use while connecting to MySQL database, defaults in system settings.
  
    This module loops forever, sleeping for SLEEP_INTERVAL seconds before each loop.

'''

import argparse
import hashlib
import os
import pickle
import shutil
import subprocess
import sys
import time
from sqlalchemy.sql import select, and_
from autonomics import settings, netutils, utility, file_io
from autonomics.file_io import FileExtensions


current_project_name = None


def file_loaded(pid, job_type, chksum, session):
    loaded = netutils.get_table_object("loaded_files", session)
    res = loaded.select(and_(loaded.c.project_id==pid,
                        loaded.c.job_type==job_type,
                        loaded.c.checksum==chksum)).execute()
    
    row = res.fetchone()
    if(row is None):
        return False
    return True


def init_project(base_cmd):
    proc = subprocess.Popen(base_cmd + " --init-project -v", shell=True)
    proc.wait()
    return proc


def load_go(base_cmd): 
    proc = subprocess.Popen(base_cmd + " --GO -v ", shell=True)
    proc.wait()
    utility.die_on_error(proc.returncode)
    proc = subprocess.Popen(base_cmd + " --GOCategories -v", shell=True)
    proc.wait()
    return proc


def load_kegg(base_cmd):
    proc = subprocess.Popen(base_cmd + " --KEGG -v", shell=True)
    proc.wait()
    return proc


def load_nr(base_cmd):
    proc = subprocess.Popen(base_cmd + " --annotation --nr -v ", shell=True)
    proc.wait()
    utility.die_on_error(proc.returncode)
    proc = subprocess.Popen(base_cmd + " --alignments --nr -v ", shell=True)
    proc.wait()
    return proc
    

def load_pfam(base_cmd):
    proc = subprocess.Popen(base_cmd + " --pfam -v", shell=True)
    proc.wait()
    return proc


def load_quantification(base_cmd):
    proc = subprocess.Popen(base_cmd + " --quantification -v ", shell=True)
    proc.wait()
    return proc


def load_sequences(base_cmd):
    #does raw data exist and need to be moved?
    raw_base = current_project_name + "/" + current_project_name
    raw_exts = ['.fasta', '.fastq']
    raws = [raw_base + ext for ext in raw_exts]
 
    for raw_file in raws:
        if(os.path.isfile(settings.NEUROBASE_DATA_PATH + raw_file)):
            shutil.move(settings.NEUROBASE_DATA_PATH + raw_file, settings.NEUROBASE_STORAGE_PATH + raw_file)
 
    proc = subprocess.Popen(base_cmd + " --seqs -v ", shell=True)
    proc.wait()
    return proc


def load_swissprot(base_cmd):
    proc = subprocess.Popen(base_cmd + " --annotation --swissprot -v ", shell=True)
    proc.wait()
    utility.die_on_error(proc.returncode)
    proc = subprocess.Popen(base_cmd + " --alignments --swissprot -v ", shell=True)
    proc.wait()
    return proc


def job_loaded_bytype(pid, job_type, session):
    loaded = netutils.get_table_object("loaded_files", session)
    res = select([loaded], and_(loaded.c.job_type==job_type,
                      loaded.c.project_id==pid)).execute()
    row = res.fetchone()
    
    if(row is None):
        return False
    else:
        return True


def mark_file_loaded(pid, job_type, file_path, session, chksum=None):
    if(chksum is None):
        chksum = file_io.sha1_file(file_path)

    loaded_files = netutils.get_table_object('loaded_files', session)
    i = loaded_files.insert().values(project_id=pid,
                                     job_type=job_type,
                                     checksum=chksum,
                                     file_path=file_path
        )
    session.conn.execute(i)
    
 
def mark_loaded(pid, job_type, session):
    jid = netutils.get_jid(pid, job_type, session)
    loaded = netutils.get_table_object("loaded_jobs", session)
    i = loaded.insert().values(job_id=jid)
    i.execute()
      

def ok_to_load(pn, jtype, project_path, file_hash, session):
    pdir = netutils.get_table_object('project_directory', session)
    #make sure there is an initiated project
    s = pdir.select(pdir.c.project_name==pn)
    res = session.conn.execute(s)
    row = res.fetchone()
    if(row is None):
        print("Warning - no pid found for: " + pn)
        return False

    pid = row.projectID
    
    if(jtype != 'assemble'):
        #check if assemble is loaded
        assem_loaded = job_loaded_bytype(pid, 'assemble', session)
        #check if quantification is loaded
        quant_loaded = job_loaded_bytype(pid, 'quantification',
                                         session)
            
        if(jtype != 'quantification'):
            if(not assem_loaded or not quant_loaded):
                print("Quantification not loaded")
                return False
        else:
            print("Quantification job.")
            if(not assem_loaded):
                print("Assembly not loaded.")
                return False

    #file operations
    for data_file, chksum in file_hash.items():
        #check if this file is already loaded
        if(file_loaded(pid, jtype, chksum, session)):
            print("File: " + data_file + " is already loaded.")
            return False

        d, f = os.path.split(data_file)
        path = project_path + "/" + f
        #check if the appropriate file exists
        if(not os.path.exists(path)):
            print("Expected input file does not exist.")
            return False
    
        #do the checksums match?
        if(not file_io.sha1_file(path) == chksum):
            print("Checksums don't work out.")
            return False
    
    return True


def main():
    
    global current_project_name

    parser = argparse.ArgumentParser()
    parser.add_argument('--sleep-interval', dest='SLEEP_INTERVAL',
                        default=settings.mainLoopSleepInterval)
    parser.add_argument('--dbhost', dest='DB_HOST',
                        default='localhost')
    parser.add_argument('--dbuser', dest='DB_USER',
                        default='zeroclick')
    parser.add_argument('--dbpass', dest='DB_PASSWD',
                        default='Whitney2011')
    parser.add_argument('--monitor-directory', dest='monitored_dir',
                        default=settings.NEUROBASE_DATA_PATH)
    args = parser.parse_args()

    out_exts = FileExtensions.output_exts
    load_order = ['assemble', 'quantification', 'blast_swissprot', 'blast_nr', 'go', 'kegg', 'pfam']
    
    load_job = {"quantification": load_quantification, "blast_swissprot": load_swissprot,
                 "blast_nr": load_nr, "go": load_go, "kegg": load_kegg, "pfam": load_pfam, 
                 "assemble": load_sequences}
    
    UPLOAD_SCRIPT = "\"" + settings.SCRIPTPATH + "upload.py\""
    UPLOAD_CMD = "".join(["python ", UPLOAD_SCRIPT, " --user ",
                          args.DB_USER, " --password ",
                          args.DB_PASSWD, " --data-dir ",
                          args.monitored_dir])


    nb_session = netutils.DBSession(args.DB_HOST, "moroz_lab",
                                    args.DB_USER, args.DB_PASSWD,
                                    driver=settings.MYSQLDB_DRIVER)

    while(True):        
        
        #set up some table objects
        proj_dir = netutils.get_table_object("project_directory",
                                           nb_session)
        
        #create a list of job types
        #results = select([default_config.c.job_type]).execute()
        #job_types = [row.job_type for row in results.fetchall()]

        '''neurobase_uploader.py
        ---------------------
            - wake up
            - iterate over directories in the landing directory for data to be uploaded
            - open the checksum.txt file in each directory
            - read the list of filenames and their checksums.
            - if checksum of corresponding file matches that in list
              - for each dependency, if job_id loaded or job_id not configured
              - load data file
              - mark job_id loaded'''
          
        #iterate over uploaded files
        projects = os.listdir(args.monitored_dir)
        for project in projects:
            
            current_project_name = project
            base_cmd = UPLOAD_CMD + " -p " + project + " --public-name " + project
            storage_path = settings.NEUROBASE_STORAGE_PATH + project + "/"
            project_path = args.monitored_dir + project + "/"
            pickle_name = project + "_checksum.pickle"
            pickle_path = project_path + pickle_name
            
            #check if this directory only has a checksum pickles
            files = os.listdir(project_path)
            if(len(files) == 1):
                if(files[0] == pickle_name and (project_path != "/" or project_path != "C:\\")):
                    shutil.rmtree(project_path)
                    continue
            elif(len(files) == 0):
                shutil.rmtree(project_path)
                continue
            
            #check if a storage directory exists for this project
            #if(not os.path.exists(storage_path)):
                #os.mkdir(storage_path)
                           
            #do we need to initialize this project?
            nb_pid_query = proj_dir.select(proj_dir.c.project_name==project).execute()
            nb_row = nb_pid_query.fetchone()
            if(nb_row is None):
                init_project(base_cmd)
                s = proj_dir.select(
                        proj_dir.c.project_name==project
                    )
                res = nb_session.conn.execute(s)
                nb_row = res.fetchone()

            pid = nb_row.projectID

            #get the project id for this project
            #results = pn_mapping.select(pn_mapping.c.project_name==project).execute()
            #pid = results.fetchone()
            #if(pid is None):
            #    print("No project entry in pn_mapping for: " + project + ", skipping.")
            #    continue

            #pid = pid.project_id
            
            #pickle-in the dictionary of filename => checksums
            if(not os.path.exists(pickle_path)):
                print("No checksum pickle found for: " + project + ", skipping.")
                continue
            
            fh = open(pickle_path, 'r')
            chk_hash = pickle.load(fh)
            fh.close()
            
            #iterate over job_types in load_order, checking the checksums of each file
            for job_type in load_order:
                print(project + " " + job_type)
                #first, check if job_type is in the checksum hash
                if(not job_type in chk_hash):
                    continue
                
                print("Checking load preconditions for " + project + " " + job_type)
                file_hash = chk_hash[job_type]
                #check if this job_is ok to load
                if(ok_to_load(project, job_type, project_path,
                              file_hash, nb_session)):
                    print("Loading data for " + project + ", job type: " + job_type)
                    #load based on the job_type
                    load_process = load_job[job_type](base_cmd)
                    utility.die_on_error(load_process.returncode) 
                    #remove the file
                    for (data_path, chksum) in file_hash.items():
                        d, f = os.path.split(data_path)
                        os.remove(project_path + f)
                        mark_file_loaded(pid, job_type, data_path,
                                         nb_session, chksum)
                            

        sys.stdout.write("Sleeping for: " + str(args.SLEEP_INTERVAL) + "\n")
        time.sleep(args.SLEEP_INTERVAL)

if __name__ == "__main__":
    main()
