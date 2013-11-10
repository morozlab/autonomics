'''
Created on Nov 28, 2012

@author: Mathew Citarella

Data_gremlin is responsible for the following tasks:

1) Setting initial configurations in the MySQL database for newly discovered sequencing runs. 
2) During configuration, determining the sequencing type (if possible) and setting the appropriate pre-assembly and assembly tasks to be run
3) Moving data from the sequence sources to wherever the manager and dispatcher processes are running when the data is ready.
4) Pushing configuration data from the web server to the compute server (if they are not the same machine)
5) Setting new runs as configured, downloaded, etc

This module configures new runs with the following logic:
    - check runname_to_pid to see if the run_name already exists in the table. If so, the run is configured. Skip.
    - check processed_runs to see if a (source_id, run_name) tuple exists in the table. If so, the run was configured at one point, and someone removed it from the system on purpose. Skip.
    
    At this point, we've determined the run is ok to configure.
    - enter a project_name for the run using run_name in pn_mapping, returning a project_id
    - use the project_id to add jobs for each job in the default_configuration table. Get a list of job_ids for those jobs.
    - use the list of job_ids and project_ids to enter a entry for the project in the configuration table for each job in default_configuration.
    - determine if the run is paired end, set the appropriate flags for all jobs in the args table
   - determine which assembly pipeline to use, and configure the pre_assemble jobs accordingly.
    
Data is moved with the move_<source_type>_data family of methods.

Configuration is pushed from the web server to the autonomics compute server with push_configuration. 

Configuration data includes (but is not limited to):
    - project_name
    - job_ids for all configured jobs
    - custom arguments for configured jobs (pipeline_args, process_args, etc)
    - whether the job is paired-end, has been configured, or downloaded
    
Configuration data is pushed to the autonomics compute server so that the MySQL instance there can be maintained as a mirror of the database intsance on the web server in the presence of firewalls. Without the firewall, the autonomics server would access the web server MySQL instance directly.

'''

import argparse
import netutils
import os
import re
import settings
import sys
import time
from autonomics.file_io import Record
from autonomics.settings import Credentials
from sqlalchemy import select, update
from sqlalchemy.sql import and_, or_


def ensure_project_paired(run_name, source_id, session):
    rn_2_pid = netutils.get_table_object('runname_to_pid', session)
    s = rn_2_pid.select(rn_2_pid.c.run_name==run_name,
                        rn_2_pid.c.source_id==source_id)
    res = session.conn.execute(s)
    row = res.fetchone()
    if(not row.paired_end):
        mark_run_paired(run_name, source_id, session)


def fix_miseq_path_4sftp(path):
    ''' path:
            Path that should be fixed for the ssh server installed on the MiSeq.
            
        Returns a path where /cygdrive is spliced out, and /d/ is replaced with /D/. Only necessary for certain cygwin-based SSH servers.
    '''
    return path.replace("/cygdrive", "").replace("/d/", "/D/")


def get_name_mapping(run_2_pid, pn_mapping, run_name):
    '''
    run_2_pid:
        Table object reflected using SQLAlchemy. Used in generative and non-generative query methods.
        
    pn_mapping:
        Same as run_2_pid, but of the pn_mapping table.
        
    run_name:
        Name of the run you wish to get the project_name for, as a string.
    
    Returns the project name for the specified run if one exists in pn_mapping, returns None otherwise.
    '''
    results = select([run_2_pid, pn_mapping], and_(run_2_pid.c.project_id == pn_mapping.c.project_id, run_2_pid.c.run_name == run_name)).execute()
    mapping = results.fetchone()
    if(mapping is None):
        return None
    else:
        return(mapping.project_name)


def get_runname(s, src_type):
    '''
    s:
        String that contains the run_name
    src_type:
        String that identifies the source type associated with the string. Directs the method how to parse the run_name from s.
    
    Splits the run_name out of a string, depending on which src_type is expected.
    
    src_types:
        ion_torrent - expecting a path to a fastq file, returns the basename of the file as run_name
    '''
    if(src_type == "ion_torrent"):
        #split the path to dir, filename; take filename and return the part without the ext as the runname
        return os.path.splitext(os.path.split(s)[1])[0]
    
def get_configuration_tables(zc_session):
    '''
    Returns a tuple of commonly-used configuration tables. 
    
    Currently tuple: (rn_2_pid, pn_mapping, jn_mapping, config, default)
    '''
    default = netutils.get_table_object("default_configuration", zc_session)
    config = netutils.get_table_object("configuration", zc_session)
    pn_mapping = netutils.get_table_object("pn_mapping", zc_session)
    jn_mapping = netutils.get_table_object("jn_mapping", zc_session)
    rn_2_pid = netutils.get_table_object("runname_to_pid", zc_session)
    return (rn_2_pid, pn_mapping, jn_mapping, config, default)


def delete_pending(run_name, source_id, session):
    ''' run_name:
            The name of the run you wish to check.
            
        source_id:
            The id of the source that produced run_name.
            
        session:
            A netutils.DBSession object with an active connection.
            
        Looks into the runname_to_pid and pn_mapping tables to determine if a delete operation is pending for the run in question.
        
        Current conditions for a delete pending: The project name assigned to the run starts with either 'delete' or "Delete'.
        
    '''
    rn_2_pid = netutils.get_table_object("runname_to_pid", session)
    pn_mapping = netutils.get_table_object("pn_mapping", session)
    
    s = select([pn_mapping.c.project_name], and_(rn_2_pid.c.project_id==pn_mapping.c.project_id, rn_2_pid.c.run_name==run_name))
    res = session.conn.execute(s)
    row = res.fetchone()
    if(row.project_name.startswith("delete") or row.project_name.startswith("Delete")):
        return True
    
    return False
    
    
def init_ftp_runs(host, source_id, zc_session, credentials, path):
    '''
    host:
        String representation of the IP address for the FTP server. Will be connected to using a netutils.SSHConnection object.
        
    source_id:
        Integer ID of the FTP server in the system.
        
    credentials:
        Credentials object created by read_credentials.
        
    path:
        Path on the FTP server to the folder containing the input data.
    
    Writes default configuration for detected runs in the ftp server specified by host. Uses credenials to connect to a remote ftp server and detect data.
    
    Only writes configurations for data that has not yet been configured, and does not have an entry in processed_runs table. This is to prevent runs that were configured and later deleted by users from being added again to the system.
    
    Searches for data files with extensions matching settings.SUPPORTED_FILETYPES in the 'path' directory on 'host', and configures them one-by-one. Updates the run as configured in runname_to_pid and adds a (source_id, run_name) entry in processed_runs.   
   
    '''

    if(not "localhost" in host):
        # open an ssh connection to the remote computer
        ssh = netutils.SSHConnection(host, username=credentials['user'], password=credentials['passwd'])
    
    else:
        # check if the path exists
        if(not os.path.exists(path)):
            print("Data source " + str(source_id) + " @ " + host + " misconfigured.")
            return
        
        # if the path exists, read the files there
        files = os.listdir(path)
        for f in files:
            filename = os.path.split(f)[1]
            #full_path = d + "/" + filename
            base, ext = os.path.splitext(filename)
            if(not ext in settings.SUPPORTED_FILETYPES):
                continue
                
            # check if this run is paired-end
            paired = False
            paired_file = path + "/" + f + settings.PAIRED_EXTENSION
            if(os.path.exists(paired_file)):
                paired = True

            # the basename is the new default project name
            if(not run_configured(base, source_id, zc_session) and not run_processed(base, source_id, zc_session)):
                
                if(not write_configuration(base, source_id, zc_session,  paired=paired)):
                    print("Did not write configuration for: " + base)
                #else:
                #    mark_processed(base, source_id, zc_session)
            else:
                #if this job is paired-end, make sure it's set as such
                if(paired):
                    ensure_project_paired(base, source_id, zc_session)
               


def init_it_runs(host, source_id, zc_session, credentials, **kwargs):
    '''
        host (str):
            String representation of the IP address for the FTP server. Will be connected to using a netutils.SSHConnection object.
        
        source_id (int):
            ID of the FTP server in the system.
        
        credentials (settings.Credentials):
            Credentials object created by read_credentials.
    
        Writes default configuration for Ion Torrent and Ion Proton runs detected on the associated Ion Servers.
    
        Default configuration is written as outlined in the write_configuration method documentation. This method performs the additional steps of checiking for the chip type used in sequencing and setting the appropriate assembler. Runs from the Ion Torrent PGM are set to use the     
    '''
    # create DB sessions/table objects
    session = netutils.DBSession(h=host, u=credentials['dbuser'], p=credentials['dbpasswd'], d='iondb', driver=POSTGRES_DRIVER)
    rundb_results = netutils.get_table_object("rundb_results", session)
    rundb_experiment = netutils.get_table_object("rundb_experiment", session)
    # seq_files = netutils.get_table_object("system_input", zc_session)

    args = netutils.get_table_object("args", zc_session)
    default_args = netutils.get_table_object("default_args", zc_session)
    default = netutils.get_table_object("default_configuration", zc_session)
    config = netutils.get_table_object("configuration", zc_session)
    pn_mapping = netutils.get_table_object("pn_mapping", zc_session)
    jn_mapping = netutils.get_table_object("jn_mapping", zc_session)
    rn_2_pid = netutils.get_table_object("runname_to_pid", zc_session)
    
    assemble_replacer = re.compile("--assembler\s\w+")

    # select the completed runs that aren't assembled or configured
    results = select([rundb_results.c.fastqLink, rundb_results.c.experiment_id]).execute()
    for row in results.fetchall():
        # check if this run is already in the system
        base_name = get_runname(row.fastqLink, "ion_torrent")
        exid = row.experiment_id
        res = rn_2_pid.select(rn_2_pid.c.run_name == base_name).execute()
        row = res.fetchone()
        if(not row is None):
            continue
        
        # check if this run has been processed
        if(run_processed(base_name, source_id, zc_session)):
            continue

        # insert an entry into pn_mapping
        s = pn_mapping.select(pn_mapping.c.project_name == base_name)
        i = pn_mapping.insert().values(project_name=base_name)
        pid = netutils.check_and_insert(s, i)
        
        if(pid is None):
            print("Project already configured with this name but not with this run, skipping.")
            continue

        # insert an entry in rn_2_pn
        s = rn_2_pid.select(and_(rn_2_pid.c.run_name == base_name, rn_2_pid.c.source_id == source_id))
        i = rn_2_pid.insert().values(source_id=source_id, run_name=base_name, project_id=pid)
        netutils.check_and_insert(s, i)

        # insert entries into jn_mapping
        res = default.select().execute()
        jids = {}
        for row in res.fetchall():
            s = config.select(and_(config.c.project_id == pid, config.c.job_type == row.job_type))
            i = config.insert().values(project_id=pid, job_type=row.job_type, code=row.code)
            netutils.check_and_insert(s, i)

            s = jn_mapping.select(and_(jn_mapping.c.project_id == pid, jn_mapping.c.job_type == row.job_type))
            i = jn_mapping.insert().values(project_id=pid, job_type=row.job_type)
            jid = netutils.check_and_insert(s, i)
            
            jids[row.job_type] = jid
            
        # update rundb_results to mark configured
        u = rn_2_pid.update().where(and_(rn_2_pid.c.source_id == source_id,
                                         rn_2_pid.c.run_name == base_name)).values(configured='Y')
        u.execute()
        
        # check the chip type and get the assembler
        chip_res = select([rundb_experiment.c.chipType], rundb_experiment.c.id == exid).execute()
        row = chip_res.fetchone()
        chipType = row.chipType.replace("\"", "")
        chipType = chipType.replace("R", "")
        assembler = settings.ION_ASSEMBLERS[chipType]
        
        update_preassemble_config(pid, assembler, zc_session)
        
        # set the appropriate pipeline_args for the assembly job
        res = default_args.select(default_args.c.job_type == 'assemble').execute()
        row = res.fetchone()
        pipe_args = row.pipeline_args
        pipe_args = assemble_replacer.sub("--assembler " + assembler, pipe_args)
        
        i = args.insert().values(job_id=jids['assemble'], loc=row.loc,
                                 executable=row.executable, process_args=row.process_args, pipeline_args=pipe_args,
                                 resources=row.resources, priority=row.priority)
        i.execute()
        
        # mark that this run has been processed
        mark_processed(base_name, source_id, zc_session)
        
        # s = seq_files.select(and_(seq_files.c.source_id==source_id, seq_files.c.base_name==base_name))
        # i = seq_files.insert().values(source_id=source_id, base_name=base_name, configured='N')
        # insert a record in sequence_files if this run hasn't been configured and isn't already in the table
        # netutils.check_and_insert(s, i)
        # fh.write(get_runname(row.fastqLink, "ion_torrent") + "\n")

    session.close()


def init_miseq_runs(host, source_id, zc_session, credentials, path):
    ''' host (str):
            String representation of the IP address for this MiSeq data source
            
        source_id (int):
            ID of the MiSeq data source in the data_sources MySQL table
            
        zc_session (netutils.DBSession):
            Active session with the Autonomics database
            
        credentials (settings.Credentials):
            Object representing the credentials necessary to connect to this MiSeq data source
            
        path (str):
            Path to the root data directory on this data source
            
        This method detects new sequencing runs on a registered MiSeq data source. Newly detected runs are added to the MySQL database as a project configured with default settings.
    ''' 
    #open an ssh connection
    ssh = netutils.SSHConnection(host, username=credentials['user'], password=credentials['passwd']) 
    
    #get a list of runs on the MiSeq to iterate over
    files = ssh.execute("ls " + path)
    run_names = []
    for line in files:
        if(line.startswith(".") or line.startswith("..") or "MiSeq" in line or "Temp" in line):
            continue
        line = line.rstrip("\n")
        run_names.append(line)

    #Cygwin SFTP resolves paths differently than ssh ls, so we need to use a slightly different path format
    #In a future release, need to convert path to the below path programmatically
    path = fix_miseq_path_4sftp(path)
    
    for run_name in run_names:
        remote_path = path + run_name + "/" 
        
        #get the needed project information
        tmp_csv = run_name + settings.MISEQ_SAMPLE_CSV
        ssh.get(remote_path + settings.MISEQ_SAMPLE_CSV, tmp_csv)
        fh = open(tmp_csv, 'r')
        project_name = None
        for l in fh:
            l = l.rstrip("\n")
            if(l.startswith("Experiment Name")):
                els = l.split(",")
                project_name = els[1].replace(" ", "_")
                    
        fh.close()
        os.remove(tmp_csv)
        
        #check if this run has already been taken care of        
        if(project_name is None or run_processed(run_name, source_id, zc_session)):
            continue
        
        if(not write_configuration(run_name, source_id, zc_session, paired=True, project_name=project_name)):
            print("Did not write configuration for: " + project_name)
        else:
            mark_configured(run_name, source_id, zc_session)
            mark_processed(run_name, source_id, zc_session)
            
    ssh.close()


def make_arg_string(row_obj):
    '''
        row_obj:
            A RowProxy object obtained via the fetchall or fetchone method of a SQLAlchemy result proxy.  Should be the result of performing a select on the args or default_args table, and have the attributes listed in row_attrs. 
    
        Converts the SQLAlchemy RowProxy object into a process_arg string or a pipeline_arg string. 
    
        row_attrs = ["executable", "loc", "process_args", "pipeline_args", "resources", "priority"]
    '''
    
    if(row_obj is None):
        return ""
        
    row_attrs = ["executable", "loc", "process_args", "pipeline_args", "resources", "priority"]
    s = []
    for attr in row_attrs:
        s.append(str(getattr(row_obj, attr)))
    
    return ";".join(s)
        
    
def make_zcdb_session():
    '''
        Returns a netutils.DBSession object with an open connection to the zero_click database installed on SRC_DB_HOST
    '''
    return netutils.DBSession(h=SRC_DB_HOST, u=SRC_DB_USER, p=SRC_DB_PASSWD, d=SRC_DB, driver=MYSQL_DRIVER)


def mark_configured(run_name, source_id, zc_session):
    ''' run_name (str):
            Name of the run to be marked as configured
            
        source_id (int):
            Source ID identifying the data source that produced the run
            
        zc_session (netutils.DBSession):
            Active session object with a connection to the Autonomics database
            
        Marks the specified run on the data source identified by source_id as configured in the runname_to_pid table.
    
    '''
    rn_2_pid = netutils.get_table_object("runname_to_pid", zc_session)
    zc_session.conn.execute(rn_2_pid.update().where(and_(rn_2_pid.c.run_name==run_name, 
                                                         rn_2_pid.c.source_id==source_id)).values(configured='Y'))
    

def mark_processed(run_name, source_id, zc_session):
    ''' run_name (str):
            The name of the run you want to mark as processed.
        
        source_id (int):
            The integer ID of the source this run was produced on
    
        zc_session (netutils.DBSession):
            Active session object with a connection to the Autonomics database
            
        Marks a run as being processed, by adding the run_name and source_id to the processed_runs table."
    '''
    processed_runs = netutils.get_table_object("processed_runs", zc_session)
    i = processed_runs.insert().values(source_id=source_id, run_name=run_name)
    i.execute()


def mark_job_paired(jid, session):
    ''' jid (int):
            Identifier for the job to mark as paired-end
            
        session (netutils.DBSession):
            Active session object with a connection to the Autonomics database
            
        Marks an individual job as having paired-end input data in the Autonomics database. 
    
    '''
    args = netutils.get_table_object('args', session)
    s = args.select(args.c.job_id==jid)
    res = session.conn.execute(s)
    row = res.fetchone()

    if(not settings.PAIRED_END_FLAG in row.pipeline_args):
        u = args.update()
        u = u.where(args.c.job_id==jid)
        u = u.values(pipeline_args=row.pipeline_args + " " + settings.PAIRED_END_FLAG)
        session.conn.execute(u) 


def mark_run_paired(run_name, source_id, session):
    ''' run_name (str): 
            Name of the run to be marked as having paired-end input
            
        source_id (int):
            Identifier for the data source that produced this run
            
        session (netutils.DBSession):
            Active session object with a connection to the Autonomics database
            
        Marks an entire run (project) as containing paired-end input data. Looks up each analysis task associated with the run that supports paired-end input and marks it as a paired-end analyis.
        
    '''
    trans = session.conn.begin()
    try:
        #update runname_to_pid
        r2pid = netutils.get_table_object('runname_to_pid', session)
        u = r2pid.update()
        u = u.where(and_(r2pid.c.run_name==run_name,
                r2pid.c.source_id==source_id))
        u = u.values(paired_end='Y')
        session.conn.execute(u)

        #get the project_id
        s = r2pid.select(and_(r2pid.c.run_name==run_name,
                              r2pid.c.source_id==source_id))
        res = session.conn.execute(s)
        row = res.fetchone()
        if(row is None):
            raise TypeError('None row received for project_id lookup')

        pid = row.project_id
    
        #update job args
        #get the jobs that support paired-end
        job_types = netutils.get_table_object('job_types', session)
        s = job_types.select(job_types.c.supports_paired_end=='Y')
        res = session.conn.execute(s)
        jobs_to_mark = [row.job_type for row in res.fetchall()]
        for job_type in jobs_to_mark:
            jid = netutils.get_jid(pid, job_type, session)
            mark_job_paired(jid, session)
        
        trans.commit()

    except:
        print("Error during mark_run_paired, rolling back.")
        trans.rollback()
    

def move_to_remote(local, remote, ssh, xfer_method='sftp'):
    '''
        local (str):
            Local path of file you wish to move.
        
        remote (str):
            Remote destination path of the file you wish to move.
        
        ssh (netutils.SSHConnection):
            A netutils.SSHConnection object connected to the host you wish to transfer data to.
            
        xfer_method (str):
            The transfer method to use when transfering local to remote. 
            Default: 'sftp'
            Supported Values: sftp, rsync, scp
    
        Moves a file with path 'local' to path 'remote', using the put method of ssh.
    '''
    
    remote_path, remote_fn = os.path.split(remote)
    # local_path, local_fn = os.path.split(local)
    remote_path += "/"
    # put the data on the zc server
    if(not netutils.remote_path_exists(ssh, remote_path)):
        ssh.execute("mkdir " + remote_path)
    
    if(xfer_method == 'sftp'):
        ssh.put(local, remote)
    else:
        remote = ssh.user + "@" + ssh.host + ":" + remote
        print(local)
        print(remote)
        netutils.scp(local, remote)


def move_ftp_data(host, source_id, zc_session, credentials, destination, src):
    '''
        host (str):
            String representation of the IP address for the FTP server containing the run data.
        
        source_id (int):
            Integer ID of the FTP server in the runname_to_pid table.
        
        credentials (settings.Credentials):
            A Credentials object produced by read_credentials.
        
        destination (str):
            The path to the ZIPline data directory (usually settings.proj_dir).
        
        src (str):
            The path to the folder on the FTP server containing the run data.
        
        Method to move data from ftp-server at 'host'.
    
        Supported filetypes are looked for in the src path at host, and moved one-by-one to destination. Skips files representing runs that have already been either configured or downloaded.
    
        After completing move, updates the run as 'downloaded' in the runname_to_pid table.
    '''
    rn_2_pid = netutils.get_table_object("runname_to_pid", zc_session)
    pn_mapping = netutils.get_table_object('pn_mapping', zc_session)
    zc_ssh = netutils.SSHConnection(settings.ZC_HOST,
                                    username=settings.ZC_USER, password=settings.ZC_PASSWD)
    loc = "local"
    if("localhost" in host):
        files = os.listdir(src)
                    
    else:
        loc = "remote"
        source_ssh = None
    
    for f in files:
        full_path = src + "/" + f
        full_paired = src + "/" + f + settings.PAIRED_EXTENSION
        dir, filename = os.path.split(full_path)
        base, ext = os.path.splitext(filename)
        if(ext in settings.SUPPORTED_FILETYPES):
            try:
                # add support for zipped files
                if(run_configured(base, source_id, zc_session) and not 
                   run_downloaded(base, source_id, zc_session) and not 
                   delete_pending(base, source_id, zc_session)):
                    #make sure the configuration is up-to-date on remote
                    if(settings.PUSH_CONFIGURATION_DATA):
                        push_configuration(zc_session)
                        '''
                        '''
                    proj_name = get_name_mapping(rn_2_pid, pn_mapping, base)
                    print(proj_name)
                    if(proj_name is None):
                        continue
                
                    # set up the remote directories and move the file
                    proj_dir = destination + proj_name + "/"
                    final_file = proj_dir + proj_name + ext
                    final_paired = proj_dir + proj_name + ext + \
                        settings.PAIRED_EXTENSION
                    dest_file = proj_dir + filename
                    dest_paired = proj_dir + filename + \
                        settings.PAIRED_EXTENSION
                
                    if(loc == "local"):
                        move_to_remote(full_path, dest_file, zc_ssh, 
                                       settings.ZC_XFER_METHOD)
                        if(os.path.isfile(full_paired)):
                            ensure_project_paired(base, source_id, zc_session)
                            move_to_remote(full_paired, dest_paired, 
                                           zc_ssh, settings.ZC_XFER_METHOD)
                            zc_ssh.execute("mv " + dest_paired + " " + \
                                           final_paired)
                    else:
                        move_remote_files(full_path, dest_file, source_ssh, 
                                          zc_ssh)
                
                    zc_ssh.execute("mv " + dest_file + " " + final_file)
                    # transfer an uploaded flag
                    cmd = "touch " + proj_dir + "SRC_UPLOADED"
                    zc_ssh.execute(cmd)
                    cmd = "chmod -R a+wrx " + proj_dir
                    zc_ssh.execute(cmd)
                    # mark run as uploaded
                    print("Updating project downloaded status.")
                    rn_2_pid.update().where(and_(rn_2_pid.c.run_name == base,
                                          rn_2_pid.c.source_id == source_id)).\
                                          values(downloaded='Y').execute()
                    print("Done.")
#                    os.remove(full_path)
#                    if(os.path.isfile(full_paired)):
#                        os.remove(full_paired)
                                          
            except Exception as e:
                raise
                sys.stderr.write("Error moving ftp file: " + f + "\n")
                sys.stderr.write(e.message + "\n") 


def move_it_data(host, source_id, zc_session, credentials, destination, src):
    '''
        host (str):
            String representation of the IP address for the FTP server 
            containing the run data.
        
        source_id (int):
            ID of the FTP server in the runname_to_pid table.
        
        credentials (settings.Credentials):
            A Credentials object produced by read_credentials.
        
        destination (str):
            The path to the ZIPline data directory (usually settings.proj_dir).
        
        src (str):
            The path to the folder on the FTP server containing the run data.
    
        Moves Ion Torrent sequencing runs on Ion Torrent Server located at 
        'host'. 
    
        The method queries the Ion Server for all finished runs. Runs that 
        have been configured but not downloaded are then moved one-by-one to 
        destination. 
    
        Expects that the Ion Torrent database is PostgreSQL, and that it 
        contains a database iondb with tables rundb_results and 
        rundb_experiment. Expects that rundb_results has a field fastqLink with
         the full path to the FASTQ-formatted output of the sequencing run.
    
        Updates runname_to_pid for the transferred run, setting downloaded to 
        'Y'.
    '''
    # get connection to destination and src
    zc_ssh = netutils.SSHConnection(settings.ZC_HOST, 
                                    username=settings.ZC_USER, 
                                    password=settings.ZC_PASSWD)
    ion_ssh = netutils.SSHConnection(host, username=credentials['user'], 
                                     password=credentials['passwd'])
    
    # connect to the ion torrent server to check for completed runs
    session = netutils.DBSession(h=host, d="iondb", u=credentials['dbuser'], 
                                 p=credentials['dbpasswd'], driver=POSTGRES_DRIVER)
    # get tables containing run data from iondb
    rundb_results = netutils.get_table_object("rundb_results", session)
    rundb_experiment = netutils.get_table_object("rundb_experiment", session)
    # get zero_click table that contains run->project mappings
    pn_mapping = netutils.get_table_object("pn_mapping", zc_session)
    run_2_pid = netutils.get_table_object("runname_to_pid", zc_session)
    results = select([rundb_results, rundb_experiment], 
                     and_(rundb_results.c.experiment_id == rundb_experiment.c.id, 
                          rundb_results.c.status == 'Completed'), 
                     use_labels=True).execute()
    for row in results.fetchall():
        
        #eventually move loop internals to new method, move_single_it
        
        fastq_path = row.rundb_results_fastqLink
        run_path = os.path.split(fastq_path)[0] + "/"
        orig_data_file = "rawlib.basecaller.bam"
        run_name = get_runname(fastq_path, "ion_torrent")
        renamed_data = run_name + ".bam"
        src_file = settings.IT_BASE_DIR + run_path + orig_data_file
        
        # check if this run has a project name in the system
        proj_name = get_name_mapping(run_2_pid, pn_mapping, run_name)
        if(proj_name is None or delete_pending(run_name, source_id, 
                                               zc_session) or 
           run_downloaded(run_name, source_id, zc_session)):
            #print("Run " + run_name + " doesn't exist in system, skipping.")
            continue
                
        proj_dir = destination + "/" + proj_name + "/"
        dest_file = proj_dir + renamed_data
        proj_fastq = proj_dir + proj_name + ".fastq"
        
        alternate = settings.IT_BASE_DIR + \
            run_path.replace("basecaller_results/", "") + "rawlib.bam"
        data_locations = [src_file, alternate]
        transfer_success = False
        for data_file in data_locations:
            try:
                netutils.get_file(data_file, renamed_data, 
                                  credentials, settings.IT_XFER_METHOD)
                transfer_success = True
                break
            
            except IOError as e:
                print("IOError: " + e.message + " | " + data_file)
        
        if(not transfer_success):
            print("Could not successfully transfer data for: " + \
                  run_name + ", skipping.")
            continue
        
        try:    
            # put the data on the zc server
            zc_ssh.execute("mkdir " + proj_dir)
            netutils.put_file(renamed_data, dest_file, settings.zc_cred, 
                              settings.ZC_XFER_METHOD)
            #zc_ssh.put(renamed_data, dest_file)
            #here we need to convert the BAM file to FASTQ 
            zc_ssh.execute("bamtools convert -format fastq -in " + dest_file +\
                            " -out " + proj_fastq)
            zc_ssh.execute("rm " + dest_file)
            # transfer an uploaded flag
            zc_ssh.execute("touch " + proj_dir + "SRC_UPLOADED")
            zc_ssh.execute("chmod -R a+wrx " + proj_dir)

            # mark run as uploaded
            zc_session.safe_exec(run_2_pid.update().\
                                 where(and_(run_2_pid.c.run_name == run_name,
                                          run_2_pid.c.source_id == source_id)).\
                                 values(downloaded='Y'))
            os.remove(renamed_data)

        except Exception as e:
            sys.stderr.write("Error during file transfer involving: " +\
                              src_file + "\n")
            sys.stderr.write(e.message + "\n")
        
    zc_ssh.close()
    ion_ssh.close()
    session.close()


def move_miseq_data(host, source_id, zc_session, credentials, destination, src):
    ''' host (str):
            Address for the MiSeq you wish to retrieve data from.
            
        source_id (int):
            Source ID for this MiSeq machine.
            
        zc_session (netutils.DBsession):
            A session object with an active connection to the autonomics 
            database.
            
        credentials (settings.Credentials):
            SSH Credentials for the MiSeq machine.
            
        destination (str):
            The destination path on the server specified by the SSH credentials.
            
        src (str):
            Source data path on the MiSeq sequencer.
            
        This function checks all runs associated with MiSeq sequencer 
        identified by source_id. Any finished run, as determined by the 
        presence of the 'FASTQ generation time:' line in the project's 
        AnalysisLog.txt file, is downloaded to the machine running 
        data_gremlin, before being pushed to the machine specified in 
        credentials.
        
    '''
    #select all of the configured runs that aren't downloaded
    ssh = netutils.SSHConnection(host, username=credentials.user, 
                                 password=credentials.passwd)
    zc_ssh = netutils.SSHConnection(settings.zc_cred.host, 
                                    username=settings.zc_cred.user, 
                                    password=settings.zc_cred.passwd)
    #iterate over the runs the system knows about
    rn_2_pid = netutils.get_table_object("runname_to_pid", zc_session)
    pn_mapping = netutils.get_table_object("pn_mapping", zc_session)
   
    sftp_path = fix_miseq_path_4sftp(src)
    log_name = "AnalysisLog.txt"
    fastq_path = "Data/Intensities/BaseCalls/"
   
    res = zc_session.conn.execute(select([pn_mapping.c.project_name, 
                                          rn_2_pid.c.run_name], 
                                         and_(rn_2_pid.c.source_id==source_id, 
                                              rn_2_pid.c.project_id==\
                                              pn_mapping.c.project_id, 
                                              rn_2_pid.c.downloaded=='N')))
    for row in res.fetchall():
        #get the log file for the run
        remote_project_dir = settings.proj_dir + row.project_name + "/"
        remote_ssh_dir = src + row.run_name + "/"
        remote_sftp_dir = sftp_path + row.run_name + "/"
        local_log_name = row.run_name + log_name
        netutils.get_file(remote_sftp_dir + log_name, local_log_name, 
                          credentials, xfer_mthd=settings.MISEQ_XFER_METHOD)
        
        seqs_generated = False
        #read the log file, seeing if the sequence generation is finished
        fh = open(local_log_name, 'r')
        for l in fh:
            if("FASTQ generation time" in l):
                seqs_generated = True
                
        fh.close()
        os.remove(local_log_name)
        
        ss_name = "SampleSheet.csv"
        local_sample_sheet = row.run_name + ss_name 
        #if the sequences are finished, do the move to the zero_click server    
        if(seqs_generated):
            #get SampleSheet.csv to read the sample name for this run
            netutils.get_file(remote_sftp_dir + ss_name, local_sample_sheet,
                               credentials, settings.MISEQ_XFER_METHOD)
            run_data = rundata_from_samplesheet(local_sample_sheet)
            os.remove(local_sample_sheet)
            
            #retrieve the sequence files for this run
            fastq_files = ssh.execute("ls " + remote_ssh_dir + fastq_path +\
                                       run_data.Sample_ID + "*.fastq.gz")
            for f in fastq_files:
                f = f.rstrip("\n")
                d, file_name = os.path.split(f)
                if("_R1_" in file_name):
                    dest_file = row.project_name + ".fastq.gz"
                elif("_R2_" in file_name):
                    dest_file = row.project_name + ".fastq.end2.gz"
                
                netutils.get_file(remote_sftp_dir + fastq_path + file_name, 
                                  dest_file, credentials, 
                                  settings.MISEQ_XFER_METHOD)
                zc_ssh.execute("mkdir " + remote_project_dir)
                netutils.put_file(dest_file, remote_project_dir + dest_file,
                                   settings.zc_cred, settings.ZC_XFER_METHOD)
                zc_ssh.execute("gunzip -f " + remote_project_dir + dest_file)
                os.remove(dest_file)
                
            zc_ssh.execute("touch " + remote_project_dir + "SRC_UPLOADED")
            zc_ssh.execute("chmod -R a+wrx " + remote_project_dir)
            zc_session.safe_exec(rn_2_pid.update().\
                                 where(and_(rn_2_pid.c.run_name==row.run_name,
                                            rn_2_pid.c.source_id==source_id)).\
                                 values(downloaded='Y'))
 
    
def move_remote_files(one, two, three, four):
    pass


def rundata_from_samplesheet(sheet_path):
    ''' sheet_path (str):
            The path to a SampleSheet.csv file.
            
        Returns a Record object with attributes and values matching those in 
        the [Data] section of the SampleSheet.csv
    '''
    fh = open(sheet_path, 'r')
    r = Record()
    WAITING_FOR_DATA = 0
    DATA_SEEN = 1
    FIELDS_PARSED = 2
    fields = []
    state = WAITING_FOR_DATA
    for l in fh:
        l = l.rstrip("\n")
        if(state == WAITING_FOR_DATA):
            if(l.startswith("[Data]")):
                state = DATA_SEEN
        elif(state == DATA_SEEN):
            els = l.split(",")
            for el in els:
                fields.append(el)
            state = FIELDS_PARSED
        elif(state == FIELDS_PARSED):
            els = l.split(",")
            for i, el in enumerate(els):
                setattr(r, fields[i], el)
            return r
        
    
def prune_projects(session):
    '''
    session (netutils.DBSession):
        A session object with an active connection to the Autonomics database. 
        
    Deletes run data from the autonomics MySQL database if they have been 
    marked for deletion via the web interface. 
    
    To mark for deletion, currently need to name the project with a prefix that
     includes either 'delete' or 'Delete'.
    '''

    pn_mapping = netutils.get_table_object("pn_mapping", session)
    rn_2_pid = netutils.get_table_object("runname_to_pid", session)
    data_sources = netutils.get_table_object("data_sources", session)
    
    # get the projects we want to prune
    results = pn_mapping.select(or_(pn_mapping.c.project_name.startswith("delete"),
                                    pn_mapping.c.project_name.startswith("Delete"),
                                    pn_mapping.c.project_name.contains("DUPLICATE"))).execute()
    
    
    for row in results.fetchall():
        print("removing project: " + str(row.project_id))
        pid = row.project_id
        # get all of the jobs for this project_id
        if(netutils.delete_project(pid, session)):
            print("Removal successful.")
        else:
            print("Removal failed")
            
    #check for proton runs that need removal
    res = select([pn_mapping.c.project_id, pn_mapping.c.project_name], 
                 and_(pn_mapping.c.project_id==rn_2_pid.c.project_id,
                      rn_2_pid.c.source_id==data_sources.c.source_id,
                      data_sources.c.source_type=='ion_proton')).execute()
    
    matcher = re.compile("_tn(_\d+)*$")
        
    for row in res.fetchall():
        if(not matcher.search(row.project_name) is None):
            if(netutils.delete_project(row.project_id, session)):
                print("Pruned project: " + row.project_name)               
                
                
def push_configuration(zc_session):
    ''' zc_session (netutils.DBSession):
            netutils.DBSession object connected to the ZIPline relational 
            database.
    
    This method pushes configuration data from the web-server and the compute 
    server, if these servers are hosted on different machines, or are behind 
    restrictive firewalls.
    
    Produces two tab-delimited files:
        <project_name>
        <project_name>_global
        
    run_name is constructed with:
        "\t".join([str(row.job_id), job.job_type, job.code, arg_str, "\n"])
        
    run_name_global:
        "\t".join([str(source_id), str(run_name), str(pid), downloaded, configured, paired])
        
    Where:
        arg_str = custom arguments to be set for the job configuration being pushed
        run_name = run_name specified by the sequencing or ftp server
        downloaded = whether or not this project is downloaded
        configured = whether or not this project is downloaded
        paired = whether or not this project includes paired-end data
    '''
    (rn_2_pid, pn_mapping, jn_mapping, config, ignore) = \
        get_configuration_tables(zc_session)
    args = netutils.get_table_object("args", zc_session)
    conn = netutils.SSHConnection(settings.ZC_HOST, 
                                  username=settings.ZC_USER,
                                  password=settings.ZC_PASSWD)
    results = select([pn_mapping, rn_2_pid], 
                     and_(pn_mapping.c.project_id == rn_2_pid.c.project_id,
                          rn_2_pid.c.downloaded == 'N', 
                          rn_2_pid.c.configured == 'Y'),
                     use_labels=True).execute()
                     
    for row in results.fetchall():
        pid = row.pn_mapping_project_id
        pn = row.pn_mapping_project_name
        source_id = row.runname_to_pid_source_id
        run_name = row.runname_to_pid_run_name
        paired = row.runname_to_pid_paired_end
        downloaded = row.runname_to_pid_downloaded
        configured = row.runname_to_pid_configured
        local_job_config = settings.CONFIG_PATH + pn
        local_project_config = settings.CONFIG_PATH + pn + "_global"
        fh_global = open(local_project_config, 'w')
        fh = open(local_job_config, 'w')
        # get each of the job_ids FROM jn_mapping, grab
        jobs = config.select(config.c.project_id == pid).execute()
        try:
            record = [str(source_id), str(run_name), str(pid), downloaded, 
                      configured, paired]
            fh_global.write("\t".join(record) + "\n")
            for job in jobs.fetchall():
                # get the job_id for each job
                jid_result = select([jn_mapping], 
                                    and_(jn_mapping.c.project_id == pid, 
                                         jn_mapping.c.job_type == job.job_type)\
                                    ).execute()
                row = jid_result.fetchone()
                # get the runname, source_id for this job
                # write the configuration of this job to file, put the args at
                # the end (if there are any)
                arg_results = args.select(args.c.job_id == row.job_id).execute()
                arg_row = arg_results.fetchone()
                arg_str = make_arg_string(arg_row)
                record = [str(row.job_id), job.job_type, job.code, 
                          arg_str, "\n"]
                fh.write("\t".join(record))
        except Exception as e:
            print("Error writing configuration files for " + run_name +\
                   " " + pn)
            print(" ".join([str(pid), str(source_id), str(run_name)]))
            print(e.message)

        fh_global.close()
        fh.close()
        
        # move the configuration files to the ZC pipeline
        config_path = settings.proj_dir + settings.CONFIG_BASE
        job_dest = config_path + pn
        project_dest = config_path + pn + "_global"
        try:
            conn.put(local_job_config, job_dest)
            conn.execute("chmod a+rwx " + job_dest)
            conn.put(local_project_config, project_dest)
            conn.execute("chmod a+rwx " + project_dest)
            os.remove(local_project_config)
            os.remove(local_job_config)
        
        except Exception as e:
            print("Error pushing configuration files for " + run_name + " "\
                   + pn)
            print(" ".join([str(pid), str(source_id), str(run_name)]))
            print("Job configuration exists: " + \
                  str(os.path.exists(local_job_config)))
            print("Project configuration exists: " + \
                  str(os.path.exists(local_project_config)))
            print(e.message)
           
    conn.close()
    

def run_configured(run_name, source_id, zc_session):
    '''
        run_name (str):
            String name of the run.
        
        source_id (int):
            Integer ID of the source that produced this run.
    
        Returns True if the run is configured in runname_to_pid, False 
        otherwise.
    '''
    rn_2_pid = netutils.get_table_object("runname_to_pid", zc_session)
    results = rn_2_pid.select(and_(rn_2_pid.c.run_name == run_name,
                                  rn_2_pid.c.source_id == source_id,
                                  rn_2_pid.c.configured == 'Y')).execute()
    configured = False
    for row in results.fetchall():
        configured = True
        
    return configured


def run_downloaded(run_name, source_id, zc_session):
    '''
        run_name (str):
            String name of the run.
        
        source_id (int):
            Integer ID of the source that produced this run.
        
        Returns True if the run is marked as downloaded in runname_to_pid, 
        False otherwise.
    '''
    rn_2_pid = netutils.get_table_object("runname_to_pid", zc_session)
    results = rn_2_pid.select(and_(rn_2_pid.c.run_name == run_name,
                                  rn_2_pid.c.source_id == source_id,
                                  rn_2_pid.c.downloaded == 'Y')).execute()
    
    for row in results.fetchall():
        return True
    
    return False

    
def run_processed(run_name, source_id, zc_session):
    '''
        run_name:
            String name of the run.
        
        source_id:
            Integer ID of the source that produced this run.
    
        Returns True if run_name, source_id is stored in the table 
        processed_runs.
    '''
    # check if this run has been processed
    processed_runs = netutils.get_table_object("processed_runs", zc_session)
    res = processed_runs.select(and_(processed_runs.c.run_name == run_name, 
                                     processed_runs.c.source_id == source_id)).\
                                     execute()
    row = res.fetchone()
    if(row is None):
        return False
    
    return True


def source_active(source_id):
    return True


def update_preassemble_config(pid, assembler, zc_session):
    '''
        assembler (str):
            The assembler to use. Supported values are Trinity and MIRA.
        
        pid (int):
            The project_id that will be used to lookup the jobs to configure.
        
        assem_config (Table):
            SQLAlchemy Table object reflection for the table containing 
            assembler-specific job configurations.
        
        config (Table):
            SQLAlchemy Table object reflection for the table containing job 
            configurations for all projects.
    
        Based on the assembler (Trinity, MIRA), configures the appropriate 
        pre-assembly steps for project with ID=pid.
    
        Current configurations:
            Trinity - quality_trim, adapter_trim, read_normalization
            MIRA - quality_trim, adapter_trim
        
        
    '''
    
    assem_config = netutils.get_table_object("assembler_config", zc_session)
    config = netutils.get_table_object("configuration", zc_session)
    # update configuration for this assembler
    res = assem_config.select(assem_config.c.assembler==assembler).execute()
    for row in res.fetchall():
        u = config.update().where(and_(config.c.project_id == pid,
                                       config.c.job_type == row.job_type)).\
                                       values(code=row.code)
        u.execute()


def update_project_assembler(pid, assembler, zc_session):
    ''' pid (int):
            Project ID off the project to associate with assembler.
        assembler (str):
            Identifier of the assembler to assign to this project's jobs.
            
        zc_session (netuils.DBSession):
            Session object with active connection to the autonomics database.
            
        Updates the pipeline_args of any 'assemble' jobs associated with the 
        supplied project ID, setting --assembler = assembler.
    '''
    default_args = netutils.get_table_object("default_args", zc_session)
    args = netutils.get_table_object("args", zc_session)
    jn_mapping = netutils.get_table_object("jn_mapping", zc_session)
    
    jid_res = jn_mapping.select(and_(jn_mapping.c.project_id==pid, 
                                     jn_mapping.c.job_type=='assemble')).\
                                     execute()
    row = jid_res.fetchone()
    jid = row.job_id
    
    assemble_replacer = re.compile("--assembler\s\w+")
    
    # set the appropriate pipeline_args for the assembly job
    res = default_args.select(default_args.c.job_type == 'assemble').execute()
    row = res.fetchone()
    pipe_args = row.pipeline_args
    pipe_args = assemble_replacer.sub("--assembler " + assembler, pipe_args)
        
    i = args.insert().values(job_id=jid, loc=row.loc,
                                 executable=row.executable, 
                                 process_args=row.process_args, 
                                 pipeline_args=pipe_args,
                                 resources=row.resources, 
                                 priority=row.priority)
    i.execute()
            

def write_configuration(run_name, source_id, zc_session, 
                        paired=False, project_name=None):
    '''
        run_name (str):
            Same as above.
        
        rn_2_pid (Table):
            Same as above.
        
        pn_mapping (Table):
            SQLAlchemy Table object for the 'pn_mapping' table.
        
        jn_mapping (Table):
            SQLAlchemy Table object for the 'jn_mapping' table.

        config (Table):
            SQLAlchemy Table object for the 'configuration' table.
        
        default (Table):
            SQLAlchemy Table object for the 'default_configuration' table.
    
        source_id (int):
            Same as above.
        
        paired (Boolean):
            Boolean, whether or not this run has paired-end input data.
    
        Writes the default configuration for the run given in run_name.
    
        Default configuration for un-configured jobs can be found in the 
        default_configuration table of the MySQL database.
    
        Performs the following steps for the run:
            - adds an entry to pn_mapping
            - adds an entry for each default job in configuration
            - creates job_ids for each default-ly configured job
            - sets paired-end flags in pipeline_args if this run was paired-end
        
        To-do: add support for writing configuration for pre-assembly jobs 
        based on assembler, currently handled within the individual configuration handlers which may call this method.
    '''
    (rn_2_pid, pn_mapping, jn_mapping, config, default) = \
        get_configuration_tables(zc_session)
    # insert an entry into pn_mapping
    if(project_name is None):
        project_name = run_name
    s = pn_mapping.select(pn_mapping.c.project_name == project_name)
    i = pn_mapping.insert().values(project_name=project_name)
    pid = netutils.check_and_insert(s, i)
    
    # if the project name is taken, try to append the source_id and submit again
    # if(pid is None): 
    #    base_name = base_name + "_" + str(source_id)
    #    s = pn_mapping.select(pn_mapping.c.project_name==base_name)
    #    i = pn_mapping.insert().values(project_name=base_name)
    #    pid = netutils.check_and_insert(s, i)
        
    # if we still didn't have any luck, return failure
    if(pid is None):
        return False

    # set the appropriate paired-end flag
    if(paired):
        paired_end = 'Y'
    else:
        paired_end = 'N'

    # insert an entry in rn_2_pn
    s = rn_2_pid.select(and_(rn_2_pid.c.run_name == run_name, 
                             rn_2_pid.c.source_id == source_id))
    i = rn_2_pid.insert().values(source_id=source_id, run_name=run_name, 
                                 project_id=pid, paired_end=paired_end)
    netutils.check_and_insert(s, i)

    # insert entries into jn_mapping, configuration
    job_types = netutils.get_table_object("job_types", zc_session)
    default_args = netutils.get_table_object("default_args", zc_session)
    args = netutils.get_table_object("args", zc_session)
    results = select([default, job_types, default_args],
                    and_(default.c.job_type == job_types.c.type,
                        default_args.c.job_type == job_types.c.type), 
                     use_labels=True).execute()
    
    for row in results.fetchall():
        jt = row.default_configuration_job_type
        jcode = row.default_configuration_code
        s = config.select(and_(config.c.project_id==pid, config.c.job_type==jt))
        i = config.insert().values(project_id=pid, job_type=jt, code=jcode)
        netutils.check_and_insert(s, i)

        s = jn_mapping.select(and_(jn_mapping.c.project_id == pid, 
                                   jn_mapping.c.job_type == jt))
        i = jn_mapping.insert().values(project_id=pid, job_type=jt)
        jid = netutils.check_and_insert(s, i)
        
        # check if this project is paired-end, and if this job supports paired-end input
        if(paired and row.job_types_supports_paired_end=='Y'):
            # insert paired-end args for this job
            i = args.insert().values(job_id=jid, 
                                     executable=row.default_args_executable, 
                                     loc=row.default_args_loc, 
                                     resources=row.default_args_resources, 
                                     process_args=row.default_args_process_args,
                                     pipeline_args="--paired_end " + \
                                     row.default_args_pipeline_args,
                                     priority=row.default_args_priority)          
            i.execute()
            
    #update the pre-assembly config, if possible
    #check for an assembly directive in the run_name
    run_base, assembler = os.path.splitext(run_name)
    assembler = assembler.lstrip(".")
    if(assembler != ""):
        update_project_assembler(pid, assembler, zc_session)
        update_preassemble_config(pid, assembler, zc_session)
            
    return True


########################### DB CONFIGURATION ###############################
SRC_DB_HOST = settings.db_cred.host
SRC_DB = "zero_click"
SRC_DB_USER = settings.db_cred.user
SRC_DB_PASSWD = settings.db_cred.passwd
if settings.plw:
    SRC_DB = "plw_db2"

MYSQL_DRIVER = "mysql+mysqldb"
POSTGRES_DRIVER = "postgresql+pg8000"

config_runs = {"ftp": init_ftp_runs, "ion_torrent": init_it_runs, 
               "ion_proton": init_it_runs, "miseq": init_miseq_runs}
file_movers = {"ftp": move_ftp_data, "ion_torrent": move_it_data, 
               "ion_proton": move_it_data, "miseq": move_miseq_data}
data_sources = ["ftp", "ion_torrent", "ion_proton", "miseq"]

def main(args):

    DEST_DIR = args.dest_dir
    SLEEP_INTERVAL = int(args.sleep_interval)
    while(True):
        # connect to the database to get a current copy of the data sources we're using
        session = netutils.DBSession(SRC_DB_HOST, SRC_DB, SRC_DB_USER, 
                                     SRC_DB_PASSWD, driver=MYSQL_DRIVER)
        
        # remove projects marked for deletion
        
        src_tbl = netutils.get_table_object("data_sources", session)
        for src_type in data_sources:
            results = src_tbl.select(and_(src_tbl.c.active == 1, 
                                          src_tbl.c.source_type == src_type))\
                                          .execute()
            for row in results.fetchall():
                print("Processing: " + row.source_name)
                # check if this data source is active
                if(source_active(row.source_id)):
                    cred = Credentials(from_file=settings.CRED_PATH + \
                                       str(row.host))
                    # write this configuration data to the ZC server
                    try:
                        print("Configuring runs on: " + row.source_name)
                        config_runs[src_type](row.host, row.source_id, 
                                              session, cred, path=row.directory)
                    except Exception as e:
                        print("Error configuring runs for source: " +\
                               row.source_name)
                        print(e)

                    # push default configuration for all runs to the db
                    if(settings.PUSH_CONFIGURATION_DATA):
                        prune_projects(session)
                        try:
                            print("Pushing configuration data.")
                            push_configuration(session)
                        except Exception as e:
                            print("Error pushing configurations for " +\
                                   row.source_name)
                            print(e)
                    # move completed runs to ZC server
                    try:
                        print("Moving run data for: " + row.source_name)
                        file_movers[src_type](row.host, row.source_id, session,
                                               cred, DEST_DIR, row.directory)
                    except Exception as e:
                        print("Error moving data files for " + row.source_name)
                        print(e)
                   

        # for each datasource, check if there's new data
        # move the data to the dest, renaming if required (check file_renames table)

        session.close()
        print("Sleeping for: " + str(SLEEP_INTERVAL))
        time.sleep(SLEEP_INTERVAL)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--dest-dir', dest='dest_dir', 
                   default=settings.proj_dir, 
                   help='Destination on the remote server where the gremlin \
                   should put files.')
    p.add_argument('--sleep-interval', dest='sleep_interval', 
                   default=60, help="How long this gremlin should sleep \
                   between checks.")
    main(p.parse_args())

