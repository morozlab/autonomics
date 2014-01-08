'''
Created on May 12, 2012

@author: Mathew Citarella


This module monitors the directory proj_config/ located in the system's main data directory for configuration files pushed there by data_gremlin.py. This module is necessary to maintain that the ZIPline server's MySQL database contains every project configured on the web server.

This module is run as a daemon and performs the following tasks:

1) Get a list of configuration files stored in settings.proj_dir/proj_config/ - configuration files for a project are named [project_name, project_name_global]; see the documentation for push_configuration() in the data_gremlin module for the format of these files.
2) Check for, and open, the global configuration file containing project and run information. If the project exists in the system, update any project-wide settings. Otherwise create a new project in the system and add the system-wide settings. System-wide settings include whether or not the project includes paired-end data, was configured on the web-server, and if it has been downloaded yet.
3) Open the job-specific configuration file. This file contains a line for each configured job for the project. If the jobs exist in the MySQL database on the ZIPline server, update job settings. Otherwise add the jobs and set the job settings.

Note: Please see the ZIPline database documentation for the runname_to_pid, pn_mapping, jn_mapping, args, and configuration tables for detailed descriptions of the various project and job settings.

'''

from autonomics import settings, netutils
from sqlalchemy.sql import select, insert, update, and_
import sqlalchemy.exc
import os, time

CHECK_PATH = settings.proj_dir + settings.CONFIG_BASE

def check_run_exists(source_id, run_name, session):
    '''
    source_id (int):
        Unique integer ID assigned to the source in the runname_to_pid table

    run_name (str):
        Name for this run as a string. source_id, run_name form the primary key of the runname_to_pid table

    session (netutils.DBSession):
        A netutils.DBSession object, which can be used to reflect database tables via netutils.get_table_object(). Can also directly execute queries through the conn attribute.

    Uses runname and source_id to check if the specified run exists in the database currently connected to by session.

    Returns True if the run exists, False otherwise.
    '''

    rn_2_pid = netutils.get_table_object("runname_to_pid", session)
    results = rn_2_pid.select(and_(rn_2_pid.c.source_id==source_id, rn_2_pid.c.run_name==run_name)).execute()
    row = results.fetchone()
    if(row is None):
        return False
    return True


def get_pid(source_id, run_name, session):
    '''
    source_id (int):
        Unique inter ID assigned to the source in runname_to_pid table

    run_name (str):
        Name of the run for which you'd like to get a project_id

    session (netutils.DBSession):
        An open session with the autonomics database

    Uses source_id and run_name to fetch and return the project_id assigned to this run.

    '''

    rn_2_pid = netutils.get_table_object("runname_to_pid", session)
    results = rn_2_pid.select(and_(rn_2_pid.c.run_name==run_name, rn_2_pid.c.source_id==source_id)).execute()
    row = results.fetchone()
    return row.project_id


def get_project_name(source_id, run_name, session):
    '''
    source_id (int):
        Unique source_id identifying the data source that produced this run

    run_name (str):
        Name of the run for which a project name should be returned

    session (netutils.DBSession):
        An open session with the autonomics database
    Uses source_id and run_name to fetch and return the project_name for the run.
    '''
    rn_2_pid = netutils.get_table_object("runname_to_pid", session)
    pn_mapping = netutils.get_table_object("pn_mapping", session)
    results = select([pn_mapping, rn_2_pid], and_(rn_2_pid.c.project_id==pn_mapping.c.project_id, rn_2_pid.c.run_name==run_name, rn_2_pid.c.source_id==source_id)).execute()
    row = results.fetchone()
    return row.project_name


def update_project_name(old_name, new_name, session):
    '''
    old_name (str):
        Name of the project you wish to change

    new_name (str):
        New name for the project

    session (netutils.DBSession):
        An open session with the autonomics database

    Changes old_name to new_name using the connection to the database specified by session.
    '''
    pn_mapping = netutils.get_table_object("pn_mapping", session)
    u = pn_mapping.update().where(pn_mapping.c.project_name==old_name).values(project_name=new_name)
    u.execute()


def main():

    session = netutils.make_db_session()
    pn_mapping = netutils.get_table_object("pn_mapping", session)
    jn_mapping = netutils.get_table_object("jn_mapping", session)
    config = netutils.get_table_object("configuration", session)
    rn_2_pid = netutils.get_table_object("runname_to_pid", session)
    args = netutils.get_table_object("args", session)

    proj_attrs = ["source_id", "run_name", "pid", "downloaded", "configured", "paired"]
    arg_attrs = ["executable", "loc", "process_args", "pipeline_args", "resources", "priority"]

    while(True):

        project_configs = {}
        for f in os.listdir(CHECK_PATH):
            if(f == "." or f == ".."):
                continue

            if("_global" in f):
                pn = f.replace("_global", "")
                project_configs[pn] = {}
                fh = open(CHECK_PATH + f, 'r')
                for line in fh.readlines():
                    #record = [str(source_id), str(run_name), str(pid), downloaded, configured, paired]
                    line = line.rstrip("\n")
                    el = line.split("\t")

                    (source_id, run_name, discard, downloaded, configured, paired) = el[:6]
                    for i, attr in enumerate(proj_attrs):
                        project_configs[pn][proj_attrs[i]] = el[i]

                pid = None
                if(check_run_exists(source_id, run_name, session)):
                    pid = get_pid(source_id, run_name, session)
                    #update the details for this project
                    update_project_name(get_project_name(source_id, run_name, session), pn, session)
                    u = rn_2_pid.update().where(rn_2_pid.c.project_id==pid).values(downloaded=downloaded,
                                                                                   configured=configured,
                        paired_end=paired)
                else:
                    #need to insert this new project
                    s = pn_mapping.select(pn_mapping.c.project_name==pn)
                    i = pn_mapping.insert().values(project_name=pn)
                    pid = netutils.check_and_insert(s, i)

                    s = rn_2_pid.select(and_(rn_2_pid.c.run_name==run_name, rn_2_pid.c.source_id==source_id))
                    i = rn_2_pid.insert().values(run_name=run_name, source_id=source_id,
                                                 project_id=pid, downloaded=downloaded, configured=configured,
                        paired_end=paired)
                    netutils.check_and_insert(s, i)

                fh.close()
                os.remove(CHECK_PATH + f)

        for pn, attrs in project_configs.items():
            try:
                pid = get_pid(attrs['source_id'], attrs['run_name'], session)
                fh = open(CHECK_PATH + pn, 'r')
                for line in fh.readlines():
                    line = line.rstrip("\n")
                    el = line.split("\t")
                    # record = [source_id, run_name, str(pid), str(jid), job.job_type, job.code, "\n"]
                    discard, jtype, jcode, arg_str = el[:4]

                    #check if this job is in the jn_mapping table
                    i = jn_mapping.insert().values(project_id=pid, job_type=jtype)
                    s = jn_mapping.select(and_(jn_mapping.c.project_id==pid, jn_mapping.c.job_type==jtype))
                    jid = netutils.check_and_insert(s, i)
                    if(jid is None):
                        jid = s.execute().fetchone().job_id

                    #update the job configuration for this job
                    netutils.update_on_exists(config.insert().values(project_id=pid, job_type=jtype, code=jcode), config.update().where(and_(config.c.project_id==pid, config.c.job_type==jtype)).values(code=jcode))

                    #update the args for this job
                    if(not arg_str == ''):
                        arg_list = arg_str.split(";")
                        i = args.insert().values(job_id=jid, executable=arg_list[0], loc=arg_list[1], process_args=arg_list[2],
                                             pipeline_args=arg_list[3], resources=arg_list[4], priority=arg_list[5])
                        u = args.update().where(args.c.job_id==jid).values(executable=arg_list[0], loc=arg_list[1],
                                                                       process_args=arg_list[2], pipeline_args=arg_list[3],
                                                                       resources=arg_list[4], priority=arg_list[5])
                        netutils.update_on_exists(i, u)

                fh.close()
                os.remove(CHECK_PATH + pn)
            except Exception as e:
                print(e.message)
        time.sleep(settings.PROCESS_MONITOR_INTERVAL)

    session.close()


if __name__ == '__main__':
    main()
