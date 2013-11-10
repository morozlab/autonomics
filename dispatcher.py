'''
Peter L. Williams
November 29, 2012 .. April 28, 2013

This module runs as a daemon and is responsible for monitoring new project
data uploaded to the autonomics server by the data_gremlin module.  New
projects found are added to the pn_mapping, run_stats and init_projects
tables; all required jobs for a project are added to the jn_mapping &
jid_dependencies tables.  In addition, projects are checked that all
of its jobs are completed

0) Loops continuously, performing the following tasks, sleeping for
settings.DISPATCHER_SLEEP_INTERVAL between loops.

1) All jobs in jn_mapping table that (a) are not in queue table, and
(b) are not marked started, and (c) have all dependencies satisfied, are
inserted into queue table and jn_mapping table set to queued = 'Y'

2) Iterates over all directories found in settings.proj_dir, treating each
directory as a project_name.

3) For each project_name, checks if the directory for that project
(settings.proj_dir + project_name) contains SRC_UPLOADED, which marks that
data_gremlin.py has finished its upload. If so, does the following for
each project.

4) Checks pn_mapping table to see if project_name has been assigned a
project_id. (project_id needed to access all other tables).  If project
in pn_mapping table either: (1) project already initiated or (2) project
configuration was submitted via web interface which makes entries in
pn_mapping, jn_mapping, config & args tables (we call these custom
projects as opposed to those that use the default job configuration;
configuration refers to what jobs to run for a project and which args
to use for a job.)

5) If project in pn_mapping table and not in init_projects table then its
a custom project and we get a list of the jobs to run from the
configuration table.  Dependencies for each job are computed from the
dependency table and job_ids for each dependency are entered into the
jid_dependencies table.  The project_id is inserted into the init_projects
table.

6) If project is not a custom project we do the same procedure as in the
above step except we use the delfault_configuration table. In addition,
the project is inserted into the pn_mapping table and the jobs are
inserted into the jn_mapping table.

7) Finally, projects that have been initiated but not completed are
checked to see if all jobs have finished, if so project is added to
completed projects table.  If project has an 'upload' job_type entry
in the jn_mapping table, this entry is set to started and finished and
the time stamps are set.
'''

from autonomics import settings, netutils
from autonomics.file_io import FileExtensions
from autonomics.queue import queue_all
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

#   Deal with database:
    session = netutils.DBSession("localhost", settings.ZC_DB_NAME,
                                settings.db_cred.user, settings.db_cred.passwd)
    completed_projects_table = netutils.get_table_object("completed_projects",
                                                                       session)
    run_stats = netutils.get_table_object('run_stats', session)
    init_projects_table = netutils.get_table_object("init_projects", session)
    configuration_table = netutils.get_table_object("configuration", session)
    dependency_table = netutils.get_table_object("dependency", session)
    default_configuration_table = netutils.get_table_object(
                                              "default_configuration", session)
    pn_mapping_table = netutils.get_table_object("pn_mapping", session)
    jn_mapping_table = netutils.get_table_object("jn_mapping", session)
    quenew_table = netutils.get_table_object("quenew", session)
    jid_dependency_table = netutils.get_table_object("jid_dependency", session)

    while True:

#       ENTER NEW PROJECTS INTO QUEUE
        queue_all(session)

#       find names of projects in pipeline_path
        pipeline_path = settings.proj_dir
        project_name = ""
        dir_names = [ f for f in os.listdir(pipeline_path)
                                              if isdir(join(pipeline_path,f)) ]

        projects_to_run = []
        for dname in dir_names:
            dname_flag = dname + "/SRC_UPLOADED"
            dname_flag_with_path = pipeline_path + '/' + dname_flag
            if os.path.exists(dname_flag_with_path):
                project_name = dname
            else:
                continue
            projects_to_run.append(project_name)

#       code from here to "END ENTER NEW PROJECTS" does not get executed 
#       unless there is a SRC_UPLOADED flag in proj dir; only if using data gremlin
#       will there be a SRC_UPLOADED flag, so if not using data gremlin, then
#       projects_to_run is empty so entire block of code skipped.

        for project_name in projects_to_run:
            # look in pn_mapping table to see if project_name has a
            # project_id assigned (need pid so can look in init_projects
            # table to see if job already initiated.)
            s = pn_mapping_table.select(
                  pn_mapping_table.c.project_id).where(
                                 pn_mapping_table.c.project_name==project_name)
            results = s.execute()
            row = results.fetchone()
            job_types = []
            if row: # job is in pn_mapping table so either: (1) job already
                    # initiated or (2) job was submitted via web interface
                    # which makes entries in pn_mapping, jn_mapping, config,
                    # & args tables.
                project_id = row.project_id
                s = init_projects_table.select(
                    init_projects_table.c.ts).where(
                                  init_projects_table.c.project_id==project_id)
                results = s.execute()
                if not results.fetchone(): # proj not in init_projects table
                    # so job must have been submitted through web interface.
                    print "Adding ", project_name, " submitted via web"
                    i = init_projects_table.insert()
                    i.execute(project_id=project_id)
                    i = run_stats.insert()
                    try:
                        i.execute(project_id=project_id)
                    except Exception as e:
                        print("Warning: " + str(project_id) +
                                         " already exists in run_stats table.")
                    s = configuration_table.select().where(
                                  configuration_table.c.project_id==project_id)
                    results = s.execute()
                    for row in results:
                        if row.code != '-':
                            job_type = row.job_type
                            job_types.append(job_type)
                    for jt in job_types:
                        s = jn_mapping_table.select(and_(
                                               jn_mapping_table.c.job_type==jt,
                                    jn_mapping_table.c.project_id==project_id))
                        results = s.execute()
                        jid = results.fetchone().job_id
                        s = dependency_table.select().where(
                                               dependency_table.c.job_type==jt)
                        results = s.execute()
                        for row in results:
                            job_type_dep = row.depends_on
                            s = jn_mapping_table.select(and_(
                                     jn_mapping_table.c.job_type==job_type_dep,
                                    jn_mapping_table.c.project_id==project_id))
                            results = s.execute()
                            row2 = results.fetchone()
                            if row2:
                                job_id_dep = row2.job_id
                                i = jid_dependency_table.insert()
                                i.execute(job_id=jid,depends_on=job_id_dep)
            else: # <project_name> not in pn_mapping, so ITS A DEFAULT JOB

                # insert project name into pn_mapping table and get back a project_id
                i = pn_mapping_table.insert()
                i.execute(project_name=project_name)
                s = pn_mapping_table.select(
                       pn_mapping_table.c.project_id).where(
                                 pn_mapping_table.c.project_name==project_name)
                results = s.execute()
                project_id = results.fetchone().project_id
                print "Adding project: ", project_name, " with default jobs"

                # insert project id into run_stats table
                i = run_stats.insert()
                try:
                    i.execute(project_id=project_id)
                except Exception as e:
                    print("Warning: " + str(projectID) +
                                         " already exists in run_stats table.")
     
                # insert proj into init_projects table
                i = init_projects_table.insert()
                i.execute(project_id=project_id)

                # find all different job_types to run for a default job
                job_types = []
                s = default_configuration_table.select()
                results = s.execute()
                for row in results:
                    if row.code == '+':
                        job_type = row.job_type
                        job_types.append(job_type)
                        job_name = project_name + '_' + job_type
                        i = jn_mapping_table.insert()
                        i.execute(job_type=job_type,project_id=project_id, job_name=job_name)

                # for each job_type, find and insert dependencies in jid_dependency_table
                # and insert jobs into queue & set queued = 'Y' in jn_mapping table.
                for jt in job_types:
                    s = jn_mapping_table.select(and_(
                                               jn_mapping_table.c.job_type==jt,
                                    jn_mapping_table.c.project_id==project_id))
                    results = s.execute()
                    jid = results.fetchone().job_id
                    s = dependency_table.select().where(
                                               dependency_table.c.job_type==jt)
                    results = s.execute()
                    for row in results:
                        job_type_dep = row.depends_on
                        s = jn_mapping_table.select(and_(
                                     jn_mapping_table.c.job_type==job_type_dep,
                                    jn_mapping_table.c.project_id==project_id))
                        results = s.execute()
                        row2 = results.fetchone()
                        if row2:
                            job_id_dep = row2.job_id
                            i = jid_dependency_table.insert()
                            i.execute(job_id=jid,depends_on=job_id_dep)
                    # insert each job_type into queue and set queued = 'Y' in jn_mapping table
                    print "i = quenew_table.insert().values(project_id = ",project_id," job_id= ",jid," job_type= ",jt,").execute()"
                    print "u = jn_mapping.update().where(jn_mapping.c.job_id== ",jid," ).values(queued=\'Y\').execute()"
                    i = quenew_table.insert().values(project_id = project_id,job_id=jid,job_type=jt).execute()
                    u = jn_mapping_table.update().where(jn_mapping_table.c.job_id==jid).values(queued='Y').execute()

#       END ENTER NEW PROJECTS

#       since init_projects table not being used, next block of code not used.
#       NOW CHECK FOR COMPLETED PROJECTS & UPLOADS
        results=session.conn.execute("select project_id from init_projects\
          where project_id not in (select project_id from completed_projects)")
        for row in results.fetchall():
            project_id = row.project_id
            s = jn_mapping_table.select(jn_mapping_table.c.project_id).where(
                                     jn_mapping_table.c.project_id==project_id)
            results = s.execute();
            if results.fetchone():
                results = session.conn.execute(
                    "select job_type,job_id from jn_mapping where(\
                    (project_id = ?) and (finished = 'N'))",(project_id,))
                found_upload_job = 0
                num_rows_not_finished = 0
                upload_job_id = 0
                for row in results.fetchall(): # all rows where finished='N'
                    num_rows_not_finished += 1
                    if row.job_type == 'upload':
                        upload_job_id = row.job_id
                        found_upload_job = 1
                if (found_upload_job and num_rows_not_finished==1) or \
                                                    num_rows_not_finished == 0:
                    i = completed_projects_table.insert()
                    i.execute(project_id=project_id)
                    u = jn_mapping_table.update().where(
                         jn_mapping_table.c.job_id==upload_job_id).values(
                        s_ts=functions.current_timestamp(), started='Y').execute()
                    u = jn_mapping_table.update().where(
                         jn_mapping_table.c.job_id==upload_job_id).values(
                    f_ts=functions.current_timestamp(), finished='Y').execute()
                    u = quenew_table.delete().where(
                               quenew_table.c.project_id==project_id).execute()

        sys.stdout.write('.')
        time.sleep(settings.DISPATCHER_SLEEP_INTERVAL)
        sys.stdout.flush()
        time.sleep(settings.DISPATCHER_SLEEP_INTERVAL)

if __name__ == '__main__':
    main()
