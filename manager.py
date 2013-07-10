'''
Name:        manager.py
Purpose:     Manage jobs submitted to the autonomics annotation system

Author:      Mathew Citarella

Created:     01/10/2012
Copyright:   (c) Mathew Citarella 2012
'''

from sqlalchemy import exc
from sqlalchemy.sql import functions, select, and_
from autonomics import netutils, settings
from autonomics.jobs import BlastJob, BlastJobNR, BlastJobSwissprot, BlatJob
from autonomics.jobs import BlastAssociationJob, PfamJob, AssemblyJob
from autonomics.jobs import JobState, Locations, Resources, AdapterTrimJob
from autonomics.jobs import ReadNormJob, PantherJob, QuantificationJob
from autonomics.jobs import QualityTrimJob, UploadJob
from autonomics.queue import remove_from_queue
import argparse
import datetime
import gc
import time
import sys


def make_upload_job(q_r,args_r):
    ''' 
        q_r: 
            A sqlalchemy row from the queue table, representing an upload job.
        args_r: 
            A sqlalchemy row from the args table for this job.

        Returns an UploadJob object.
    '''
    return UploadJob(q_r.project_id, q_r.job_id, q_r.job_type)


def make_hpc_job(q_r,args_r):
    '''
        q_r: 
            An sqlalchemy row from the queue table for any of the HPCJob  
            subclasses.
    
        args_r: An sqlalchemy row from the args table for this job.

    
        Returns either a BlastJobNR, BlastJobSwissprot, BlastJob, PfamJob, 
        BlatJob, or QuantificationJob object.
    '''
    jobs = {"blast_nr": BlastJobNR, "blast_swissprot": BlastJobSwissprot, 
            "blast": BlastJob, "pfam": PfamJob, "blat": BlatJob, 
            "quantification": QuantificationJob}
    return jobs[q_r.job_type](q_r.project_id, q_r.job_id, q_r.job_type, 
                args_r.executable, args_r.resources, args_r.pipeline_args, 
                args_r.process_args)


def make_blast_association_job(q_r,args_r):
    '''
    q_r: An sqlalchemy row from the queue table representing a 
         BlastAssociation job.
    args_r: An sqlalchemy row from the args table for this job.

    Returns a BlastAssocationJob object.
    '''
    return BlastAssociationJob(q_r.project_id, q_r.job_id, q_r.job_type, 
                               args_r.pipeline_args, args_r.process_args)


def make_assembly_job(q_r,args_r):
    '''
    q_r: 
        An sqlalchemy row from the queue table representing an Assembly job.
    args_r: 
        An sqlalchemy row from the args table for this job.

    Returns an AssemblyJob object.
    '''
    return AssemblyJob(q_r.project_id, q_r.job_id, q_r.job_type, 
                       args_r.executable, args_r.resources, 
                       args_r.pipeline_args, args_r.process_args)


def make_preassembly_job(q_r,args_r):
    '''
    q_r: An sqlalchemy row from the queue table with data for one of the 
         three pre-assembly job types.
    args_r: 
        An sqlalchemy row from the args table for this job.

    Returns either a QualityTrimJob, AdapterTrimJob, or ReadNormJob object.
    '''
    jobs = {"quality_trim": QualityTrimJob, "adapter_trim": AdapterTrimJob, 
            "read_normalization": ReadNormJob}
    return jobs[q_r.job_type](q_r.project_id, q_r.job_id, q_r.job_type, 
                args_r.executable, args_r.resources, args_r.pipeline_args, 
                args_r.process_args)


def make_panther_job(q_r,args_r):
    '''
    q_r: 
        An sqlalchemy row from the queue table with data for a Panther job.
    args_r: An sqlalchemy row from the args table for this job.

    Returns a PantherJob object.
    '''
    return PantherJob(q_r.project_id, q_r.job_id, q_r.job_type,
                      args_r.executable, args_r.resources, 
                      args_r.pipeline_args, args_r.process_args)

job_constructors = {"blast": make_hpc_job, "blast_nr": make_hpc_job, 
         "blast_swissprot": make_hpc_job, "blat": make_hpc_job,
         "go":make_blast_association_job, "kegg": make_blast_association_job,
         "pfam": make_hpc_job, "assemble": make_assembly_job, 
         "quality_trim": make_preassembly_job, "panther": make_panther_job,
         "adapter_trim": make_preassembly_job, 
         "read_normalization": make_preassembly_job, "upload": make_upload_job,
         "quantification": make_hpc_job}


def mark_error(jid, session):
    '''
        jid (int): 
            The identifier to mark an error status for
            
        session (netutils.DBSession):
            A session object with an active connection to the Autonomics
            database
            
        Marks the provided job as having a status of 'E' in the jn_mapping table. 
    '''
    jn_mapping = netutils.get_table_object('jn_mapping', session)
    session.conn.execute(jn_mapping.update().where(
              jn_mapping.c.job_id==jid).values(finished='E',
              f_ts=functions.current_timestamp()))


def start_job(job, job_list, pipe_resources, mysql_session, queue):
    ''' 
        job: 
            The object to be started. Must be a subclass of Job.
    
        job_list: 
            A list of jobs that this job should be added to, once started.
    
        pipe_resources: 
            A Resources object that this job will be taking resources
                    from.
    
        mysql_session: 
            A netutils.DBSession object with an active DB connection.
    
        queue: 
            The queue in which this job was scheduled.

        This method checks whether or not the job has any unfinished dependencies.
        If it does not, the job takes its resources from the Resources object and
        is added to the job_list.
    '''
    try:
        dependency = None
        if (queue==settings.normal_queue):  # alternative is que_special which has no dependencies
            #check if this job has unfinished dependencies in the queue
            depends = netutils.get_table_object("jid_dependency", mysql_session)
            jn = netutils.get_table_object("jn_mapping", mysql_session)

            print "checking dependency of job_id: ", job.jid

            s = select([jn, depends], and_(depends.c.job_id==job.jid,
                                           jn.c.job_id==depends.c.depends_on,
                                           jn.c.finished=='N',))
            dependency = s.execute().fetchone()
            if(dependency is None):
              print "dependency is None"
            else:
              print "depends on: ", dependency.depends_on

        if(dependency is None):
            print " calling job.start() for job_id: ", job.jid, " job_type: ", job.job_type
            job.start()
            pipe_resources.give_to(job)
            job_list.append(job)
            mysql_session.conn.execute("UPDATE jn_mapping SET started='Y', \
                   s_ts=CURRENT_TIMESTAMP() WHERE job_id='" + str(job.jid) + "'")
            remove_from_queue(job.jid, queue, mysql_session)
    except Exception as e:
        print "manager in exception after trying to start job"
        sys.stderr.write(e.message + "\n")
        sys.stderr.write("Error starting job: " + job.job_name + "\n")

class Unbuffered:
    '''
    This class is an unbuffered in/out stream which can be used to make 
    print unbuffered.
    '''
    def __init__(self, stream):
        self.stream = stream
    def write(self, data):
        self.stream.write(data)
        self.stream.flush()
    def __getattr__(self, attr):
        return getattr(self.stream, attr)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--name", dest="my_name", default=None)
    parser.add_argument("-u", "--user", dest="user", default=None)
    parser.add_argument("-mu", "--mail-user", dest="mailUser", default=None)
    parser.add_argument("-mh", "--mail-host", dest="mailHost", 
                                                          default="gmail.com")
    parser.add_argument("-du", "--database-user", dest="dbUser", default=None)
    parser.add_argument("-d", "--directory", dest="localDir", 
                      default="/srv/data2/pipeline/", 
                      help="Directory this manager should watch for new jobs.")
    parser.add_argument('--remote-dir', dest='remoteDir', 
                                                default='/scratch/hpc/mcitar/')

    parser.add_argument("-q", "--queue", dest="queue", default="bio", 
                                      help="Queue for submission of HPC jobs.")
    parser.add_argument("--sleep-interval", dest="mainLoopSleepInterval", 
                        default=settings.mainLoopSleepInterval, 
     help="How often this manager should sleep between checking job statuses.")

    manager_queues = [settings.special_queue,settings.normal_queue]
    args = parser.parse_args()
    sys.stdout = Unbuffered(sys.stdout)
    if(args.my_name is None):
        t = datetime.datetime.now()
        my_name = "manager" + str(t.year) + str(t.month) + str(t.day) + \
                                    str(t.hour) + str(t.minute) + str(t.second)
    else:
        my_name = args.my_name

    PIPE_RESOURCES = Resources()
    PIPE_RESOURCES.add_resource("cpu", 61)
    PIPE_RESOURCES.add_resource("a8r", 5)

    HPC_RESOURCES = Resources()
    HPC_RESOURCES.add_resource("cpu", 1024)
    HPC_RESOURCES.add_resource("blast_nr", 2)

    resources_at_location = {Locations.LOCAL: PIPE_RESOURCES, 
                             Locations.HPC: HPC_RESOURCES}

    SLEEP_INTERVAL = float(args.mainLoopSleepInterval)

    if(not args.user is None):
       passwd = raw_input("Enter HPC password:")
       settings.hpc_cred.update(args.user, passwd)
    if(not args.mailUser is None):
        mail_passwd = raw_input("Enter mail account password: ")
        settings.mail_cred.update(args.user, mail_passwd)
    if(not args.dbUser is None):
        db_passwd = raw_input("Enter database passwod:")
        settings.db_cred.update(args.user, db_passwd)
    session = netutils.DBSession("localhost", settings.ZC_DB_NAME, 
                                settings.db_cred.user, settings.db_cred.passwd)
    job_list = []
    finished = []
    lloop_num = 0
    print_res = 0

    while(True):
        lloop_num = lloop_num + 1
        try:
            for queue in manager_queues:
                q = netutils.get_table_object(queue, session)
                results = ""
                if (queue == 'quenew'):
                    results=session.conn.execute("select * from quenew order by priority desc")
                else:
                    s = q.select()
                    results = s.execute()
                for q_row in results.fetchall():
                    if(q_row.job_type == 'upload'):
                        continue
                    args = netutils.get_table_object('args', session)
                    results = args.select(args.c.job_id==q_row.job_id).execute()
                    args_row = results.fetchone()
                    if(args_row is None):
                        default_args = netutils.get_table_object(
                                                          'default_args', session)
                        results = default_args.select(
                                default_args.c.job_type==q_row.job_type).execute()
                        args_row = results.fetchone()
                        if (args_row is None):
                            print "ERROR: can not find args for job_id: ", \
                                       q_row.job_id, " job_type: ", q_row.job_type
                            sys.exit()
                    job = job_constructors[q_row.job_type](q_row, args_row)
                    if(resources_at_location[job.location].has_enough_free(
                                                                  job.resources)):
                        print_res = 1
                        start_job(job, job_list, 
                              resources_at_location[job.location], session, queue)
            finished = []

            for j in job_list:
                state = j.check()
                if(state == JobState.FINISHED):
                    print "manager .. ", j.job_name, " FINISHED"
                    finished.append(j)
                    j.complete()
                    resources_at_location[j.location].take_from(j)
                    print_res = 1
                elif(state == JobState.ERROR):
                    print("manager .. " + j.job_name + " ERROR")
                    finished.append(j)
                    mark_error(j.jid, session)
                    resources_at_location[j.location].take_from(j)
                    print_res = 1
                else:
                    if (lloop_num == 10):
                        print "manager .. ", j.job_name, " RUNNING"

            job_list = [j for j in job_list if not j in finished]

            if print_res == 1:
                print_res = 0
                for loc, resources in resources_at_location.items():
                    print('--- Resources ---')
                    print("free: ", resources.free)
                    print("total: ", resources.totals)
                    print('-----------------')

            gc.collect()

            if (lloop_num == 5):
              print("Sleeping for: " + str(SLEEP_INTERVAL) + " printing this every fifth sleep")
              lloop_num = 0
            SLEEP_INTERVAL = 60
            time.sleep(SLEEP_INTERVAL)
        except exc.OperationalError as e:
            if("MySQL server has gone away" in e.message):
                print("MySQL server has gone away, restarting session.\n")
                session = netutils.DBSession("localhost", settings.ZC_DB_NAME,
                                settings.db_cred.user, settings.db_cred.passwd)

if __name__ == '__main__':
    main()
