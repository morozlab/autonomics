import argparse
import os
import subprocess
from sqlalchemy.sql import and_
from autonomics import settings, netutils
from autonomics.file_io import Record
from autonomics.utility import row2dict
from autonomics.queue import remove_from_queue, configured_to_run


'''
UPCOMING FEATURES

1) Rename project option.
2) Workflow and Location Management
3) Update configuration to allow specifying which jobs a given job depends on

'''


def add_adapters(adapters, session):
    ka = netutils.get_table_object("known_adapters", session)
    if(os.path.isfile(adapters)):
        fh = open(adapters, 'r')
        for line in fh:
            line = line.rstrip("\n")
            el = line.split("\t")
            if(len(el) > 1):
                i = ka.insert().values(adapter_name=el[0], adapter_sequence=el[1])
                i.execute()
    else:
        el = adapters.split(":")
        if(len(el) >  1):
            i = ka.insert().values(adapter_name=el[0], adapter_sequence=el[1])
            i.execute()


def add_job(project_name, job_type, session):
    pn = netutils.get_table_object('pn_mapping', session)
    jn = netutils.get_table_object('jn_mapping', session)
    pid = netutils.get_pid(project_name, session)
    if(pid is None):
        print("Could not add job: " + job_Type + ", no project_id found for: " + project_name)
        return False
    job_name = project_name + "_" + job_type
    session.conn.execute(jn.insert().values(project_id=pid, job_type=job_type, job_name=job_name))

    '''
    try:
        if(job_name is None):
            job_name = project_name + "_" + job_type
        session.conn.execute(jn.insert().values(project_id=pid, job_type=job_type, job_name=job_name))
    except Exception as e:
        print(e.message)
        return False
    '''
    return True


def add_project(project_name, session):
    pn = netutils.get_table_object('pn_mapping', session)
    i = pn.insert().values(project_name=project_name)
    session.conn.execute(i)
    return True


def assign_adapters(adapters, project_names, session):
    #get a list of adapters
    pa = netutils.get_table_object("project_adapters", session)
    ka = netutils.get_table_object("known_adapters", session)
    pn_mapping = netutils.get_table_object("pn_mapping", session)

    #first, get a list of the adapter names
    adapter_names = [adapters]
    if(os.path.isfile(adapters)):
        adapter_names = listify_file(adapters)

    #get the adapter ids
    adapt_ids = []
    for adapt in adapter_names:
        res = ka.select(ka.c.adapter_name==adapt).execute()
        row = res.fetchone()
        if(not row is None):
            adapt_ids.append(row.adapter_id)

    projects = [project_names]
    if(os.path.isfile(project_names)):
        projects = listify_file(project_names)

    #iterate over the projects and assign adapters
    for project in projects:
        res = pn_mapping.select(pn_mapping.c.project_name==project).execute()
        row = res.fetchone()
        if(row is None):
            print("Could not look up pid for project: " + project + ", skipping.")
            continue

        pid = row.project_id
        for adapt_id in adapt_ids:
            i = pa.insert().values(project_id=pid, adapter_id=adapt_id)
            i.execute()


def do_to_projects(project_names, func, session=None):
    if(session is None):
        session = netutils.make_db_session()

    for pn in project_names:
        pid = netutils.get_pid(pn, session)
        if(pid is None):
            print("Could not find project_id for: " + pn + ", skipping.")
            continue

        func(pid, session)


def get_project_names(project_args):
    ret = [project_args]
    if(os.path.exists(project_args) and os.path.isfile(project_args)):
        ret = []
        fh = open(project_args, 'r')
        for line in fh:
            line = line.rstrip('\n')
            ret.append(line.strip())
        fh.close()

    return ret


def listify_file(file_name):
    ret = []
    fh = open(file_name, 'r')
    for line in fh:
        line = line.rstrip('\n')
        ret.append(line.strip())
    fh.close()

    return ret


def set_for_reschedule(pid, q_name, session):
    initp = netutils.get_table_object('init_projects', session)
    submitp = netutils.get_table_object('submitted_projects', session)
    jn_mapping = netutils.get_table_object('jn_mapping', session)
    queue = netutils.get_table_object(q_name, session)

    trans = session.conn.begin()

    try:
        session.conn.execute(initp.delete(initp.c.project_id==pid))
        session.conn.execute(submitp.delete(submitp.c.project_id==pid))
        session.conn.execute(queue.delete(queue.c.project_id==pid))
        u = jn_mapping.update().where(jn_mapping.c.project_id==pid).values(started='N', finished='N')
        trans.commit()
        print("Successfully removed project " + str(pid) + " from queue.")

    except Exception as e:
        print("Error removing project " + str(pid) + " from queue.")
        trans.rollback()


def remove_project(pid, session, remove_files=False):
    #first, open the appropriate table objects
    args = netutils.get_table_object('args', session)
    config = netutils.get_table_object('configuration', session)
    non_mirrored=True
    #these projects aren't mirrored currently, and may not exist
    try:
        initp = netutils.get_table_object('init_projects', session)
        submitp = netutils.get_table_object('submitted_projects', session)
        padapts = netutils.get_table_object('project_adapters', session)
        queue = netutils.get_table_object('quenew', session)
    except Exception as e:
        print("Non-mirrored tables not found, not making objects for them.")
        non_mirrored=False

    jn_mapping = netutils.get_table_object('jn_mapping', session)
    pn_mapping = netutils.get_table_object('pn_mapping', session)
    runname_to_pid = netutils.get_table_object('runname_to_pid', session)

    trans = session.conn.begin()
    try:
        #get the project name
        res = session.conn.execute(pn_mapping.select(pn_mapping.c.project_id==pid))
        row = res.fetchone()
        pn = row.project_name

        #delete stuff for jobs first
        res = session.conn.execute(jn_mapping.select(jn_mapping.c.project_id==pid))
        for row in res.fetchall():
            jid = row.job_id
            #delete args
            session.conn.execute(args.delete(args.c.job_id==jid))
            #delete from jn_mapping
            session.conn.execute(jn_mapping.delete(jn_mapping.c.job_id==jid))
            #delete from queue
            if(non_mirrored):
                session.conn.execute(queue.delete(and_(queue.c.project_id==pid, queue.c.job_id==jid)))

        #delete stuff at the project level
        if(non_mirrored):
            session.conn.execute(initp.delete(initp.c.project_id==pid))
            session.conn.execute(submitp.delete(submitp.c.project_id==pid))
            session.conn.execute(padapts.delete(padapts.c.project_id==pid))

        session.conn.execute(pn_mapping.delete(pn_mapping.c.project_id==pid))
        session.conn.execute(runname_to_pid.delete(runname_to_pid.c.project_id==pid))
        if(remove_files):
            project_dir = settings.home_dir + pn + "/"
            while(True):
                answer = raw_input("Are you sure you want to delete project files located in " + project_dir + "? (y/n)")
                if(answer == 'y'):
                    shutil.rmtree(settings.home_dir + pn + "/")
                    break
                elif(answer == 'n'):
                    break

        trans.commit()

    except Exception as e:
        trans.rollback()
        print("Error during delete of project: " + str(pid))
        print(e.message)


def set_args(project_name, job_type, arg_list, session):
    pid = netutils.get_pid(project_name, session)
    jid = netutils.get_jid(pid, job_type, session)
    if(jid is None):
        print("No job entry found for: " + project_name + " " + job_type)
        return

    def_args = netutils.get_table_object('default_args', session)
    args = netutils.get_table_object('args', session)

    def_arg_res = session.conn.execute(def_args.select(def_args.c.job_type==job_type))
    def_dict = row2dict(def_arg_res.fetchone())

    passed_args = {}
    for arg in arg_list:
        el = arg.split(";")
        passed_args[el[0]] = el[1]

    def_dict.update(passed_args)
    args_for_insert = def_dict
    args_for_insert_jid = def_dict.copy()
    args_for_insert_jid['job_id'] = jid

    i = args.insert().values(args_for_insert_jid)
    u = args.update().where(args.c.job_id==jid).values(args_for_insert)
    netutils.update_on_exists(i, u)


def set_configuration(project_name, job_type, code, session):
    pid = netutils.get_pid(project_name, session)
    if(pid is None):
        print("Couldn't find a project_id for:" + project_name + ", can't set job configurations")
        return False

    conf = netutils.get_table_object('configuration', session)
    i = conf.insert().values(project_id=pid, job_type=job_type, code=code)
    u = conf.update().where(and_(conf.c.project_id==pid,
                                 conf.c.job_type==job_type)).values(code=code)
    netutils.update_on_exists(i, u)


def set_workflow(project_name, session):
    deps = netutils.get_table_object('dependency', session)
    job_deps = netutils.get_table_object('jid_dependency', session)
    pid = netutils.get_pid(project_name, session)

    res = session.conn.execute(deps.select())
    for row in res.fetchall():
        jtype = row.job_type
        dep = row.depends_on
        jid = netutils.get_jid(pid, jtype, session)
        depid = netutils.get_jid(pid, dep, session)

        if(not jid is None and
           not depid is None and configured_to_run(pid, jtype, session)):
            i = job_deps.insert().values(job_id=jid,
                                     depends_on=depid)
            u = job_deps.update().where(job_deps.c.job_id==jid).values(depends_on=depid)
            netutils.update_on_exists(i, u)


def update_pipe_args(args, project_names, job_type, session):
    args_table = netutils.get_table_object("args", session)
    queue = netutils.get_table_object("queue", session)

    for pn in project_names:
        pid = netutils.get_pid(pn, session)
        jid = netutils.get_jid(pid, job_type, session)
        if(jid is None or pid is None):
            print("Error finding jobtype: " + job_type + " for project: " + pn + ", skipping.")
            continue

        u = args_table.update().where(args_table.c.job_id==jid).values(pipeline_args=args)
        u.execute()
        u = queue.update().where(and_(queue.c.project_id==pid, queue.c.job_id==jid)).values(pipeline_args=args)
        u.execute()


def restart_projects(project_args, session):
    pn_mapping = netutils.get_table_object("pn_mapping", session)
    queue = netutils.get_table_object("queue", session)
    completed_projects = netutils.get_table_object("completed_projects", session)
    project_names = get_project_names(project_args)

    for pn in project_names:
        #get the pid for this project
        res = pn_mapping.select(pn_mapping.c.project_name==pn).execute()
        row = res.fetchone()
        if(row is None):
            print("Project name: " + pn + " not found in system, skipping.")
            continue

        pid = row.project_id

        #update completed jobs
        d = completed_projects.delete().where(completed_projects.c.project_id==pid)
        d.execute()

        #update the queue
        u = queue.update().where(queue.c.project_id==pid).values(started='N', finished='N')
        u.execute()


def main():


    session = netutils.make_db_session()
    parser = argparse.ArgumentParser()
    parser.add_argument("--add-adapters", dest='add_adapters', default=None, help='Add known adapters to the system. Value can either be a single adapter to add in the form name:sequence or the path to a list of adapters. If a file is specified, must be a tab-delimited file where the first column in the adapter name and the second is the adapter sequence.')
    parser.add_argument('--add-jobs', dest='add_jobs', default=None, nargs='+', help='Space-separated list of job_types you wish to add to the project(s) given by --project-names.')
    parser.add_argument('--add-project', dest='add_project', default=False, const=True, action='store_const', help='Add the projects given by --project-names to the system, if they don\'t already exist.')
    parser.add_argument("--assign-adapters", dest='assign_adapters', default=None, help='Assign adapters to the supplied projects. Value can either be a single adapter name, or the path to a file of adapter names.')
    parser.add_argument("--assign-workflow", dest='assign_workflow', default=False, const=True, action='store_const', help='Assign the specified projects a given workflow.')
    parser.add_argument('--mark-for-dispatch', dest='mark_for_dispatch', default=False, const=True, action='store_const', help='Mark the project(s) specified in --project-names as ready for dispatcher.py to add to the queue.')
    parser.add_argument("--full-project-remove", dest="remove_project", default=False, const=True, action='store_const', help="Delete the project from system tables.")
    parser.add_argument("--restart-projects", dest="restart_projects", default=False, const=True, action='store_const', help="Restarts projects in the queue.")
    parser.add_argument('--remove-from-queue', dest='remove_from_queue', default=False, const=True, action='store_const', help='Removes a project from the specified queue.')
    parser.add_argument('--set-args', dest='set_args', default=None, nargs='+', help='Space-separated list of string arguments consisting of a job_type and arguments for the job. Arguments are set for each job of job_type for the projects specified by --project-names. Format for each string is: job_type|arg_name1;arg_val1|arg_name2;arg_val2.')
    parser.add_argument('--set-config', dest='set_config', default=None, nargs='+', help='Set configuration for the jobs supplied with this flag. Format should be --set-config job_type1:code1 job_type2:code2 etc.')
    parser.add_argument("--update-pipe-args", dest="update_pipe_args", default=None, help="Update pipeline_args for projects in the args table. Value is the new value to be set for queue.pipeline_args")
    parser.add_argument("--job-type", dest="job_type", help="Which job type do you want this operation to affect?")
    parser.add_argument("--project-names", dest="project_names", default=None,  help="Either the name of a project or the path to a list of project names, one per line. Operations will be carried out on these projects.")

    args = parser.parse_args()

    pns = None
    if(not args.project_names is None):
        pns = get_project_names(args.project_names)
    else:
        print("Warning, no project names supplied on command line.")

    if(not args.add_adapters is None):
        add_adapters(args.add_adapters, session)

    if(args.add_project):
        for name in pns:
            add_project(name, session)

    if(not args.add_jobs is None):
        for name in pns:
            for job in args.add_jobs:
                add_job(name, job, session)

    if(not args.set_config is None):
        for name in pns:
            for job_code in args.set_config:
                els = job_code.split(":")
                set_configuration(name, els[0], els[1], session)

    if(args.assign_workflow):
        for name in pns:
            set_workflow(name, session)

    if(not args.set_args is None):
        for name in pns:
            for arg_str in args.set_args:
                els = arg_str.split("|")
                if(len(els) > 0):
                    job_type = els[0]
                    arg_lst = els[1:]
                    set_args(name, job_type, arg_lst, session)

    if(args.mark_for_dispatch):
        for name in pns:
            #make a SRC_UPLOADED file in each project directory specified
            p = subprocess.Popen("touch " + settings.home_dir + name + "/SRC_UPLOADED", shell=True)
            p.wait()

    if(not args.assign_adapters is None):
        assign_adapters(args.assign_adapters, args.project_names, session)

    if(args.remove_from_queue):
        if(not args.project_names is None):
            project_names = get_project_names(args.project_names)
            pn_mapping = netutils.get_table_object('pn_mapping', session)

            for name in project_names:
                res = session.conn.execute(pn_mapping.select(pn_mapping.c.project_name==name))
                row = res.fetchone()
                if(row is None):
                    print("Could not find a project_id for project name: " + name + ", skipping")
                    continue

                remove_from_queue(row.project_id, session)

    if(args.remove_project):
        if(not args.project_names is None):
            project_names = get_project_names(args.project_names)
            do_to_projects(project_names, remove_project, session)


    if(args.restart_projects):
        restart_projects(args.project_names, session)

    if(not args.update_pipe_args is None):
        update_pipe_args(args.update_pipe_args, get_project_names(args.project_names), args.job_type, session)


if __name__ == "__main__":
    main()
