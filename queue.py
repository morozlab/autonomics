from sqlalchemy.sql import and_, select, insert, delete
from zeroclick import netutils, settings


def configured_to_run(pid, jtype, session):
    config = netutils.getTableObject('configuration', session)
    s = config.select(and_(config.c.project_id==pid,
                           config.c.job_type==jtype))
    res = session.conn.execute(s)
    row = res.fetchone()
    if(not row is None):
        if(row.code != '-'):
            return True

    return False


def dependencies_satisfied(jid, session):
    jn = netutils.getTableObject('jn_mapping', session)
    depends = netutils.getTableObject('jid_dependency', session)

    s = select([jn, depends], and_(depends.c.job_id==jid,
                                   jn.c.job_id==depends.c.depends_on,
                                   jn.c.finished!='Y'))

    dependency = s.execute().fetchone()
    if(dependency is None):
        return True
    else:
        return False


def elligible_jobs(session):
    query = ''' SELECT project_id, jn_mapping.job_id, job_type
                FROM jn_mapping
                WHERE queued='N' and started='N'
            '''

    res = session.conn.execute(query)
    for row in res.fetchall():
        if(dependencies_satisfied(row.job_id, session) and
           configured_to_run(row.project_id, row.job_type, session)):
            yield row


def insert_in_queue(pid, jid, job_type, queue, session):
    q_table = netutils.getTableObject(queue, session)
    jn_mapping = netutils.getTableObject('jn_mapping', session)
    i = q_table.insert().values(project_id = pid,
                                job_id=jid,
                                job_type=job_type)
    u = jn_mapping.update().where(jn_mapping.c.job_id==jid).values(queued='Y')
    t = session.conn.begin()
    try:
        session.conn.execute(i)
        session.conn.execute(u)
        t.commit()
    except Exception as e:
        t.rollback()
        l = ['Warning, unable to insert job: ',
             str(jid), ' into queue: ' , queue,
             ', already exists.']
        print(''.join(l))
        print(e.message)


def queue_tracking(jid, queue, session):
    queued_jobs = netutils.getTableObject('queued_jobs', session)

    i = queued_jobs.insert().values(job_id=jid, queue=queue)
    session.conn.execute(i)


def queue_all(session):
    #basic logic for deciding where to queue
    #can make more complicated as we improve queuing functionality

    #ADD: Method to ensure that all dependencies satisfied

    for job in elligible_jobs(session):
        q = settings.normal_queue
        if(job.project_id == 0):
            q = settings.special_queue

        insert_in_queue(job.project_id, job.job_id,
                        job.job_type, q, session)


def queue_job(jid, queue, session):
    jn_mapping = netutils.getTableObject('jn_mapping', session)
    res = session.conn.execute(jn_mapping.select(jn_mapping.c.job_id==jid))
    job_row = res.fetchone()
    insert_in_queue(job_row.project_id, jid,
                    job_row.job_type, queue, session)


def remove_from_queue(jid, q_name, session):
    queue = netutils.getTableObject(q_name, session)
    d = queue.delete(queue.c.job_id==jid)
    session.conn.execute(d)


def remove_queue_tracking(jid, session):
    queued_jobs = netutils.getTableObject('queued_jobs', session)
    d = queued_jobs.delete(queued_jobs.c.job_id==jid)
    session.conn.execute(d)
