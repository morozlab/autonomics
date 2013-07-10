import argparse
import sys
from sqlalchemy import func
from sqlalchemy.sql import and_, select, or_
from autonomics import netutils

def project_duration(pid, session, count_nr=False):
    jn_mapping = netutils.get_table_object('jn_mapping',
                                         session)

    if(not count_nr):
        s = select([
        func.unix_timestamp(func.min(jn_mapping.c.s_ts)),
        func.unix_timestamp(func.max(jn_mapping.c.f_ts))
               ], and_(jn_mapping.c.project_id==pid,
                       func.unix_timestamp(
                           jn_mapping.c.s_ts)!=0,
                       jn_mapping.c.job_type!='blast_nr')
        )

    else:

    res = session.conn.execute(s)
    row = res.fetchone()
    if(row is None):
        sys.stderr.write("Warning - no jobs set for \
                          project: " +  str(pid) + "\n"
        )

        return 0

    elif(row[1] is None or row[1] is None):
        sys.stderr.write("Warning - received None for \
                          either max finish time or min \
                          start time.\n"
        )
        return 0

    return row[1] - row[0]


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--project-ids', dest='pids',
                        default=None, nargs='+',
                        help='The project_ids of the \
                        projects you wish to work on.'
    )

    parser.add_argument('--pid-range', dest='pid_range',
                        default=None, help='A range of \
                        project ids specified in the \
                        form: start:end'
    )

    parser.add_argument('--durations', dest='durations',
                        default=False,
                        const=True, action='store_const',
                        help='Calculate run time durations\
                        for an entire project.'
    )

    args = parser.parse_args()

    s = netutils.make_db_session()

    if(args.durations):
        pids = args.pids
        if(not args.pid_range is None):
            split = args.pid_range.split(":")
            start = int(split[0])
            end = int(split[1])
            pids = [i for i in xrange(start, end + 1)]


        for pid in pids:
            print("\t".join(
                       [str(pid),
                        str(project_duration(int(pid), s))
                       ]
                       )
            )



if __name__ == '__main__':
    main()
