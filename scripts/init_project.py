from zeroclick import settings
import argparse

def make_job_insert(vals):

    return {"job_type": vals[0], "executable": vals[1], "pipeline_args": vals[2], "process_args": vals[3], "loc": vals[4], "resources": vals[5]}

def make_default_insert():

    return {"blast_nr": make_blast_nr(), "blast_swissprot": make_blast_swissprot(), "pfam": make_pfam(), "go": make_go(), "kegg": make_kegg()}


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('-p', '--project-name', dest='project_name', required=True)
    parser.add_argument('--pipeline-args', dest='pipe_args', nargs='+', default=None)
    parser.add_argument('--process-args', dest='proc_args', nargs='+', default=None)
    parser.add_argument('--default', dest='default', const=True, default=False, action='store_const')
    parser.add_argument('--translate', dest='translate', const=True, default=False, action='store_const')


    args = parser.parse_args()

    if(args.translate):
        translate_str = "--translate"
    else:
        translate_str = ""

    #check if individual jobs are specified
    if(args.default):
        #add a row in the queue for every existing job
    else:



if __name__ == '__main__':
    main()
