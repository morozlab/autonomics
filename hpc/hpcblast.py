from datetime import datetime
from files.sequences import Utils
from hpc import qsub
from network import access
from network import ssh
import argparse
import os.path
import shutil
import sys


parser = argparse.ArgumentParser(description="This script splits a BLAST query and submits the required jobs to the HPC cluster for computation.")
parser.add_argument("-q", "--query-file", dest="queryFile", required=True)
parser.add_argument("-d", "--database", dest="db", required=True)
parser.add_argument("--host", dest='host', default="submit.hpc.ufl.edu")
parser.add_argument("--username", dest="user", required=True)
parser.add_argument("--passwd", dest="passwd", required=True)
parser.add_argument("--mail-server", dest="emailHost", default="imap.gmail.com")
parser.add_argument("--mail-user", dest="mailUser", required=True)
parser.add_argument("--mail-password", dest="mailPass", required=True)
parser.add_argument("--mailbox", dest='mailBox', default='INBOX')
parser.add_argument("--seqs-per-file", dest="seqsPerFile", default=10000)
parser.add_argument("-n", "--job-name", dest="jobName", default="parallelBlast")
parser.add_argument("--queue", dest='targetQ', default='bio')
parser.add_argument("--directory", dest="baseDir", default="/scratch/hpc/mcitar/")
parser.add_argument("--memory-request", dest="memReq", default="500mb")
parser.add_argument("--qsub-only", dest="qsubOnly", default=False, const=True, action='store_const')
parser.add_argument("-o", "--output-file",  dest="output", default="blast.out")
parser.add_argument("-m", "--output-format", dest="outFormat", default=0)
parser.add_argument("-b", "--num-align", dest="numAlign", default=1)
parser.add_argument("-v", "--num-oneline", dest="numOneLine", default=1)
parser.add_argument("-e", "--evalue", dest="evalue", default=10)
parser.add_argument("-p", "--program", dest="prog", default="blastn")



def main():
    #parse in the input parameters
    args = parser.parse_args()

    #set up the target BLAST database
    if(args.db == "nr" or args.db == "swissprot"): dbParams = "-db $BLASTDB/" + args.db
    else: dbParams = "-db " + args.db

    #set up result limiting
    if("6" in args.outFormat): numAlnParams = " -max_target_seqs " + str(args.numAlign)
    else: numAlnParams = " -num_alignments " + str(args.numAlign) + " -num_descriptions " + str(args.numOneLine)

    #create credential objects
    sshCredentials = access.Credentials(args.host, args.user, args.passwd)
    mailCredentials = access.Credentials(args.emailHost, args.mailUser, args.mailPass)

    #get the current time
    thisTime = datetime.now()

    #make a temporary local directory for this project
    tmpDir = "hpcblast_" + args.jobName + str(thisTime.year) + str(thisTime.month) + str(thisTime.day) + "_" + str(thisTime.hour) + str(thisTime.minute) + str(thisTime.second)
    os.mkdir(tmpDir)
    #move the query file to the temp directory
    shutil.copy(args.queryFile, tmpDir)
    #change to the temporary directory
    os.chdir(tmpDir)

    #split the query sequences
    queryDir, queryFile = os.path.split(args.queryFile)
    files = Utils.splitSeqFile(queryFile, args.seqsPerFile)

    #create list of jobs
    jobs = []

    #create the remote working directory for this blast
    c = ssh.Connection(args.host, username=args.user, password=args.passwd)
    runDir = args.baseDir + "/" + args.jobName + str(thisTime.year) + str(thisTime.month) + str(thisTime.day) + "_" + str(thisTime.hour) + str(thisTime.minute) + str(thisTime.second)
    command = "mkdir " + runDir
    c.execute(command)


    jobNo = 1
    #create a qsub file for each query
    for f in files:
        #create a new qsub job for this query fragment (job name, working directory, input file, output file)
        q = qsub.job(args.jobName + "_" + str(jobNo), runDir, f, f + ".out")
        #set the header for this qsub job, containing the PBS directives
        q.createPBSHeader("192:00:00", "morozhpc@gmail.com", args.targetQ, args.memReq, 1, 1)
        #add the directive for blast to the qsub script
        q.appendQsubCommand("module load ncbi_blast")
        q.appendQsubCommand(args.prog +  " -query " + runDir + "/" + f + " "  + dbParams + " -evalue " + str(args.evalue) + " -outfmt " + str(args.outFormat) + numAlnParams + " -out " + runDir + "/" + f + ".out")
        #check if we're just writing the qsub file
        if(args.qsubOnly):
            q.writeQsub()
            sys.exit()

        #submit the job to the hpc cluster
        q.submit(c)
        #add to the running jobs
        jobs.append(q)
        jobNo += 1

    #close the ssh connection
    c.close()

    #monitor the running jobs until all have finished execution
    monitor = qsub.monitor(jobs, sshCredentials, mailCredentials, args.mailBox, True)
    monitor.start(300)

    #combine the output into a single file
    combined = open("../" + args.output, 'w')
    c = ssh.Connection(args.host, username=args.user, password=args.passwd)
    for q in jobs:
        #retrieve the file from the cluster and store the filename
        singleOut = q.retrieveOutput(c)
        #open the output file for reading
        tmp = open(singleOut, 'r')
        #write the contents of the file to the combined file, close, and remove the file
        combined.write(tmp.read())
        tmp.close()
        os.remove(singleOut)
    c.close()
    combined.close()

    #clean up temporary files
    c = ssh.Connection(args.host, username=args.user, password=args.passwd)
    c.execute("rm -r " + runDir + "/")
    c.close()
    os.chdir("../")
    shutil.rmtree(tmpDir)




if __name__ == '__main__':
    main()
