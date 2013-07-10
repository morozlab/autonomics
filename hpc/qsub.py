#-------------------------------------------------------------------------------
# Name:        Qsub
# Purpose:     This module contains various functions for creating Qsub scripts.
#
# Author:      mat
#
# Created:     20/07/2012
# Copyright:   (c) mat 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from network import ssh
import imaplib
import os
import ssl
import sys
import time

class job:

    def __init__(self, n, d, inputFile, outputFile = None):
        self.qsub = ""
        self.commands = ""
        self.name = n
        self.memReq = 0
        self.walltime = ""
        self.email = ""
        self.q = ""
        self.nodes = 0
        self.ppn = 0
        self.jobID = ""
        self.localDir, self.localInput = os.path.split(inputFile)
        if(self.localDir == ""):
            self.localDir = "./"
        else:
            self.localDir = self.localDir + "/"
        self.remoteInput = self.localInput
        self.remoteDir = d
        self.status = "idle"
        self.finishedCount = 0;
        if(outputFile == None):
            self.outputFile = self.name + ".out"
        else:
            self.outputFile = outputFile


    def _getConnection(self, h, u, p):
        return ssh.Connection(h, username=u, password=p)

    def appendQsubCommand(self, command):
        self.commands += command + "\n"

    def appendToPBS(self, directive):
        self.qsub += directive + "\n"

    def createPBSHeader(self, walltime, email, q, mem, nodes, ppn):
        self.qsub = "#! /bin/bash\n";
        self.qsub += "#PBS -r n\n";
        self.qsub += "#PBS -N " + self.name + "\n";
        self.qsub += "#PBS -o " + self.name + ".out\n";
        self.qsub += "#PBS -e " + self.name + ".err\n";
        if(q == 'billed'): self.qsub += "#PBS -W group_list=billed\n"
        else: self.qsub += "#PBS -q " + q + "\n";
        self.qsub += "#PBS -m a\n";
        self.qsub += "#PBS -M " + email + "\n";
        self.qsub += "#PBS -l walltime=" + walltime + "\n";
        self.qsub += "#PBS -l pmem=" + mem + "\n";
        self.qsub += "#PBS -l nodes=" + str(nodes) + ":ppn=" + str(ppn) + "\n";
        self.qsub += "\n"
        self.memReq = int(mem[:-2])
        self.walltime = walltime
        self.email = email
        self.q = q
        self.nodes = nodes
        self.ppn = ppn

    def retrieveOutput(self, c):
        sys.stdout.write("Retrieving: " + self.remoteDir + "/" + self.outputFile + "\n")
        try:
            c.get(self.remoteDir + "/" + self.outputFile, self.localDir + self.outputFile)
        except IOError:
            sys.stdout.write("Error " + str(IOError.errno) + " retrieving: " + self.remoteDir + "/" + self.outputFile  + "\n")
        return self.outputFile

    def resubmit(self, c, memIncrease = 0):
        self.createPBSHeader(self.walltime, self.email, self.q, str(self.memReq + memIncrease) + "mb", self.nodes, self.ppn)
        self.submit(c)

    def setDir(self, directory):
        self.remoteDir = directory

    def setMemRequest(self, req):
        self.memReq = req

    def submit(self, c):
        #write the qsub submit file
        self.writeQsub()
        #put this qsub script on the server
        c.put(self.localDir + self.name + ".qsub", self.remoteDir + "/" + self.name + ".qsub")
        #put input files on the server
        c.put(self.localDir + self.localInput, self.remoteDir + "/" + self.remoteInput)
        #remove the temporary qsub file
        #os.remove(self.name + ".qsub")
        #remove the temporary input file
        #os.remove(self.inputFile)
        #start the job
        command = "cd " + self.remoteDir + ";"
        command += "qsub " + self.name + ".qsub"
        jobID = c.execute(command)
        self.jobID = jobID[0].strip()
        self.status = "running"
        self.finishedCount = 0

    def writeQsub(self):
        out = open(self.localDir + self.name + ".qsub", 'w')
        out.write(self.qsub)
        out.write(self.commands)
        out.close()

class monitor:

    def __init__(self, jobs, sshCred, mailCred, mailbox = "INBOX", restart = False,):
        self.jobs = jobs
        self.restart = restart
        self.mailHost = mailCred.host
        self.mailUser = mailCred.user
        self.mailPass = mailCred.passwd
        self.mailbox = mailbox
        self.host = sshCred.host
        self.user = sshCred.user
        self.passwd = sshCred.passwd
        self.mailServer = None
        self.interval = 0
        self.retryInterval = 10
        self.retries = 5
        self.finished = False

    def _mailConnect(self):
        retries = self.retries
        while(True):
            try:
                mailConnect = imaplib.IMAP4_SSL(self.mailHost, '993')
                mailConnect.login(self.mailUser, self.mailPass)
                return mailConnect
            except imaplib.IMAP4_SSL.abort:
                if(retries > 0):
                    sys.stderr.write("IMAP abort exception")
                    retries -= 1
                    time.sleep(self.retryInterval)
                else:
                    raise
            except ssl.SSLError:
                if(retries > 0):
                    sys.stderr.write("SSLError exception!")
                    retries -= 1
                    time.sleep(self.retryInterval)
                else:
                    raise

    def _sshConnect(self):
        try:
            return ssh.Connection(self.host, username=self.user, password=self.passwd)
        except AttributeError:
            sys.stderr.write("SSH server disappeared, sleeping for two minutes and trying again.\n")
            time.sleep(120)


    def _exceededMemory(self, job):
        retries = self.retries
        while(True):
            try:
                typ, msg_ids = self.mailServer.search(None, '(SUBJECT "PBS JOB ' + str(job.jobID) + '" BODY "exceeded MEM usage hard limit")')
                if(typ == "OK"):
                    if(msg_ids[0] != ''):
                        typ, data = self.mailServer.fetch(msg_ids[0], '(BODY.PEEK[TEXT])')
                        return int(self._parseEmailBody(data, "mem"))
                return 0

            except(ssl.SSLError, imaplib.IMAP4_SSL.abort):
                if(retries > 0):
                    sys.stderr.write("Exception raised during mail access!")
                    retries -= 1
                    #sleep, then start a new connection to mail server
                    time.sleep(self.retryInterval)
                    self.mailServer = self._mailConnect()
                    self.mailServer.select(self.mailbox)
                else:
                    raise


    def _exitWithError(self, job):
        return False

    def _parseEmailBody(self, response, directive):
        if(directive == "mem"):
            for response_part in response:
                if(isinstance(response_part, tuple)):
                    text = response_part[1]
                    it = iter(text.splitlines())
                    #iterate over the email text
                    for line in it:
                        #if this is the memory usage line, parse memory usage
                        if("exceeded MEM usage hard limit" in line):
                            return line.split("(")[1].split(">")[0].strip()

        return ""


    def checkJobStatus(self, job):

        #check if this job terminated due to exeeding memory requests
        exceeded = self._exceededMemory(job)
        if(exceeded > 0):
            job.status = "exceedMem"
            if(self.restart):
                #set the new memory request to the amount used when the job was killed
                job.memReq = exceeded
            return job.status

        #check if this job exited with an error status ######IMPLEMENT THIS##############
        elif(self._exitWithError(job)):
            job.status = "error"
            return job.status

        command =  "qstat -u mcitar | grep -c " + job.jobID[:16]
        try:
            running = int(self.serverConnect.execute(command)[0])
        except ValueError:
            #if we get to this section of code, grep has failed and we have no way of knowing anything about the process, abandon it
            sys.stderr.write("Error retrieving status of running job " + job.name + " - abandoning.\n")
            running = -1

        if(running == 0):
            job.status = "finished"
        elif(running == -1):
            job.status = "error"
        else:
            job.status = "running"

        return job.status

    def checkJobsFinished(self):
        finished = 0
        self.mailServer = self._mailConnect()
        self.serverConnect = self._sshConnect()
        typ, data = self.mailServer.select(self.mailbox)
        index = 0
        errors = []
        for job in self.jobs:
            #update the status of this job
            s = self.checkJobStatus(job)
            if(s == 'finished'):
                job.finishedCount += 1
                #give error messages time to reach the mail server
                if(job.finishedCount * self.interval >= 300):
                    finished += 1
                    job.status ="processed"

            elif(s == 'exceedMem'):
                #check if we should restart this job
                if(self.restart):
                    print("Job " + job.name + " exceeded its memory allowance, restarting with " + str(job.memReq + 200))
                    job.resubmit(self.serverConnect, 200)
                else:
                    finished += 1
                    job.status = "error"
            elif(s == "error"):
                errors.append(index)
            index += 1

        for i in reversed(errors):
            self.jobs.pop(i)
        self.mailServer.close()
        self.mailServer.logout()
        self.serverConnect.close()

        if(finished == len(self.jobs)):
            self.finished = True

    def restartJob(self, job):
        pass

    def start(self, t = 60, blocking = True):
        self.interval = t
        #dont start a new thread
        if(blocking):
            while(not self.finished):
                self.checkJobsFinished()
                time.sleep(t)
