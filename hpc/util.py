#-------------------------------------------------------------------------------
# Name:        util
# Purpose:     Utilities for submitting, checking, and retrieving hpc jobs
#
# Author:      Mathew Citarella
#
# Created:     18/09/2012
# Copyright:   (c) mat 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------


class JobServer:

    jobs = {}
    finished = []

    def __init__(self, serverName, user, passwd, mailUser, mailPasswd, host, mailHost):
        self.name = serverName
        self.u = user
        self.p = passwd
        self.mailU = mailUser
        self.mailP = mailPasswd
        self.h = host
        self.mailH = mailHost

    def submit(self, j):
        pass

    def _startMonitor():



class Job:

    def __init__(self, name, t, ex, params = {}, mem = "500mb", procs = 300, q = "bio", wall="192:00:00", email="morozhpc@gmail.com"):
        self.name = name
        self.type = t
        self.program = ex
        self.params = params
        self.memory = mem
        self.procs = procs
        self.queue = q
        self.walltime = wall
        self.email = email
        self._initModule()

    def _initModule(self):
        if(self.type == "blast"):
            self.module = "module load ncbi_blast"




def main():
    pass

if __name__ == '__main__':
    main()
