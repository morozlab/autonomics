import os, xmlrpclib
from zeroclick import netutils, utility, settings

debug = False

class Location:

    def __init__(self, host, sandbox,
                 cred_file=None, xfer_mthd='sftp'):
        self.host = host
        self.xfer_mthd = xfer_mthd
        self.sandbox = sandbox
        self.sandbox = os.path.normpath(self.sandbox)
        if(not self.sandbox.endswith("/")):
            self.sandbox = self.sandbox + "/"
            
        if(cred_file is None):
            try:
                self.credentials = settings.creds[host]
            except Exception as e:
                if(debug):
                    print('''Could not get credentials for host:
                          ''' + str(host))
                self.credentials = None
                    
        else:
            self.credentials = settings.Credentials(from_file=cred_file)

        self.task_count = 0
            
    def _res_connect(self):
        return xmlrpclib.ServerProxy('http://localhost:' + str(settings.RESOURCE_SERVER_PORT))

    def available_resources(self):
        s = self._res_connect()
        return s.current_resources()

    def get(self, remote_path, local_path):
        netutils.get_file(remote_path, local_path, self.credentials,
                          self.xfer_mthd)

    def ls(self, path):
        pass

    def make_task(self, cmd_template, arg_dict, resources):
        tname = "task_" + utility.generate_time_string() + str(self.task_count)
        self.task_count += 1
        tdir = self.sandbox + tname + "/"
        return Task(tname, cmd_template, arg_dict,
                 tdir, resources)

    def mkdir(self, path):
        pass

    def put(self, local_path, remote_path):
        netutils.put_file(local_path, remote_path, self.credentials,
                          self.xfer_mthd)

    def remove(self, path):
        pass

    def request_resources(self, res):
        resource_server = self._res_connect()
        return(resource_server.request(res, self.host))
    
    def relinquish_resources(self, res):
        resource_server = self._res_connect()
        return(resource_server.relinquish(res, self.host))
            
    def start_task(self, task):
        pass

    def task_status(self, task):
        pass

    def stop_task(self, task):
        pass


class SingleMachine(Location):

    def __init__(self,
                 host,
                 sandbox,
                 cred_file=None,
                 xfer_mthd='sftp'):
        Location.__init__(self, host, sandbox,
                          cred_file, xfer_mthd)

    def _make_exec_script(self, task):
        fn = task.id + ".sh"
        fh = open(fn, 'w')
        j = "".join
        #fh.write(j(["mkdir ", task.working_dir, "\n"]))
        fh.write(j(["cd ", task.working_dir, "\n"]))
        fh.write(j(["nohup ", task.cmd_str, " > ", task.id, ".out",
                    " &", "\n"]))
        fh.write(j(["echo", " $!", "\n"]))
        
        fh.close()

        return fn

    def start_task(self, task):
        conn = netutils.SSHConnection(self.credentials.host,
                                      self.credentials.user,
                                      self.credentials.passwd)
        #write a bash script
        script_path = self._make_exec_script(task)
        script_name = os.path.split(script_path)[1]
        remote_script = self.sandbox + script_name

        #make the task's working directory
        conn.execute("mkdir " + task.working_dir)

        #put and execute the script
        self.put(script_path, remote_script)
        conn.execute("chmod u+x " + remote_script)
        out = conn.execute("sh " + remote_script)
        #script only prints the process id of the task
        task.proc_id = out[0].strip()

        #remove the script file
        conn.execute("rm " + script_path)
        
        
    def stop_task(self, task):
        pass

    def task_status(self, task):
        pass


class Task:

    def __init__(self,
                 iden,
                 cmd_template,
                 arg_dict,
                 working_dir,
                 resources,
                 outputs = set()):
        self.cmd_template = cmd_template
        self.cmd_str = ""
        self.arg_dict = arg_dict
        self.id = iden
        self.proc_id = None
        self.used_resources = resources
        self.working_dir = working_dir
        self.local_files = {}
        self.remote_paths = {}
        self.outputs = output_dict
        self._handle_argfiles()
        self._handle_outputs()

    def _handle_argfiles(self):
        '''
            Converts any local files contained in self.arg_dict
            into their working path counterparts. Updates the value
            in the dictionary for the file's flag to the working path.

            Adds the original file path to the self.local_files dict.
            
        '''
        for key, val in self.arg_dict.items():
            if(os.path.isfile(val)):
                self.local_files[key] = val
                d, fn = os.path.split(val)                
                #convert the path to a remote version
                self.arg_dict[key] = self.working_dir + fn
                    
    def build_cmd(self):
        self.cmd_str = self.cmd_template
        for key, val in self.arg_dict.items():
            self.cmd_str = self.cmd_str.replace("<" + key + ">",
                                                str(value))
        
    
