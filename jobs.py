'''


Author: Mathew Citarella

jobs.py: Holds classes and methods for creating and manipulating autonomics jobs

'''

from sqlalchemy.sql import functions, select, and_, select, or_
from multiprocessing import Process
from autonomics.file_io import AlignmentReader, FileExtensions
from autonomics.file_io import translate_seq_file, make_record, sha1_file
from autonomics import statistics
from autonomics.statistics import do_stats_update
from autonomics.utility import attr_or_default, die_on_error
from autonomics import settings, netutils
from utility import convert_if_int
import datetime
import imaplib
import os
import pickle
import shutil
import ssl
import subprocess
import sys
import threading
import time
import glob
import redis

session = netutils.make_db_session()


JOB_INPUT_FLAGS = {
    'adapter_trim': {'input': '.fastq', 'end2': '.fastq.end2'},
    'quality_trim': {'input': '.fastq', 'end2': '.fastq.end2'},
    'read_normalization': {'input': '.fastq',
                           'end2': '.fastq'},
    'quantification': {'query': '.fastq',
                       'query2': '.fastq.end2',
                       'db': '_project.fasta'},
    'assemble': {'input': '.fasta',
                 'end1': '.fastq',
                 'end2': '.fastq.end2'},
    'blast_swissprot': {'input': '_project.fasta'},
    'blast_nr': {'input': '_project.fasta'},
    'pfam': {'input': '_project.fasta'},
    'go': {'input': '_blast_swissprot.txt'},
    'kegg': {'input': '_blast_swissprot.txt'},
    'panther': {'input': '.fasta'}
}


def generate_split_filename(fn, index, ext):
    ''' file (str):
            Base path to the file you are splitting.

        index (int):
            Index of this particular split file.

        ext (str):
            Extension of the original input file.

        Creates a split file name. Returns file + str(index) + ext.
    '''
    return fn + str(index) + ext


def get_input_size(infile):
    ''' infile (str):
            Path to the f for which you want to determine the input size.

        Returns the number of input elements in the specified infile. 
        Currently supports counting the number of input elements in fasta and 
        fastq files.
    '''
    f = open(infile, 'r')
    base, ext = os.path.splitext(infile)
    if(ext == '.fasta' or ext == '.fa'):
        count = 0
        for l in f:
            if(">" in l):
                count += 1
        return count
    elif(ext == '.fastq' or ".fastq" in base):
        count = 0
        lc = 0
        for l in f:
            lc += 1
            if(lc % 4 == 0):
                count += 1
        return count

    else:
        i = -1
        for i, l in enumerate(infile):
            pass
        return i + 1


def get_job_name(jid, special_run,  retries=5):
    ''' jid (int):
            Identifier for the job.

        retries (int):
            Number of times to attempt to get the job name before raising
             an error. Default is 5.

        Returns the job_name for the job specified by jid.
    '''
    global session

    try:
        jn_mapping = netutils.get_table_object("jn_mapping", session)
        pn_mapping = netutils.get_table_object("pn_mapping", session)
        q = jn_mapping.join(
                            pn_mapping,
                            jn_mapping.c.project_id==pn_mapping.c.project_id
                            ).select(jn_mapping.c.job_id==jid)
        results = q.execute()
        row = results.fetchone()
        results.close()
        if special_run:
            return row.job_name
        else:
            return row.project_name + "_" + row.job_type
    except:
        raise
        if(retries <= 0):
            raise
        session.close()
        session = netutils.make_db_session()
        return get_job_name(jid, special_run, retries - 1)


def get_project_name(pid, retries=5):
    '''
        Same functionality as get_job_name, but returns the project_name for
         the project given by pid instead.
    '''
    global session
    try:
        pn_mapping = netutils.get_table_object("pn_mapping", session)
        results = pn_mapping.select(pn_mapping.c.project_id==pid).execute()
        row = results.fetchone()
        results.close()
        return row.project_name
    except:
        raise
        if(retries <= 0):
            raise
        session.close()
        session = netutils.make_db_session()
        return get_project_name(pid, retries - 1)


def make_blast_index(db_path, db_type, blocking=True):
    ''' db_path (str):
            Path to the FASTA-formatted sequence file you wish to make into a
            BLAST index.

        db_type (str):
            Whether the database file contains amino acid or nucleotide
            sequences. NT for nucleotides, AA for amino acids.

        blocking (boolean):
            If True, the method waits for makeblastdb to return before
            returning. If False, this method returns immediately.

        Indexes a FASTA-formatted sequence file for use with the BLAST+
        suite of sequence alignment tools.
    '''
    if(db_type == "NT"):
        t_string = "nucl"
    else:
        t_string = "prot"
    proc = subprocess.Popen(
                            "makeblastdb -in " + db_path + " -dbtype " +
                            t_string, shell=True
                            )
    if(blocking):
        proc.wait()


def make_bowtie_index(db_path, out_name, blocking=True):
    ''' db_path (str):
            Path to the FASTA-formatted sequence file you wish to index for
            bowtie.

        out_name (str):
            Base name for the bowtie indexes.

        blocking (boolean):
            Same as in make_blast_index.

        Indexes a FASTA-formatted sequence file for use with the bowtie
        aligner. Only supports RNA/DNA database sequences.
    '''
    cmd = "bowtie-build " + db_path + " " + out_name + " >& /dev/null"
    proc = subprocess.Popen(cmd, shell=True)
    if(blocking):
        proc.wait()


def make_remote_dir(dir_name, c):
    ''' dir_name (str):
            Path to mkdir on.

        c (netutils.SSHConnection):
            Active ssh connection that will be used to connect and execute
            the mkdir.

        Executes mkdir on the remote server to which c is connected.
    '''

    c.execute("mkdir " + str(dir_name))


def next_split_file(base_name, ext, current_num):
    ''' base_nane (str):
            Base name of the file you are splitting.

        ext (str):
            Extension of the file you are splitting.

        current_num (int):
            Current index of the file you are splitting.

        Creates a new split filename by incrementing current_num and
        concatenating base_name, str(new_count), ext. Returns a tuple
        (new_name, new_count)
    '''
    new_count = current_num + 1
    new_name = base_name + "_" + str(new_count) + ext
    return (new_name, new_count)


def parse_resources(resource_str):
    ''' resource_str (str):
            The 'resources' field from a row in the args table, corresponding
             to the resources for a given job.

        Returns a dictionary of resource->amount pairs. Assumes the
        resource_str is in the format resource1:amount1, resource2:amount2,
        etc.

    '''
    if(len(resource_str) > 0):
        return {(element.split(":")[0],
                 convert_if_int(element.split(":")[1])
                 ) for element in resource_str.split(",")
                }
    else:
        return {}


def split_input(f, job_type, proc_units):
    ''' f (str):
            Path to the file you are going to split.

        job_type (str):
            Job type of the job requsting the split.

        proc_units(int):
            Number of processing units used to run the job requesting the
            split. This method will try to split the file into this many pieces.

        Returns a list of file names that correspond to the files produced by
        the split. Currently knows how to split fastq and fasta files.
    '''
    base, ext = os.path.splitext(f)
    if(ext == ".end2"):
        base, ext = os.path.splitext(base)
    if(ext == '.fasta' or ext =='.fa'):
        ret = split_fasta(f, job_type, proc_units)
    elif(ext == '.fastq'):
        ret = split_fastq(f, job_type, proc_units)

    return ret


def split_fastq(f, job_type, proc_units):
    '''
        Arguments are the same as split_input.

        Splits a FASTQ-formatted sequence file into proc_units number of files. Returns a list of names for those files.
    '''
    fh = open(f, 'r')
    d, fname = os.path.split(f)
    base, ext = os.path.splitext(fname)
    seqs_per = (get_input_size(f)/proc_units) + 1
    seq_count = 0
    file_count = 1
    filenames = []
    base = base + "_" + job_type
    curfile = base + "_" + str(file_count) + ext
    curpath = d + "/" + curfile
    fw = open(curpath, 'w')
    filenames.append(curpath)
    lc = 0
    for line in fh:
        lc += 1
        line = line.strip()
        fw.write(line + "\n")
        if(lc % 4 == 0):
            seq_count += 1
        if(seq_count == seqs_per):
            fw.close()
            curfile, file_count = next_split_file(base, ext, file_count)
            curpath = d + "/" + curfile
            fw = open(curpath, 'w')
            filenames.append(curpath)
            seq_count = 0

    return filenames


def split_fasta(q_file, job_type, proc_units):
    '''
        Arguments are the same as split_input.

        Splits a FASTA-formatted file into proc_units number of files. Returns a list of names for those files.
    '''
    fh = open(q_file, 'r')
    directory, file = os.path.split(q_file)
    file_base, file_ext = os.path.splitext(file)
    seqs_per_file = (get_input_size(q_file)/proc_units) + 1
    seq_counter = 0
    file_counter = 1
    filenames = []
    file_base += "_" + job_type
    cur_file = file_base + "_" + str(file_counter) + file_ext
    curpath = (directory + "/" + cur_file)
    fw = open(curpath, 'w')
    filenames.append(curpath)
    for line in fh:
        line = line.strip()
        if(">" in line):
            seq_counter += 1
            if(seq_counter > int(seqs_per_file)):
                #close the open partition file
                fw.close()
                #increment the file counter and open the next file
                file_counter += 1
                cur_file = file_base + "_" + str(file_counter) + file_ext
                curpath = directory + "/" + cur_file
                fw = open(curpath, 'w')
                filenames.append(curpath)
                #reset the seq_counter to 1
                seq_counter = 1

        fw.write(line + "\n")

    fw.close()
    return filenames


class Arguments:
    '''
        Description: This class is a container for arguments passed to a job through the system, either as arguments to the pipeline (pipeline_args), or to a process the job is going to spawn (process_args).
    '''

    def __init__(self, arg_string, job_type=None):
        self.has_required = set()
        self.arg_string = arg_string
        self.positional = []
        self.job_type = job_type
        self._mark_required(self.arg_class)
        self._init_attributes(self.arg_class)


    def parse(self):
        ''' This method parses an argument string that was passed
            to Arguments(). Attributes are assigned to the object
            with names matching the flags given in the arg string.

            example:
            arg string '--query foo -db bar foobar'
            Will be parsed into an Arguments() object with the
            following attributes:

                obj.query = 'foo'
                obj.db = 'bar'
                obj.positional = ['foobar']

        '''

        class ParseState:
            ''' This class acts as an enum, providing valid states for the
                parser.

                states:
                    SEEK: 1
                    SINGLE_DASH: 2
                    DOUBLE_DASH: 3
                    ARG_NAME: 4
                    VALUE: 5
                    AFTER_ARG: 6
                    POSITIONAL: 7

            '''

            SEEK = 1
            SINGLE_DASH = 2
            DOUBLE_DASH = 3
            ARG_NAME = 4
            VALUE = 5
            AFTER_ARG = 6
            POSITIONAL = 7

        state = ParseState.SEEK
        arg_name = ""
        value = ""
        quote_open = False
        has_required_value = False

        if(len(self.arg_string) == 0):
            return

        for char in self.arg_string:
            if(state == ParseState.SEEK):
                if(char == '-'):
                    state = ParseState.SINGLE_DASH
                elif(char == '"'):
                    state = ParseState.POSITIONAL
                    value = char
                    quote_open = True
                elif(char != " "):
                    value = char
                    state = ParseState.POSITIONAL

            elif(state == ParseState.POSITIONAL):
                if(char == '"'):
                    quote_open = not quote_open
                elif(char == ' '):
                    if(quote_open):
                        value += char
                    else:
                        self.positional.append(value)
                        arg_name, value = self._reset_names()
                        state = ParseState.SEEK
                else:
                    value += char

            elif(state == ParseState.SINGLE_DASH):
                if(char == '-'):
                    state = ParseState.DOUBLE_DASH
                elif(char == ' '):
                    state = ParseState.SEEK
                else:
                    state = ParseState.ARG_NAME
                    arg_name += char

            elif(state == ParseState.DOUBLE_DASH):
                if(char == ' '):
                    state = ParseState.SEEK
                else:
                    state = ParseState.ARG_NAME
                    arg_name += char

            elif(state == ParseState.ARG_NAME):
                if(char == ' '):
                    arg_name = arg_name.replace("-", "_")
                    has_required_value = self.is_value_required(arg_name)
                    state = ParseState.AFTER_ARG

                else:
                    #keep parsing until there's a space
                    arg_name += char

            elif(state == ParseState.AFTER_ARG):
                if(char == ' '):
                    self._setattr(arg_name, True)
                    arg_name, value = self._reset_names()
                    state = ParseState.SEEK
                elif(char == '-'):
                    self._setattr(arg_name, True)
                    arg_name, value = self._reset_names()
                    state = ParseState.SINGLE_DASH
                else:
                    if(has_required_value):
                        if(char == '"'):
                            quote_open = True
                            state = ParseState.VALUE
                        else:
                            value = char
                            state = ParseState.VALUE
                    else:
                        self._setattr(arg_name, value)
                        state = ParseState.POSITIONAL
                        if(char == '"'):
                            quote_open = True
                        else:
                            value = char

            elif(state == ParseState.VALUE):
                if(char == ' '):
                    if(quote_open):
                        value += char
                    else:
                        self._setattr(arg_name, value)
                        arg_name, value = self._reset_names()
                        state = ParseState.SEEK
                elif(char == '"'):
                    quote_open = not quote_open
                else:
                    value += char

        arg_name = arg_name.replace("-", "_")
        if(state == ParseState.POSITIONAL):
            self.positional.append(value)
        elif(state != ParseState.VALUE):
            self._setattr(arg_name, True)
        else:
            self._setattr(arg_name, value)

    def _init_attributes(self, from_table, overwrite=False, retries=5):
        ''' from_table (str):
                Table to use when checking for default flag values. Currently
                accepted values are: 'proc_options', 'pipe_options'

            overwrite (boolean):
                Do you want to overwrite values already assigned to this object
                with default values? Defaults to no.

            Looks up default arguments in from_table, and assigns their values
            to this arguments object. If overwrite is True, will overwrite
            existing attributes with default values from the tables.
            here, get all of the supported options and set them to the default 
            values, if not already set
        '''

        session = netutils.make_db_session()
        options = netutils.get_table_object(from_table, session)
        res = session.conn.execute(options.select(options.c.job_type==self.job_type))
        value_map = {'False': False, 'True': True}
        for row in res.fetchall():
            if(hasattr(self, row.flag) and not overwrite):
                continue
            if(row.arg_required=='Y'):
                setattr(self, row.flag, row.default_value)
            else:
                if(row.default_value in value_map):
                    setattr(self, row.flag, value_map[row.default_value])
                else:
                    setattr(self, row.flag, False)

    def _mark_required(self, table_name):
        '''
            Marks flags as having required args, based on whether or not they
            are said to have arguments in the proc_options or pipe_options
            table.

            Uses self.job_type to select the default flags supported by the job
            type these arguments are attached to.
        '''

        session = netutils.make_db_session()

        if(self.job_type is None):
            return
        table = netutils.get_table_object(table_name, session)
        results = table.select(and_(table.c.job_type==self.job_type,
                                    table.c.arg_required=='Y')).execute()
        for row in results.fetchall():
            self.mark_value_required(row.flag)

    def _reset_names(self):
        '''
            Returns a two-tuple of two empty strings.
        '''
        return ("", "")

    def _setattr(self, name, value):
        '''
            name (str):
                Name of the attribute you are setting.

            value (obj):
                Value to assign to the attribute specified by name.

                If value is an empty string, assigns a value of True, as the 
                flag being assigned does not have a value.
            
        '''
        if(value == ''):
            setattr(self, name, True)
        else:
            setattr(self, name, value)

    def lmark_value_required(self, list_of_names):
        ''' Marks all flag names in list_of_names as requiring an argument.
        '''
        self.has_required.update(list_of_names)

    def mark_value_required(self, name):
        ''' name (str):
                The name of the flag  to mark as having a required value.

            Marks the flag specified by name as requiring an argument. Used
            when parsing the arguments used to construct this object.
        '''
        self.has_required.add(name)

    def is_value_required(self, name):
        if(name in self.has_required):
            return True
        else:
            return False

    def print_me(self):
        print(self.positional)
        for key, value in self.__dict__.items():
            print(str(key) + " " + str(value))


class PipeArgs(Arguments):
    ''' Subclass of Arguments that sets up default Pipeline Args for the
        instantiating job.
    '''

    arg_class = 'pipe_options'

    def __init__(self, arg_string, job_type):
        Arguments.__init__(self, arg_string, job_type)


class ProcArgs(Arguments):
    ''' Subclass of Arguments that sets up default Process Args for the
        instantiating job.
    '''

    arg_class = 'proc_options'

    def __init__(self,arg_string, job_type):
        Arguments.__init__(self, arg_string, job_type)


class Resources:
    ''' Represents resources used/free at a given location.

        Attributes:
            totals: a dictionary of resource -> total amount
            free: a dictionary of resource-> free amount

        Manager.py uses this class to determine how many total resources
        are free at a given location known to the system.

        In a future release, this object will be replaced by a stand-alone
        resource server holding global information about resource levels.
  
        amount of resources used for a job found in args/default_args table

job_type     | executable | loc   | process_args        | pipeline_args                 | resources | priority |
adapter_trim              | local |                     | --adapter-trimmer cutadapt    | cpu:1     |      100 |
quality_trim              | local | --quality-cutoff 20 | --trimmer cutadapt            | cpu:1     |      100 |
read_normalization        | local |                     |                               | cpu:1     |      100 |
assemble                  | local |                     | --assembler trinity           | cpu:20    |        1 |
go                        | local |                     |                               | cpu:1     |        1 |
kegg                      | local |                     |                               | cpu:1     |        1 |
pfam         | NA         | HPC   | -fasta <input>      | --translate                   | cpu:150   |        1 |
quantifica   | bowtie     | HPC   | --query <query> --db <db> | --aligner bowtie --db-type NT | cpu:50 |   100 |
blast_nr     | blastx     | HPC   | -query <input> -evalue 1e-04 -num_alignments 5 -num_descriptions 5 -out <output> |   | cpu:200 | 1 |
blast_swiss  | blastx     | HPC   | -query <input> -evalue 1e-04 -num_alignments 5 -num_descriptions 5 -out <output> |   | cpu:100 | 1 |

    '''


    def __init__(self):

        self.totals = {}
        self.free = {}

    def add_resource(self, res_type, res_max_amt):
        ''' res_type (str):
                The type of resource to update.

            res_max_amt (int):
                The total amount of the resource available at the location
                represented by this Resource object.

            Adds a new resource 'res_type' with total available 'res_max_amt'
            to this Resource object.
        '''
        self.totals[res_type] = res_max_amt;
        self.free[res_type] = res_max_amt;

    def give_to(self, job):
        '''
            job (Job):
                A subclass of the jobs.Job class.

            Checks the required resources of job. If this object has enough
            resources to give, remove the resources from this object and
            returns True.
        '''
        enoughResources = self.has_enough_free(job.resources)

        if(enoughResources):
            for key in job.resources:
                if(key in self.totals):
                    self.free[key] -= job.resources[key]

        return enoughResources

    def has_enough_free(self, resources):
        '''
            resources (dict{res:amount}):
                A dictionary of resources to check against the available
                resources in this object.

            Checks if this object has enough of all resources in the supplied
            dictionary, if so returns True. Returns False otherwise.
        '''
        enough = True
        for key, value in resources.items():
            if(key in self.totals):
                if(value > self.free[key]):
                    enough = False
        return enough

    def take_from(self, job):
        '''
            job (jobs.Job):
                Subclass of jobs.Job from which resources are collected.

            Returns all resources used by the job to the pool of free resources
            for this object, using the jobs.resource attribute.
        '''
        resources = [key for key in job.resources.iterkeys() if key in self.totals]
        for key in resources:
            if(self.free[key] + job.resources[key] < self.totals[key]):
                self.free[key] += job.resources[key]
            else:
                self.free[key] = self.totals[key]

class Locations:

    LOCAL = 1
    HPC = 2


class JobState:

    RUNNING = 0
    FINISHED = 1
    ERROR = 2


class Job:
    ''' This class is the ancestor of all analysis tasks.

        The main responsibility of the class is to provide the manager with
        an interface for starting, managing, and completing various analyses
        across the system.

        As such, it has the following attributes:

            - pid (int):
                The ID of the project this job belongs to.

            - pn (str):
                The name of the project this job belongs to.

            - jid (int):
                The unique job ID for this analysis task.

            - job_name (str):
                The name of this analysis task.
                default: project_name + " " + job_name

            - location (jobs.Locations):
                One of Locations attributes, specifying where this job should
                be executed.

            - files ([]):
                A list of file paths for files associated with this job.

            - processes ([]):
                A list of Popen objects that constitute the subprocesses
                this job uses to complete its work. A mixture of local and
                remote processes.

            - executable (str):
                The name of the executable that this job is currently using
                to perform analysis tasks. Can be changed throughout the
                life of a job.

                HPCJobs use this value to select a command template from the
                executables table.

            - local_dir (str):
                The working directory for this job. All intermediary files
                and output are placed in this directory.

                Defaults to: settings.proj_dir + pn

                If this is a command-line submitted job, this attribute will
                be overridden by: PATH_TO_COMMAND_LINE_JOBS/job_name/

            - run_name (str):
                The unique name for this particular instance of the analyis
                task represented by the job_type

                For instance, a project may have more than one BLAST annotation
                running in paralell. The run_name is used to name temporary
                files and remote directories to avoid name collisions that
                would occur by simply using the job_name attribute.

                Currently generated by taking the local system time and
                appending the job name to it.

            - resources ({resource_name: resource_amt}}:
                A dictionary of resources required to run this job. Initially
                created by a call to the _default_resources() method.

                Each subclass of Job is expected to override _default_resources
                in order to specify the resources it would like to request
                in the default run scenario.

            - input_files ({file_flag: file_path}):
                A dictionary of input files used by this job class. file_flag
                is the name of the flag that is accepted in the pipeline_args
                field of the args table to specify this particular input file.

                Each job has default input files specified in the
                JOB_INPUT_FLAGS dictionary. If flags for input are not found
                in pipeline_args for a given job, the default values in
                JOB_INPUT_FLAGS are used to populate input_files instead.

            - output_files (set):
                A set of the paths to all output files produced over the
                lifetime of the job.

            - current_output (str):
                The full path to the file currently being written to by this
                job's active analysis process.

            - DO_NOT_SPLIT ([str]):
                A list of file_flags which should not be split when submitting
                concurrent analyses to remote computational locations.
    '''

    def __init__(self, pid, jid, t, executable="NA", resources="",
                 pipeline_arg_str='', process_arg_str=''):
        self.job_type = t
        self.pid = pid
        self.jid = jid
        self.pn = get_project_name(self.pid)
        self.special_run = 0
        if self.pn == 'user_job':
            self.special_run = 1;
        self.job_name = get_job_name(self.jid, self.special_run)
        self.location = Locations.LOCAL
        self.errors = []
        self.files = []
        self.processes = []
        self.executable = executable
        self.local_dir = settings.proj_dir + "/" + self.pn + "/"
        if self.special_run:
            self.local_dir = settings.proj_dir + "/special_runs/" + \
                            self.job_name + "/"
        self.run_name = self.generate_run_name()
        self.resources = self._default_resources()
        self.resources.update(parse_resources(resources))
        self.input_file = ""
        #A dict of flag->input file pairs, where the input file is either automatically detected or
        #provided via a flag in pipeline_args. Each job will provide its own dict of default_in_files.
        self.input_files = {}
        #each output file created during job execution should be added to this list
        self.output_files = set()
        self.DO_NOT_SPLIT = []
        self.outputFile = ""
        out_modifier = self.job_name + ".txt"
        if(t in FileExtensions.output_exts):
            out_modifier = self.pn + FileExtensions.for_output(self.job_type)
        self.current_output = self.local_dir + out_modifier
        self.pipeline_args = PipeArgs(pipeline_arg_str, self.job_type)
        self.process_args = ProcArgs(process_arg_str, self.job_type)
        self._handle_args()
        self.input_files = self._determine_input()
        self.current_output = self._determine_output()
        self.output_files.add(self.current_output)

    def _determine_input(self):
        '''
            Interates over the input_flags set in JOB_INPUT_FLAGS for this
            job_type. For each found flag, it checks to see if that flag was
            set in the pipeline_args field of the args table for this job.

            If the flag was found, input_files[input_flag] is set to the
            value in pipeline_args. Otherwise the default value from
            JOB_INPUT_FLAGS is used.

            Each job type must have an entry in JOB_INPUT_FLAGS to
            automatically determine input during class instantiation. Jobs
            that do not do so must handle their input within the run() method
            of the job process class.

            Returns a dictionary of input files for this job in the form:
            {input_flag: path_to_input}

        '''
        #get the list of required input files for this job type
        ret = {}
        if(self.job_type in  JOB_INPUT_FLAGS):
            for flag, extension in JOB_INPUT_FLAGS[self.job_type].items():
                ret[flag] = self.local_dir + self.pn + extension
                if(hasattr(self.pipeline_args, flag)):
                    if(os.path.isfile(getattr(self.pipeline_args, flag))):
                        ret[flag] = getattr(self.pipeline_args, flag)

        return ret

    def _determine_output(self):
        '''
            Checks to see if the --output-file flag was set in the
            pipeline_args for this job. If so, returns the value of that flag.
            Otherwise returns self.current_output.
        '''
        out_flag = 'output_file'
        if(hasattr(self.pipeline_args, out_flag)):
            return getattr(self.pipeline_args, out_flag)

        return self.current_output

    def _handle_args(self):
        '''
            Provides custom argument handling. Each job should override this
            method depending on its own needs for using the flags set in
            pipeline_args and process_args. By default, marks that all
            default input file flags, ouput file flags, and job name flags
            have a required value. Additionally determines if this job is
            using paired-end input.
        '''
        command_line_inputs = ['input_file', 'output_file', 'job_name']
        self.pipeline_args.lmark_value_required(command_line_inputs)
        default_inputs = JOB_INPUT_FLAGS[self.job_type]
        self.pipeline_args.lmark_value_required(default_inputs.keys())
        self.pipeline_args.parse()
        self.process_args.parse()

        #see if we need to get defaults
        if(hasattr(self.pipeline_args, 'paired_end')):
            self.paired_end = self.pipeline_args.paired_end
        else:
            self.paired_end = False

    def _input_exists(self):
        '''
            Determines if the file path given in self.input_files['input']
            exists. Returns True if so, False otherwise.
        '''
        return os.path.exists(self.input_files['input'])

    def _default_resources(self):
        '''
            Returns a dictionary of default resources this job class requests
            before being run by the manager process. Each class should override
            this class to provide custom resource requirements.

            Default:
                {'cpu': 1}
        '''
        return {"cpu": 1}

    def _exec_cmd(self, cmd, blocking=True, capture_output=False):
        '''
            Executes the command given by 'cmd' locally. If blocking is True
            the thread calling _exec_cmd waits for the command to exit. If
            capture_output is True then stdout and stderr are captured to
            internal variables and will not show in the console running
            the system.
        '''
        redirect = None
        if(capture_output):
            redirect = subprocess.PIPE
        p = subprocess.Popen(cmd, shell=True, stdout=redirect, stderr=redirect)
        if(blocking):
            p.wait()
        return p

    def _rm_working_dir(self):
        pass

    def _set_resources(self, resource_str):
        '''
            Given the string resource_str, parses that string and sets
            self.resources equal to rhe dictionary returned.

            Deprecated, and will be removed in an upcoming version.
        '''
        self.resources = self.parse_resources(resource_str)

    def check(self):
        '''
            Iterates over the list of objects in self.processes and determines
            if they are still running.

            If one or more processes have exited with an error status, returns
            JobState.ERROR

            If the call to is_alive() on at least one process returns True,
            returns JobState.RUNNING

            If no process objects fall into either of the above criteria,
            returns JobState.FINISHED

        '''
        for p in self.processes:
            if(p.is_alive()):
                return JobState.RUNNING
            elif(p.exitcode != 0):
                print(p.exitcode)
                print"Got an unexpected exit code.  for jid: ", self.jid
                return JobState.ERROR
            elif(self.error_status()):
                print "Got an error status for jid: ", self.jid
                print(self.error_status())
                return JobState.ERROR

        return JobState.FINISHED

    def cleanup(self):
        '''
            Method tasked with removing undeeded intermediatry files and
            cleaning up the database tables after job execution. Each subclass
            should override this method.
        '''
        pass

    def complete(self):
        '''
            Marks the job as complete in the Autonomics database. Job status
            is held in the jobs table. 
            
            By default, this method also cleans up intermediary and temporary
            files and will upload the output of this job to a data destination
            if upload is confingured for this job in the system. This method
            will also calculate statistics for the analysis task, if possible.
            Statistics are stored in the run_stats table.  

            Currently, the job_type 'upload' is not configured. In the past it
            was used to copy finished results to our only neurobase server.  Now
            we have multiple neurobase servers on different machines.  After 
            confirming the finished results are valid and complete, we copy the 
            results by hand (rsync) to the appropriate neurbase server.

        '''
        session = netutils.make_db_session()

        results = session.conn.execute(
            "select project_id from jn_mapping where ((project_id = ?) and \
            (job_type = 'upload'))",
            (self.pid))
        xxx = results.fetchone()
        if(not xxx is None): # then we upload this job
            self.upload()

        self.cleanup()
        self.mark_complete()

    def error_status(self):
        '''
            Returns:
                True if this job did not produce its specified output 
                files, False otherwise. 
        '''
        for f in self.output_files:
            if(not os.path.isfile(f)):
                return True

        return False

    def generate_run_name(self):
        '''
            Generates a run name for use in naming temporary files and folders.
            
            Returns:
                self.job_name + str(t.year) + str(t.month) + str(t.day) + \
                str(t.hour) + str(t.minute) + str(t.second)
                
                Where t is set to datetime.datetime.now().
        '''
        t = datetime.datetime.now()
        name = self.job_name + str(t.year) + str(t.month) + str(t.day) + \
                str(t.hour) + str(t.minute) + str(t.second)
        return name.strip()

    def mark_complete(self):
        '''
            Marks this job as complete in the jobs table. This is currently
            accomplished by setting the finished field = 'Y' and f_ts =
            CURRENT_TIMESTAMP().
        '''
        session = netutils.make_db_session()
        jn = netutils.get_table_object('jn_mapping', session)
        jn.update(
                  ).where(
                          jn.c.job_id==self.jid
                          ).values(
                                   f_ts=functions.current_timestamp(),
                                   finished='Y'
                                   ).execute()
        if self.special_run:
            d = netutils.get_table_object('quenew_special', session)
            d.delete().where(d.c.job_id==self.jid).execute()
        else:
            d = netutils.get_table_object('quenew', session)
            d.delete().where(d.c.job_id==self.jid).execute()

    def replace_in_proc_args(self, placeholder, repl):
        '''
            Deprecated by the switch to the command template system.
        '''
        self.process_args.arg_string = self.process_args.arg_string.\
            replace("<" + placeholder + ">" , repl)

    def replace_in_pipe_args(self, placeholder, repl):
        '''
            Same as replace_in_proc_args.
        '''
        self.pipeline_args.arg_string = self.pipeline_args.arg_string.\
            replace("<" + placeholder + ">" , repl)


class HPCJob(Job):
    '''
        Class that represents analyses run on an HPC cluster utilizing
        a PBS system. Extends the base Job class.
        
        This job class performs the following actions:
            
            1) Determines input files for the job
            2) Parses pipeline and process arguments
            3) Gets a command template for the executable that is to be run
            4) Converts local paths for all input and auxiliary files to 
               remote paths
            5) Creates a remote working directory
            6) Moves auxiliary files to the working directory
            7) Splits the input files, one file created for each CPU (node) 
               requested
            8) For each concurrent analysis:
                - creates a fully-qualified command by replacing process_args
                  and split input files into the command template
                - creates a qsub script to call the fully-qualifed command on 
                  the split input file
                - moves the individual input file and the qsub script to the 
                  remote machine
                - executes the qsub script to submit the analysis to the HPC 
                  cluster
            9) Monitors the concurrent analyses for this job, retreiving
               output as it is created
            10) Resubmits failed jobs
            11) Marks the job as complete when all concurrent analysis tasks 
                have successfully completed
            
        
        Additional attributes:
            hpc_job_files ([str]):
                A list of file paths. The files in this list are necessary
                for the remote analysis task to perform its work, but are 
                not input files and thus should not be split. 
                
                An example would be BLAST databases that need to be moved 
                to HPC to perform the alignment.
                
            hpc_qsub files ([str]):
                A list of file paths. The files in this list are input files
                to the remote analysis task. These files will be split into
                individual inputs to the concurrent analysis task at HPC.
                
            mail_credentials (settings.Credentials):
                A credentials object used to store the means of connecting
                to an email account. This account should be set to receive
                messages from the high-performance cluster about the status
                of running/failed/aborted jobs.
                
            ssh_credentials (settings.Credentials):
                A credentials object used to connect to the submission node
                for the HPC cluster via ssh.
                
            mem_increment (int)=200:
                Value (in megabytes) used to increment memory requests when 
                resubmitting jobs that fail to due to exceeding their initial
                 memory allocation.
                
            retry_interval (int)=30:
                How long to wait (in seconds) between attempts to perform
                an operation via ssh.
            
            restart_jobs (Boolean)=True:
                Whether or not jobs that fail or are aborted should be 
                automatically resubmitted.
                
            remote_dir (str):
                Path to the working directory for this job on the HPC cluster.
                
                Default:
                    settings.remote_dir + self.run_name
                    
            hpc_command (str):
                Command string to be executed on the HPC cluster.
                
                This command is executed once per requested CPU in this job's
                resource dictionary.
                
            constructor (PipeProcess):
                A PipeProcess subclass that is used to fork a new process
                in which the actual analysis starting and monitoring takes 
                place, freeing the calling thread to continue starting 
                additional analyses.
    '''

    mem_increment = 200
    retry_interval = 30
    retries = 5
    restart_jobs = settings.restart_jobs
    location = Locations.HPC

    def __init__(self, pid, jid, job_type, executable, resources="",
                 pipeline_arg_str="", process_arg_str=""):
        Job.__init__(self, pid, jid, job_type, executable, resources, 
                     pipeline_arg_str, process_arg_str)
        self.executable = executable
        self.hpc_job_files = {}
        self.hpc_qsub_files = {}
        self.processes = []
        self.mail_credentials = settings.mail_cred
        self.ssh_credentials = settings.hpc_cred
        self.location = HPCJob.location
        self.mail_connect = None
        self.remote_connect = None
        self.mail_enabled = True
        self.jid = jid
        self.constructor = HPCProcess
        self._resourceStr = resources
        self.remote_dir = settings.remote_dir + "/" + self.run_name + "/"
        self.hpc_command = None
        self._handle_args()

    def _default_resources(self):
        '''
            Returns a dictionary of default resources for this job. These
            values are overridden with resources specified in the args table 
            for this particular analysis task.
            
            Returns:
                {"mem": "500mb", "wall": "24:00:00", "ppn": 1, "cpu": 1, 
                "nodesPerSubJob": 1, "queue": "bio"}
        '''
        return {"mem": "500mb", "wall": "24:00:00", "ppn": 1, "cpu": 1, 
                "nodesPerSubJob": 1, "queue": "bio"}

    def start(self):
        '''
            Instantiates a subclass of PipeProcess and calls its start() 
            method. The PipeProcess forks a system process which handles
            the splitting of input files, command creation, file transfer,
            analysis monitoring, and output retrieval.
        '''
        self.p = self.constructor(self.input_files, 
                                  self.hpc_job_files, 
                                  self.hpc_qsub_files, 
                                  self.job_type, 
                                  self.resources, 
                                  self.remote_dir, 
                                  self.job_name, 
                                  self.local_dir, 
                                  self.executable, 
                                  self.process_args, 
                                  self.hpc_command, 
                                  self.current_output, 
                                  self.special_run, 
                                  self.jid, 
                                  self.outputFile)
        self.p.start()
        self.processes.append(self.p)

    def stop(self):pass


class LocalJob(Job):
    '''
        Subclass of Job providing basic functionality to start analyses local
        with to the Autonomics system server. Specifically, jobs subclassing
        this class will run on the same machine that manager.py is running.
    '''

    def __init__(self,
                 pid,
                 jid,
                 t,
                 executable="NA",
                 resources="",
                 pipeline_arg_str="",
                 process_arg_str=""):
        Job.__init__(self,
                     pid,
                     jid,
                     t,
                     executable=executable,
                     resources=resources,
                     pipeline_arg_str=pipeline_arg_str,
                     process_arg_str=process_arg_str)
        self.location = Locations.LOCAL

    def stop(self):
        for p in self.processes:
            p.terminate()

    def complete(self):
        Job.complete(self)


class AdapterTrimJob(LocalJob):
    '''
        Job class representing the task of trimming adapters from DNA or RNA
        sequences. Adapters are trimmed from reads using an already-installed
        third party package. 
        
        The system currently supports adapter trimming with either cutadapt or
        fastx_clipper.
        
        All adapters known to the system are contained with the known_adapters
        table. Adapters to be removed from a given project are assigned in the
        project_adapters table.
        
        For a given project, this job will iterate over all of the adapters
        assigned to that project, removing each one in turn.
        
        Process Constructor:
            AdapterTrimProcess
            
        input flags:
            input - A single FASTQ file, or the forward file from paired-end
                sequencing
            end2 (optional) - A FASTQ file representing the reverse file from
                paired-end sequencing
                      
        output:
            Files (with the same name as the input) that have the 
            adapter sequences removed from each input sequence.
            
        
    '''

    def __init__(self,
                 pid,
                 jid,
                 jtype,
                 executable,
                 resources,
                 pipeline_args,
                 process_args
                 ):
        LocalJob.__init__(self,
                          pid,
                          jid,
                          jtype,
                          executable,
                          resources,
                          pipeline_args,
                          process_args
                          )
    
        self.constructor = AdapterTrimProcess
        if(self.paired_end):
            self.output_files.add(self.input_files['end2'])

    def complete(self):
        self.cleanup()
        self.mark_complete()

    def start(self):
        p = self.constructor(self.pn,
                             self.input_files,
                             self.current_output,
                             self.pipeline_args,
                             self.process_args,
                             self.resources,
                             paired_end=self.paired_end)
        p.start()
        self.processes.append(p)


class AssemblyJob(LocalJob):
    '''
        This job class represents the task of locally assembling raw sequence
        into transcripts.
        
        The class supports assembly using two different assemblers: 
            mira Useful for low-coverage sequencing projects, including 
                sequencing performed on the 454 or Ion Torrent platforms
            trinity: More suited for high-coverage projects, including those
                produced by the Illumina HiSeq, Illumina MiSeq, and Ion Proton
                
            The assembler can be specified by setting the --assembler flag v
            in the pipeline_args field in the args table.
            
        Default resource request:
            cpu -> 20
            
        Constructor:
            AssemblyProcess
        
        input_flags:
            input - FASTQ file containing the forward reads
            end2 - FASTQ file containing the reverse reads (optional)
            
        output:
            FASTA file of assembled transcript + transcript fragments. If the
            assembler used produces quantification data (MIRA only), this
            job also produces the project quantification file.
    '''

    def __init__(self,
                 pid,
                 jid,
                 jtype,
                 executable,
                 resources,
                 pipeline_args,
                 process_args):

        LocalJob.__init__(self,
                          pid,
                          jid,
                          jtype,
                          executable,
                          resources,
                          pipeline_args,
                          process_args)

        self.constructor = AssemblyProcess

    def _default_resources(self):
        return {"cpu": 20}

    def start(self):
        '''
            Determines the assembler type from the assembler attribute of 
            this object's pipeline args. 
            
            Creates an AssemblyProcess object and calls its start method, 
            before appending the AssemblyProcecss to its process list 
            (self.processes).
        '''
        self.pipeline_args.parse()
        if (self.pipeline_args.assembler == "trinity" or
            self.pipeline_args.assembler == "mira" or
            self.pipeline_args.assembler == "miranewbler"):
            self.assembler = self.pipeline_args.assembler

            paired_end = 0
            if(self.paired_end):
                paired_end = 1
            p = self.constructor(self.pn,
                                 self.input_files['input'],
                                 self.current_output,
                                 self.process_args,
                                 paired_end,
                                 self.assembler,
                                 self.resources['cpu'],
                                 self.special_run,
                                 self.job_name
                                 )
            p.start()
            self.processes.append(p)

    def complete(self):
        '''
            Adds a file of unused reads to the set of output files for this 
            job, if such a file is present.
            
            Calculates assembly statistics for the overall assembly and stores
            them in the run_stats table. 
            
            If the assembler used creates a quantification file, this method
            adds that file to the set of output for this job and marks any
            downstream quantification jobs for this project as complete.
            
            The method then calls Job.complete().
        '''
        unused = self.local_dir + self.pn + "_unused_reads.fasta"
        if self.special_run:
            unused = self.local_dir + self.job_name + "_unused_reads.fasta"
        if(os.path.isfile(unused)):
            self.output_files.add(unused)
        fn = self.current_output
        rdict = statistics.CalculateAssemblyStatistics(fn)
        statistics.do_stats_update(rdict, self.pid,
                                    self.special_run, self.local_dir
                                    )
        if not self.special_run:
            fn = self.input_files['input']
            fn = fn.replace(".fasta", ".fastq")
            paired = 0
            if(self.input_files.has_key('end2') and
               os.path.exists(self.input_files['end2'])):
                paired = 1
                stats_temp_file = self.local_dir + 'stats.temp'
                cmd = 'cat ' + fn + ' ' + self.input_files['end2'] + \
                     ' > ' + stats_temp_file
                os.system(cmd)
                stats_temp_file = stats_temp_file.replace(".fasta", ".fastq")
                fn = stats_temp_file
            rdict = statistics.CalculateReadsStatistics(fn)
            statistics.do_stats_update(rdict, self.pid,
                                        self.special_run, self.local_dir
                                        )
            if paired:
                os.remove(stats_temp_file)

        if(self.assembler in settings.QUANTIFICATION_ASSEMBLERS):
            session = netutils.make_db_session()
            #add the quantification fn to upload
            self.output_files.add(self.local_dir + self.pn + \
                                   "_quantification.txt")
            jid = netutils.get_jid(self.pid, 'quantification', session)
            if(not jid is None):
                jn = netutils.get_table_object('jn_mapping', session)
                u = jn.update().where(
                                      jn.c.job_id==jid
                                      ).values(started='Y',
                                               finished='Y',
                                               f_ts=functions.\
                                               current_timestamp()
                                               )
                session.conn.execute(u)
        Job.complete(self)


class BlastAssociationJob(LocalJob):
    '''
        This class handles analysis tasks that associate annotation to a 
        transcript by using that transcript's annotation against a public 
        database, such as SwissProt or NR. 
        
        This class currently supports the following type of associations:
        
            BLAST to Gene Ontology Terms
            BLAST to KEGG Pathways
                  
        Constructor:
            KEGGProcess - if KEGG annotation
            GOProcess - if GO annotation
            
            arguments
                pid (int):
                    ID of the project this job belongs to
                jid (int):
                    Unique job identifer assigned to this job
                annot_type (str):
                    The type of the BLAST association annotation 
                    
                    Supports either 'kegg' or 'go' for this argument
                resources (str):
                    Comma separated list of resource:value pairs, in the form:
                    'cpu:1, mem:500mb'
                pipeline_arg_str:
                    String representation of the pipeline arguments from args
                    table
                process_arg_str:
                    String representation of the process arguments from the
                    args table
                
        input_flags:
            input: the BLAST anotation file for these transcripts against 
                SwissProt
                 
        output:
            A tab-delimited file of annotations of the input transcripts
            against GO or SwissProt. 
    '''

    def __init__(self, 
                 pid, 
                 jid, 
                 annot_type, 
                 pipeline_arg_str, 
                 process_arg_str
                 ):
        LocalJob.__init__(self, 
                          pid, 
                          jid, 
                          annot_type, 
                          resources="cpu:1", 
                          pipeline_arg_str=pipeline_arg_str, 
                         process_arg_str=process_arg_str
                          )
        if(annot_type == "kegg"):
            self.constructor = KEGGProcess
        elif(annot_type == "go"):
            self.constructor = GOProcess

    def complete(self):
        '''
            Calculates GO or KEGG statistics and stores them in the run_stats
            table. 
            
            Calls Job.complete()
        '''
        if(self.job_type == "go"):
            self.output_files.add(self.local_dir + self.pn + "_gocats.txt")
            fn = self.current_output
            rdict = statistics.CalculateGOStatistics(fn)
            statistics.do_stats_update(rdict, self.pid,
                                        self.special_run, self.local_dir)
        if(self.job_type == "kegg"):
            fn = self.current_output
            rdict = statistics.CalculateKEGGStatistics(fn)
            statistics.do_stats_update(rdict, self.pid,
                                       self.special_run, self.local_dir)
        Job.complete(self)

    def start(self):
        '''
            Creates either a KEGGProcess or GOProcess object, calls the 
            object's start method, and finally appends the object ot the list
             of processes for this job. 
        '''
        arg1 = self.pn
        out_f = self.current_output

        p = self.constructor(
                             arg1,
                             self.input_files['input'],
                             out_f, self.process_args,
                             self.special_run,
                             self.job_name
                             )
        p.start()
        self.processes.append(p)


class PantherJob(LocalJob):

    def __init__(
                 self,
                 pid,
                 jid,
                 jtype,
                 executable,
                 resources,
                 pipeline_arg_str,
                 process_arg_str
                 ):
        LocalJob.__init__(
                          self,
                          pid,
                          jid,
                          jtype,
                          executable,
                          resources,
                          pipeline_arg_str,
                          process_arg_str
                          )
        self.constructor = PantherProcess

    def start(self):
        print("Starting PantherJob " + self.job_name)
        if(self._input_exists()):
            ofile = self.current_output
            p = self.constructor(self.pn, self.input_files['input'],
                                 ofile, self.process_args,)
            p.start()
            self.processes.append(p)
        else:
            print "panther input DOES NOT exist."

    def complete(self):
        self.cleanup()
        self.mark_complete()


class QualityTrimJob(LocalJob):
    ''' 
        This class is responsible for implementing the functionality to trim
        low-quality regions from next-generation sequence data.
        
        Default resources:
            cpu -> 1
        
        Constructor:
            QualityTrimProcess
        
        input_flags:
            input - a FASTQ file of either shotgun reads or the forward reads
                from paired end sequencing
            end2 - a FASTQ file of the reverse reads from paired-end sequencing
                (optional)
                
        output:
            FASTQ file(s) named after the original sequence files, in which 
            each sequence has had low-quality regions trimmed.

        The choice of quality-trimming tool is made with the --quality_trimmer
        flag in pipeline_args. 
        
        Supported quality-trimming tools:
            cutadapt
            fastx_trimmer
    '''

    def __init__(self,
                 pid,
                 jid,
                 jtype,
                 executable,
                 resources,
                 pipeline_args,
                 process_args):
        LocalJob.__init__(self,
                          pid,
                          jid,
                          jtype,
                          executable,
                          resources,
                          pipeline_args,
                          process_args)

        if(self.paired_end):
            self.output_files.add(self.input_files['end2'])

        self.constructor = QualityTrimProcess


    def complete(self):
        self.cleanup()
        self.mark_complete()

    def start(self):
        p = self.constructor(self.pn, self.input_files, self.current_output,
                             self.pipeline_args, self.process_args,
                             self.resources, paired_end=self.paired_end)
        p.start()
        self.processes.append(p)


class ReadNormJob(LocalJob):
    '''
        This class represents the analysis step of performing read 
        normalization on short reads before assembly. 
        
        Please note that this step is memory-intensive, as the job invokes
        a script that must keep a kmer dictionary of all reads in memory during
        execution. On HiSeq data, expect to use ~20GB of memory per instance
        of this job.
        
        Default resourecs:
            cpu -> 1
            
        Constructor:
            ReadNormProcess
            
        input_flags:
            input - a FASTQ file of either shotgun sequencing reads or the
                forward reads from paired-end sequencing
            end2 - a FASTQ file of reverse end reads from paired-end sequencing
                (optional)
                
        output:
            A FASTA representation of each input sequence file after applying
            read normalization.
    '''
    def __init__(self, 
                 pid,
                 jid, 
                 jtype, 
                 executable, 
                 resources, 
                 pipeline_args, 
                 process_args):

        LocalJob.__init__(self,
                          pid,
                          jid,
                          jtype,
                          executable,
                          resources,
                          pipeline_args,
                          process_args)
        self.constructor = ReadNormProcess

    def start(self):
        p = self.constructor(self.pn,
                             self.input_files['input'],
                             self.current_output,
                             self.pipeline_args,
                             self.process_args,
                             self.resources,
                             paired_end=self.paired_end)
        p.start()
        self.processes.append(p)


class UploadJob(LocalJob):

    def __init__(self,
                 pid,
                 jid,
                 t):
        LocalJob.__init__(self,
                          pid,
                          jid,
                          t)

    def complete(self):
        self.mark_complete()

    def start(self):
        pass


class PipeProcess(Process):
    '''
        This is the base class for all analysis processes running within
        the Autonomics system. Each job class instantiates and starts one or
        more of these objects to perform analyses. Each PipeProcess runs in a 
        separate system process. 
        
        PipeProces objects starting remote analyses spend most of their 
        time sleeping, while the remote process executes, and as such
        are not considered by the system to utlize any CPU time. 
        
        PipeProcess objects that perform local analyses do need to request
        resources from the system and release them upon completion.
    '''

    def __init__(self):
        Process.__init__(self)
        self.location = Locations.LOCAL
        self.input_files = {}


    def _exec_cmd(self, cmd, blocking=True, capture_out=False):
        '''
            cmd (str):
                The fully-formed command you wish to execute. This command
                is sent to the underlying system via the subproces.Popen
                constructor. 
                
            blocking (Boolean):
                Whether or not this method should block, waiting for the 
                command executed to return. 
                
                This ethod will hang if True, returns immediately if False.
                
            capture_out (Boolean):
                Whether or not input/output/error from the child process
                is captured.
                
                If True, stdin, stderr, and stdout are redirected to 
                subprocess.PIPE
                
            Executes the command 'cmd' by instantiating a Popen object. 
            Optionally blocks the calling thread until execution finishes. 
            Optionally captures stdin/stdout/stderr of the cmd.
            
            Returns the Popen object created within the method.
        '''
        redirect = None
        if(capture_out):
            redirect = subprocess.PIPE
        p = subprocess.Popen(cmd, shell=True, stdin=redirect, stdout=redirect)

        if(blocking):
            p.wait()

        return p

    def _get_cmd_template(self, prog, job_type, session):
        '''
            prog (str):
                The name of the program for which a command template is
                returned
                
            job_type (str):
                The job_type that uitlizes the specified program
                
            session (netutils.DBSession):
                A session object with an active connection to the Autonomics
                database
                
            Uses the program name and job_type to look up a command template
            in the executables table of the Autonomics database. This 
            command template has the form:
            
                executable <argument placeholder>* <hard-coded arguments>* 
                    <redirect>
                    
            Returns a command template if an entry is found in the executables
            table, None otherwise.
        '''
        location_str = ''
        if(self.location == Locations.LOCAL):
            location_str = 'local'
        elif(self.location == Locations.HPC):
            location_str = 'HPC'
        else:
            raise Exception("Unknown location value: " + str(self.location))

        #retrive how this command is run
        execs = netutils.get_table_object("executables", session)
        command_res = execs.select(and_(execs.c.exec_name==prog,
                                        execs.c.location==location_str,
                                        execs.c.job_type==job_type)).execute()
        template = command_res.fetchone()
        if(template is None):
            #no executable found, can't automatically construct command string
            return None

        return template.command + " " + template.args

    def _prepare_cmd(self, cmd_template, replacement_dict):
        '''
            Prepares a fully-formed command  by replacing argument
            placeholders in a command template with the values stored in a 
            replacement dictionary. 
            
            The method uses the keys of the dictionary to search for 
            placeholders within the command template. When a match is found,
            the method replaces the matched placeholder with the value from the
            dictionary.
        '''

        ret_cmd = cmd_template
        for key, value in replacement_dict.items():
            ret_cmd = ret_cmd.replace("<" + key + ">", str(value))
        return ret_cmd

    def _proc_options_dict(self, job_type, session):
        '''
            Looks up the default process_options for the given job_type in 
            the proc_options table and returns a dictionary of 
            placeholder->value pairs for replacement into a command template.
        '''
        #get the list of possible process options
        proc_options = netutils.get_table_object("proc_options", session)
        proc_opt_res = proc_options.select(
                                           proc_options.c.job_type==job_type
                                           ).execute()

        ret = {}
        #set user-supplied flags
        for row in proc_opt_res.fetchall():
            val = row.default_value
            if(hasattr(self.process_args, row.flag)):
                val = getattr(self.process_args, row.flag)
            ret[row.flag] = val

        return ret

    def run(self):
        pass
        #this method needs to do the following: check if the job_type supports num_threads or job_types.embarassingly_parallel = 'Y'
        #if it does not support num_threads, and is embarassingly parallel, split the input and start multiple commands based on the cpu resources for the process
        #if either of the above is true, simply use the process's _prepare_cmd method and execute the returned command


class AdapterTrimProcess(PipeProcess):
    '''
        This process performs the actual work of removing adapters from 
        sequences prior to downstream analysis.
        
        This process currently supports adapter trimming with either the 
        cutadapt or fastx_trimmer softwarepackges. 
        
        input:
            One or two FASTQ files to have adapter sequences removed
            
        output:
            One or two FASTQ files (named identically to the input), that
            have all adapter subsequences removed
    '''
    
    def __init__(
                 self,
                 project_name,
                 input_files,
                 output_file,
                 process_args,
                 pipeline_args,
                 resources,
                 paired_end=False
                 ):
        PipeProcess.__init__(self)
        self.project_name = project_name
        self.input_files = input_files
        self.current_output = output_file
        self.process_args = process_args
        self.pipeline_args = pipeline_args
        self.resources = resources
        self.paired_end = paired_end

    def trim_adapters(self, f, args):
        '''
            f (str):
                Path to the file from which adapters will be removed
            
            args (Arguments):
                An Arguments object containing the pipeline args to control
                this adapter trimming.
                
            THis method gets a list of adapters to remove from the sequence
            file, f, using the known_adapters and project_adapters tables.
            
            These adapters are then removed from the sequences in the file f 
            using the program specified by the adapter_trimmer attribute
            of args.
            
            When the run method of this class exits, a file with the same
            name of the input file should be created, but the sequences in that
            file should have all adapter subsequences removed.
        '''

        tmp = f + ".trimming.tmp"
        prog = "fastx_clipper"
        if(hasattr(args, "adapter_trimmer")):
            prog = args.adapter_trimmer

        #get the adaptors sequences
        session = netutils.make_db_session()
        adapts = netutils.get_adapter_rows(netutils.get_pid(self.project_name,
                                                        session), session)
        if len(adapts) ==0: print "There are no known_adapters for this project, so skipping adapter_trim"
        for adapt in adapts:
            adapt_flg = adapt.adapter_sequence
            if(prog == 'cutadapt'):
                if(adapt.end == 3):
                    adapt_flg = '-a ' + adapt_flg
                else:
                    adapt_flg = '-g '+ adapt_flg
            self.process_args.adapter = adapt_flg
            cmd = self._get_cmd_template(prog, "adapter_trim", session)
            param_dict = self._proc_options_dict("adapter_trim", session)
            param_dict['input'] = f
            param_dict['output'] = tmp
            fc = f + '.before.adapter_trim'
            shutil.copyfile(f,fc )
            cmd = self._prepare_cmd(cmd, param_dict)
            p = self._exec_cmd(cmd)
            die_on_error(p.returncode, cmd_str=cmd)
            #overwrite the original file with the trimmed file
            shutil.move(tmp, f)
        session = None

    def run(self):
        '''
            Iterates over input files, removing all adapters assigned to this
            project in project_adapters from those files.
            
            Produces files with the same names as the AdapterTrimJob's input
            files, without adapter sequences in them.
        '''
        #parse the process args
        files = [self.input_files['input']]
        if(self.paired_end):
            #add the second input file
            files.append(self.input_files['end2'])

        for f in files:
            self.trim_adapters(f, self.process_args)


class AssemblyProcess(PipeProcess):
    '''
        This class represents the task of locally assembling short reads
        into either transcripts or genomic regions.
        
        Default resources:
            cpu -> 20
        
        input_flags:
            input - Either a FASTQ file of shotgun reads or the forwrad end
                reads from paired-end sequencing
            
            end2 - A FASTQ file of the reverse reads from paired-end sequencing
                (optional)
                
        output:
            A FASTA file of assembled sequences, and a quantification file, if
            the assembly is performed on transcriptomic sequence and the 
            assembler supports automatic estimation of transcript abundance.
            
        The underlying assembler can be changed by specifying the --assembler
        flag in the pipeline_args field in the args table for this job. 
        
        Supported values for --assembler are: mira, trinity
    '''

    def __init__(
                 self,
                 project_name,
                 input_file,
                 output_file,
                 process_args,
                 paired_end,
                 assembler,
                 cpus,
                 special_run,
                 job_name):

        PipeProcess.__init__(self)
        self.pn = project_name
        self.input_file = input_file
        self.current_output = output_file
        self.process_args = process_args
        self.paired_end = paired_end
        self.assembler = assembler
        self.cpus = cpus
        self.special_run= special_run
        self.job_name = job_name

    def run(self):

        out_dir = settings.proj_dir + self.pn + "/"
        if (self.paired_end):
            if (self.assembler == "trinity"):
                in_file1 = ""
                in_file2 = ""
            if (self.special_run):
                in_file1 = self.input_file
            else:
                in_file1 = out_dir + self.pn + ".fasta"
                in_file2 = out_dir + self.pn + ".end2.fasta"
            cmd = "python " + settings.SCRIPTPATH + "run_trinity.py \
                    -in_file1 " + in_file1 + " -in_file2 " + in_file2 + \
                    " -cpus " + str(self.cpus) 
        else:
            prefix = ""
            if (self.special_run):
                out_dir = settings.special_runs_dir + self.job_name + '/'
                prefix = self.job_name
            else:
                prefix = self.pn
                self.input_file = out_dir + self.pn + ".fasta"
            if (self.assembler == "trinity"):
                cmd = "python " + settings.SCRIPTPATH + "run_trinity.py \
                    -in_file1 " + self.input_file + " -cpus " + str(self.cpus)
            elif (self.assembler == "mira"):
                cmd = "python " + settings.SCRIPTPATH + "assembly_pipeline.py \
                    -mira -o " + prefix + " -cpu " + str(self.cpus) + " fastq \
                    " + out_dir
            elif (self.assembler == "miranewbler"):
                cmd = "python " + settings.SCRIPTPATH + "assembly_pipeline.py\
                     -o " + prefix + "-cpu " + str(self.cpus) + " fastq "\
                      + out_dir
        p = self._exec_cmd(cmd)
        die_on_error(p.returncode, cmd_str=cmd)


class BlastAssociationProcess(PipeProcess):
    '''
        A PipeProcess object for performing BLAST-association annotations. 
        In these analyses, annotations from known, public proteins are 
        associated with input sequences by utilizing similarity captured in 
        BLAST ouput and mapping the annotations from the known proteins to
        your input sequences.
        
        All analyses tasks utilizing Blast association will subclass this class.
        
        Currently subclassing objects: GOProcess, KEGGProcess
    '''

    def __init__(
                 self,
                 base_name,
                 blast_file,
                 out_file,
                 process_args,
                 num_associations=10
                 ):
        PipeProcess.__init__(self)
        self.base_name = base_name
        self.local_dir = settings.proj_dir + "/" + base_name + "/"
        self.blast_file = blast_file
        self.current_output = out_file
        self.process_args = process_args
        self.key_prefix = None

    def run(self):
        '''
            This method maps annotation from a database of sequences to a set
            of query sequences by utilizing similarity information presented 
            as a textual alignment file. Currently, this method only supports
            mapping annotations from SwissProt using BLAST-style annotation 
            files.
            
            To perform the mapping, the method inspects each alignment between
            a query sequence and the SwissProt database sequence. If the two
            sequences share enough similarity, the annotation from the 
            SwissProt sequence is mapped to the query sequence, by looking up
            the SwissProt accession number for the database sequence in a redis
            keystore. This keystore maps the accession number plus an 
            annotation source identifier to the annotations for that sequence
            in that annotation source.
            
            For example, the key zero_click:goa:AA1451 stores all of the Gene
            Ontology terms mapped to the SwissProt sequence with accession
            AA1451.
        '''
        rServer = redis.Redis(host = settings.REDIS_HOST, 
                              port = settings.REDIS_PORT, db=0)
        r = AlignmentReader(self.blast_file, "hpc-blast")
        r.read()
      
        out = open(self.current_output, 'w')
        for query, hits in r.annotations.items():
            for hit in hits:
                if(hit.significance > float(self.process_args.evalue)):
                    continue
                reference = hit.reference_id.split("|")[1]
                reference = reference.split(".")[0]
                #for each hit asccesion, get all of the associated records in the db
                records = rServer.lrange(self.key_prefix + reference, 0, -1)
                for record in records:
                    out.write(query + "\t" + record + "\t" + str(hit.score) + \
                               "\t" + str(hit.significance) + "\n")
        out.close()


class GOProcess(BlastAssociationProcess):
    '''
        Class that represents performing BLAST-association annotation of Gene
        Ontology terms to input sequences. 
        
        Additional attributes:
            - self.key_prefix = "zero_click:goa:"
                Key prefix that allows for selection of Gene Ontology terms
                for a given SwissProt sequences from the redis keystore. The 
                full key is formed by appending the SwissProt accession to the
                key_prefix
        
        input:
            A BLAST-style alignment output between a set of query sequences
            and the SwissProt database
            
        output:
            A tab-delimited file containing Gene Ontology anotations for the
            query sequences contained in the input BLAST file
    '''
    def __init__(
                 self,
                 base_name,
                 blast_file,
                 out_path,
                 process_args,
                 special_run,
                 job_name
                 ):
        BlastAssociationProcess.__init__(
                                         self,
                                         base_name,
                                         blast_file,
                                         out_path,
                                         process_args
                                         )
        self.key_prefix = "zero_click:goa:"
        self.special_run = special_run
        self.job_name = job_name

    def run(self):
        '''
            Performs basic mapping of Gene Ontology terms from known sequences
            in SwissProt to a set of query sequences, based on alignment 
            similarity.
            
            After creating the initial Gene Ontology annotaion file, limits the
            number of GO annotations if the total annotations in the file
            exceeds a given threshold, currently set at 200,000. This
            limiting is done using the flatten functionality present in the
            filetools.py script. The script will limit the maximum number of 
            annotations assigned to a given input sequence to 10. 
            
            This method also creates a GO category file containing 1) the 
            number of unique sequences placed into each GO category 2) the
            number of unique annotations present in each category and 3) the
            abundance of all sequences placed into that category, if 
            quantification data for the input sequences is detected. 
        '''
        BlastAssociationProcess.run(self)
        #run two rounds of file-flattening
        #first, flatten based on seq_id, go_term
        #then flatten based on seq_id, limit 10 records by default (maybe allow pipeline parameter to change this?)

        p = self._exec_cmd("python " + settings.SCRIPTPATH + "/filetools.py  \
                --fields seq_id:1 sp_acc:2 go_term:3 go_ref:4 comaprtment:5 \
                gene_id:6 score:7 evalue:8:float --key-col 1 3 --flatten \
                --filter evalue:lt:1e-04 --flatten-depth 1 " + \
                self.current_output)
        die_on_error(p.returncode)
        p = self._exec_cmd("wc -l " + self.current_output + "_flattened.txt",
                           capture_out=True)

        die_on_error(p.returncode)
        output = p.stdout.readline()
        num = convert_if_int(output.split(" ")[0])
        tmp_outs = self.current_output + "_flattened.txt"
        if(num > 200000):
            #having more than 200000 records triggers flattening

            p = self._exec_cmd("python " + settings.SCRIPTPATH + "/filetools.py --fields seq_id:1 sp_acc:2 go_term:3 evalue:8:float --key-col 1 --flatten --flatten-depth 10 " + self.current_output  + "_flattened.txt")
            die_on_error(p.returncode)
            os.remove(tmp_outs)
            tmp_outs = self.current_output + "_flattened.txt_flattened.txt"

        p = self._exec_cmd("mv " + tmp_outs + " " + self.current_output)
        die_on_error(p.returncode)

        #run the GOAnnotator software on the annotations to get category-level summaries
        #check if the quantification file exists for this project (it should, if not categories will lack abundance data)

        quant_file = ""
        mid_arg = ""
        if self.special_run:
            quant_file = self.current_output + "_quantification.txt"
            mid_arg = self.job_name
        else:
            quant_file = self.local_dir + self.base_name + "_quantification.txt"
            mid_arg = self.base_name
        java_cmd = "java -Xmx10g goannot8r.GOAnnotator " + mid_arg + " "  + \
            self.current_output

        if(os.path.exists(quant_file)):
            java_cmd += " " + quant_file

        if self.special_run:
            mid_arg = self.current_output
        else:
            mid_arg = self.local_dir + self.base_name

        java_cmd += " > " + mid_arg + "_gocats.txt"

        proc = subprocess.Popen(java_cmd, shell=True)
        proc.wait()
        die_on_error(proc.returncode, cmd_str=java_cmd)


class KEGGProcess(BlastAssociationProcess):
    '''
        Class that represents performing BLAST-association annotation of KEGG
        pathway terms to input sequences. 
        
        Additional attributes:
            - self.key_prefix = "zero_click:koa:"
                Key prefix that allows for selection of Gene Ontology terms
                for a given SwissProt sequences from the redis keystore. The 
                full key is formed by appending the SwissProt accession to the
                key_prefix
        
        input:
            A BLAST-style alignment output between a set of query sequences
            and the SwissProt database
            
        output:
            A tab-delimited file containing KEGG anotations for the
            query sequences contained in the input BLAST file
    '''
    def __init__(
                 self,
                 base_name,
                 blast_file,
                 out_file,
                 process_args,
                 special_run,
                 job_name
                 ):
        BlastAssociationProcess.__init__(
                                         self,
                                         base_name,
                                         blast_file,
                                         out_file,
                                         process_args)
        self.key_prefix = "zero_click:koa:"


class PantherProcess(PipeProcess):

    def __init__(
                 self,
                 project_name,
                 input_file,
                 output_file,
                 process_args
                 ):

        PipeProcess.__init__(self)
        self.project_name = project_name
        self.input_file = input_file
        self.current_output = output_file
        self.process_args = process_args

    def run(self):
        self.input_file = translate_seq_file(self.input_file, "fasta")
        cmd = "perl " + settings.INSTALL_DIR + "scripts/pantherScore.pl -V -l \
                " + settings.panther_data_path + " -D A -V -n -i " + \
                self.input_file + " -o " + self.current_output
        process_obj = subprocess.Popen(cmd, shell=True)
        process_obj.wait()


class PfamProcess(PipeProcess):
    '''
        This class presents the functionality to annotate transcripts with 
        Pfam domains.
        
        input:
            A set of transcript or protein sequences
            
        output:
            A tab-delimited file of Pfam annotations for the given input
            sequences
    '''

    def __init__(self, seq_file, output_file):
        Process.__init__(self)
        self.input_file = seq_file
        self.current_output = output_file

    def run(self):
        '''
            This method annotates a set of input sequences with Pfam domains.
            
            If the input sequences are not in protein space, the method first
            translates them to protein in all six reading frames. These 
            translated sequences are placed in a file named 
            <original_file>_translated.fa.
            
            pfam_scan.pl is then run on the translated sequences, producing
            the Pfam annotation file.
        '''
        #translate the input file
        self.translated = translate_seq_file(self.input_file, "fasta")
        os.chdir(settings.pfam_exec_path)
        #run pfam
        p = subprocess.Popen("perl pfam_scan.pl -cpu 1 -outfile " + \
                             self.current_output + " -dir " + \
                             settings.pfam_data_path + " -fasta " + \
                             self.translated, shell=True)
        #wait for pfam to finish
        p.wait()

class QualityTrimProcess(PipeProcess):
    '''
        This PipeProcess performs the  work of removing low quality regions from 
        sequences prior to downstream analysis.
        
        This process currently supports quality trimming with either cutadapt 
        software package. 
        
        input:
            One or two FASTQ file
            
        output:
            One or two FASTQ files (named identically to the input), that
            have low-quality regions removed
    '''
    def __init__(
                 self,
                 project_name,
                 input_files,
                 output_file,
                 pipeline_args,
                 process_args,
                 resources,
                 paired_end=False
                 ):
        PipeProcess.__init__(self)
        self.project_name = project_name
        self.input_files = input_files
        self.current_output = output_file
        self.pipeline_args = pipeline_args
        self.process_args = process_args
        self.resources = resources
        self.paired_end = paired_end

    def trim_quality(self, f, args):
        '''
            f (str):
                Path to a sequence file
            args (Arguments):
                An Arguments object representing the pipeline_args field 
                initially provided to the Job class creating this PipeProcess
                
            This method determines which quality trimmer is going to be used, 
            and then generates a fully-formed command from the template 
            command for that quality trimmer and the supplied process_args from
            the args table. 
            
            The command is then executed locally, trimming low-quality regions
            from the input sequences.
        '''
        session = netutils.make_db_session()
        tmp = f + ".quality.tmp"
        prog = "cutadapt"
        if(hasattr(args, "quality_trimmer")):
            prog = args.quality_trimmer

        cmd = self._get_cmd_template(prog, 'quality_trim', session)
        param_dict = self._proc_options_dict('quality_trim', session)
        param_dict['input'] = f
        param_dict['output']= tmp
        fc = f + '.before.quality_trim'
        shutil.copyfile(f,fc )
        cmd = self._prepare_cmd(cmd, param_dict)
        #use the trimmer to trim quality

        p = self._exec_cmd(cmd)
        die_on_error(p.returncode)
        #re-write the original file as the trimmed file
        shutil.move(tmp, f)

    def run(self):
        #parse the process args
        files = [self.input_files['input']]
        if(self.paired_end):
            #add the second input file
            files.append(self.input_files['end2'])
        for f in files:
            self.trim_quality(f, self.process_args)


class ReadNormProcess(PipeProcess, Job):
    '''
        This PipeProcess provides the functionality to perform read 
        normalization on input sequences. 
        
        Read normalization is currently performed using the khmer package, 
        available at: https://github.com/ctb/khmer
        
        input:
            One or more files in FASTQ format. If two files are provided, 
            the 'input' file is assumed to be the file of forward reads from 
            paired-end sequencing, while the 'end2' file is assumed to be the
            file of reverse reads.
            
        output:
            One or more files in FASTA format, where the sequences in each file
            are those from the original files that passed read normalization.
    '''
    def __init__(
                 self,
                 project_name,
                 input_file,
                 output_file,
                 pipeline_args,
                 process_args,
                 resources,
                 paired_end=False
                 ):
        PipeProcess.__init__(self)
        self.project_name = project_name
        self.input_file = input_file
        self.current_output = output_file
        self.pipeline_args = pipeline_args
        self.process_args = process_args
        self.resources = resources
        self.paired_end = paired_end

    def run_khmer(self):
        '''
            This method runs khmer in either single-ended or paired-ened mode,
            depending on the input.
        '''
        if(self.paired_end):
            cmd = "python " + settings.SCRIPTPATH + "run_read_normalization.py -in_file\
                 " + self.input_file + " -paired_end 1"
        else:
            cmd = "python " + settings.SCRIPTPATH + "run_read_normalization.py -in_file\
                 " + self.input_file
        p = self._exec_cmd(cmd)
        die_on_error(p.returncode, cmd)

    def run(self):
        #parse the process args
        self.pipeline_args.parse()
        self.process_args.parse()
        self.run_khmer()


class BlastJob(HPCJob):
    '''
        Base class that represents a BLAST-based annotation job on a remote
        high-performance computing cluster.
        
        This class is designed to utilize the BLAST+ family of alignment tools
        available from NCBI.
        
        All BLAST jobs require that you specify the number of cpu's to use
        in the 'resource' field of the args table for the job.
        
        Default resources:
            mem -> 500mb
            wall -> 24:00:00
            modules -> module load ncbi_blast
            ppn -> 4
            
        input_flags:
            input - a FASTA formatted sequence file containing the query
                sequences
                
        output:
            A BLAST annotation file resulting from the alignment of the
            input sequences against the specified database    
    '''

    def __init__(
                 self,
                 pid,
                 jid,
                 job_type,
                 executable,
                 resources,
                 pipeline_args,
                 process_args
                 ):
        HPCJob.__init__(
                        self,
                        pid,
                        jid,
                        job_type,
                        executable,
                        resources,
                        pipeline_args,
                        process_args
                        )
        self.process_args.arg_string += " -num_threads " + \
                                        str(self.resources['ppn'])
        self.hpc_qsub_files['input'] = self.input_files['input']

    def _default_resources(self):
        ret = HPCJob._default_resources(self)
        ret['mem'] = "500mb"
        ret['wall'] = "24:00:00"
        ret['modules'] = "module load ncbi_blast"
        ret['ppn'] = 4
        return ret

    def _handle_args(self):
        '''
            Determines which alignment program to use by checking the --aligner
            flag of the provided pipeline_args. Sets self.executable equal to 
            the found value, if any.
        '''
        HPCJob._handle_args(self)
        if(hasattr(self.pipeline_args, 'aligner')):
            self.executable = self.pipeline_args.aligner


class BlastJobNR(BlastJob):
    '''
        An analysis job class for the alignment of query sequences to the 
        non-redundant (NR) protein database from NCBI.
        
        Default resources:
            mem -> 5000mb
            wall -> 192:00:00
            blast_nr -> 1
            
            Manager is expected to set the number of blast nr jobs that can 
            run in parallel by setting a resource level for blast_nr. 
            These jobs should be limited as otherwise they will request 
            dedicated cores on the remote HPC cluster and use them for days,
            preventing other, smaller jobs from running. 
            
        This class automatically appends --db nr to it's process_args, 
        which is then replaced into the command string at the <db> placeholder.
    '''

    def __init__(
                 self,
                 pid,
                 jid,
                 job_type,
                 executable,
                 resources,
                 pipeline_args,
                 process_args
                 ):
        BlastJob.__init__(
                          self,
                          pid,
                          jid,
                          job_type,
                          executable,
                          resources,
                          pipeline_args,
                          process_args
                          )
        self.process_args.arg_string += " -db nr"

    def _default_resources(self):
        ret = BlastJob._default_resources(self)
#        ret['mem'] = "12000mb"
        ret['mem'] = settings.BLAST_NR_MAX_MEM
#        ret['wall'] = "192:00:00""99:00:00"
        ret['wall'] = settings.BLAST_NR_MAX_WALL_TIME
        ret['blast_nr'] = 1
        return ret

    def complete(self):
        f = self.current_output
        rdict = statistics.CalculateBlastNRStatistics(f)
        statistics.do_stats_update(rdict, self.pid,
                                   self.special_run, self.local_dir)
        BlastJob.complete(self)


class BlastJobSwissprot(BlastJob):
    '''
        An analysis job class for the alignment of query sequences to the 
        SwissProt protein database.
        
        Default resources:
            mem -> 1000mb
            wall -> 48:00:00
            
        This class automatically appends --db swissprot to it's process_args, 
        which is then replaced into the command string at the <db> placeholder.
    '''

    def __init__(self, pid, jid, jobType, executable, resources,
                 pipeline_args, process_args):
        BlastJob.__init__(self, pid, jid, jobType, executable, resources,
                          pipeline_args, process_args)
        self.process_args.arg_string += " -db swissprot"

    def _default_resources(self):
        ret = BlastJob._default_resources(self)
#        ret['mem'] = "2000mb"
        ret['mem'] = settings.BLAST_SWISSPROT_MAX_MEM
#        ret['wall'] = "48:00:00"
        ret['wall'] = settings.BLAST_SWISSPROT_MAX_WALL_TIME
        return ret

    def complete(self):
        f = self.current_output
        rdict = statistics.CalculateBlastSwissprotStatistics(f)
        statistics.do_stats_update(rdict, self.pid,
                                   self.special_run, self.local_dir)
        BlastJob.complete(self)


class BlatJob(HPCJob):
    '''
        Not fully implemented
    '''
    def __init__(
                 self,
                 pid,
                 jid,
                 job_type,
                 executable,
                 resources,
                 pipeline_args,
                 process_args
                 ):

        HPCJob.__init__(
                        self,
                        pid,
                        jid,
                        job_type,
                        executable,
                        resources,
                        pipeline_args,
                        process_args
                        )

    def _default_resources(self):
        ret = HPCJob._default_resources(self)
        ret['wall'] = "86:00:00"
        ret['modules'] = "module load blat"
        return ret


class PfamJob(HPCJob, threading.Thread):
    '''
        This class represents a remote Pfam annotation job.
        
        Default resources:
            wall -> 48:00:00
            module -> module load hmmer
            
        input_flags:
            input - a FASTA file of sequences to be supplied to pfam_scan.pl
                for annotation
                
        output:
            A tab-delimited file of Pfam annotations
            
        This class behaves much the same as the local PfamJob class. The main
        difference is that this class leverages query splitting to achieve
        simplistic concurrency on a remote computational resource.
        
        The main query file is split and each split file is translated to 
        protein space, if necessary. Files are then transferred to the remote
        machine along with a qsub script, required to start the job.
        
        The process of splitting input files, creating qsubs, transferring data,
        and monitoring status is described within the HPCJob/HPCProcess class
        documentation. 
    '''
    def __init__(
                 self,
                 pid,
                 jid,
                 job_type,
                 executable,
                 resources,
                 pipeline_args,
                 process_args
                 ):
        threading.Thread.__init__(self)
        HPCJob.__init__(
                        self,
                        pid,
                        jid,
                        job_type,
                        executable,
                        resources,
                        pipeline_args,
                        process_args
                        )
        self.executable = "pfam_scan.pl"
        self.process_args.arg_string = " -dir " + settings.pfam_data_path + \
            " -cpu " + str(self.resources["ppn"]) + " " + process_args
        self.pipeline_args.parse()
        self.hpc_qsub_files['input'] = self.input_files['input']

    def _default_resources(self):
        ret = HPCJob._default_resources(self)
        ret['wall'] = settings.PFAM_MAX_WALL_TIME
#        ret['wall'] = "48:00:00"
        ret['modules'] = "module load hmmer"
        return ret

    def start(self):
        threading.Thread.start(self)
        self.join()

    def run(self):
        #check if we need to translate the input
        if(hasattr(self.pipeline_args, "translate")):
            self.input_files['input'] = translate_seq_file(self.input_files[\
                                            'input'], "fasta")
            self.hpc_qsub_files['input'] = self.input_files['input']

        self.executable = "pfam_scan.pl"
        HPCJob.start(self)

    def complete(self):
        f = self.current_output
        rdict = statistics.CalculatePfamStatistics(f)
        statistics.do_stats_update(rdict, self.pid,
                                   self.special_run, self.local_dir)
        fa_file = self.local_dir + "*.fa"
        for filename in glob.glob(fa_file):
            os.remove(filename)
        HPCJob.complete(self)


class QuantificationJob(HPCJob):
    '''
        This class provides functionality to assign expression level estimates
        to a set of reference sequences. Specifically, this class is designed
        to determine expression level estimates on a remote computational 
        resource.
        
        Default resources:
            ppn -> 4
            wall -> 24:00:00
            
        Constructor:
            HPCProcess
            
        Attributes:
            query_files ({input_flag:file_path}):
                A dictionary of input_flags to file paths, where each file
                present as a value in the dictionary is used as a query file
                in the alignment step of the quantification job
                
            local_db_path (str):
                Local path to the database file (sequences to be quantified)
                
            remote_db_path (str):
                Path to the database file on the remote resource
                
            db_files ([str]):
                A list of all database files. All files in this list are
                renamed relative to the remote working directory and moved
                to the remote resource
                
            aligner_modules ({str: str}):
                 A dictionary of module strings. The class decides which 
                 modules it needs to load on the remote computational resource
                 by looking into this dictionary, using the quantification type
                 as the key
                 
            prepare_for_aln ({str: method}):
                A dictionary of methods used to prepare the files used in the 
                quantification alignment. The class determines which method
                to use by looking into this dictionary with the quantification
                type.
                
            hpc_qsub_files ({str: str}):
                A dictionary of string to string mappings, representing the
                input files to this analysis that need to be split to achieve
                concurrency. Each file stored in the dictionary is split 
                and moved to the remote resource, after renaming each of the 
                split files to it's correct remote path. 
                         
        input_flags:
            query - a FASTA/Q file of shotgun reads or the forward end reads 
                from paired-end sequencing
            query2 - A FASTA/Q file of the reverse reads from paired-end 
                sequencing (optional)
            db - a FASTA file of reference sequences
            
        output:
            A tab-delimited file of transcript expression estimates
    '''
    
    class QTypes:
        '''
            Class that emulates an 'enum' type and holds attributes for the 
            supported types of quantification jobs.
            
            Attributes:
                BLAST (int)=0:
                    Attribute representing the BLAST quantification type
                BOWTIE (int)=1:
                    Attribute representing the bowtie quantification type
        '''

        BLAST = 0
        BOWTIE = 1

        types = [BLAST, BOWTIE]

    def __init__(
                 self,
                 pid,
                 jid,
                 t,
                 executable = "NA",
                 resources = {},
                 pipeline_arg_str = "",
                 process_arg_str = ""
                 ):
        self.query_files = {}
        self.remote_db_path = ""
        self.db_files = []
        self.paired_end = False
        HPCJob.__init__(
                        self,
                        pid,
                        jid,
                        t,
                        executable=executable,
                        resources=resources,
                        pipeline_arg_str=pipeline_arg_str,
                        process_arg_str=process_arg_str
                        )
        self.complete_quant_type = {
                                    QuantificationJob.QTypes.BLAST: \
                                    self.complete_blast,
                                    QuantificationJob.QTypes.BOWTIE: \
                                    self.complete_bowtie
                                    }
        self.aligner_modules = {self.QTypes.BLAST: "module load ncbi_blast",
                                self.QTypes.BOWTIE: "module load bowtie"}
        self.prepare_for_aln = {self.QTypes.BLAST: self.prep_blast,
                                self.QTypes.BOWTIE: self.prep_bowtie}
        self.hpc_qsub_files['query'] = self.input_files['query']
        self.local_db_path = self.input_files['db']
        self.db_file = os.path.split(self.local_db_path)[1]
        self.remote_db_path = self.remote_dir + self.db_file
        self.process_args.db = self.remote_db_path

    def _cleanup_local_files(self):
        '''
            Removes temporary files created during the class's execution.
        '''
        for f in self.db_files:
            if(f != self.local_db_path and os.path.isfile(f)):
                os.remove(f)
        HPCJob._cleanup_local_files(self)

    def _default_resources(self):
        resources = HPCJob._default_resources(self)
        resources['ppn'] = 4
#        resources['wall'] = "24:00:00"
        resources['wall'] = settings.QUANTIFICATION_MAX_WALL_TIME
        return resources

    def _handle_args(self):
        '''
            Determines the type of the quantification task, depending on the
            value of the --aligner flag in pipeline_args.
        '''
        Job._handle_args(self)

        self.db_type = attr_or_default(self.pipeline_args, "db_type", '')

        BLAST_TYPE = QuantificationJob.QTypes.BLAST
        BOWTIE_TYPE = QuantificationJob.QTypes.BOWTIE

        aligner = attr_or_default(self.pipeline_args, "aligner", '')

        self.quant_type = BOWTIE_TYPE

        if('blast' in aligner):
            self.quant_type = BLAST_TYPE

    def complete(self):
        '''
            Calls the appropriate 'complete' method stored in 
            complete_quant_type, depending on the type of this quantification
            job. 
            
            Calculates quantification statistics for the job, cleans up local
            files, uploads output, and marks the job as complete. 
        '''
        self.complete_quant_type[self.quant_type]()
        f = self.current_output
        rdict = statistics.CalculateQuantificationStatistics(f)
        statistics.do_stats_update(rdict, self.pid,
                                   self.special_run, self.local_dir)
#        self.upload()
        self.cleanup()
        self.mark_complete()

    def complete_blast(self):
        '''
            Performs complete actions related to a BLAST quantification type.
            
            1) Moves the original output file to a temporary file.
            2) Calls the seqtools.py script under the INSTALL_DIR + scripts
                directory to calculate quantification from the BLAST output
        '''
        tmp = self.local_dir + "tmp_quant.txt"
        shutil.move(self.current_output, tmp)
        cmd = "".join(["python ", settings.SCRIPTPATH,
                       "seqtools.py -at blast-tab -ft fasta -f ",
                        self.query_files[self.quant_type], " -a ",
                        tmp, " --quantification ",
                        "> ", self.current_output])
        p = self._exec_cmd(cmd)
        die_on_error(p.returncode, cmd)
        os.remove(tmp)

    def complete_bowtie(self):
        '''
            Moves the original bowtie output to a temporary file before 
            counting the number of mappings to each reference sequence and
            creating a new quantification file with the results.
        '''
        #move the quantification file
        tmp_quant = self.local_dir + "quant_tmp.txt"
        #change this to use shutil

        proc = subprocess.Popen("mv " + self.current_output + " " + \
                                tmp_quant, shell=True)
        proc.wait()
        #open the temporary quantification file
        results = {}
        fh = open(tmp_quant, 'r')
        for line in fh:
            line = line.rstrip("\n")
            if(line.startswith("@")):
                continue
            r = make_record(line, "sam")
            #no alignment if the flags are cleanly divisible by 4
            if(int(r.flags) & 4 != 0):
                continue
            if(r.reference in results):
                results[r.reference] += 1
            else:
                results[r.reference] = 1
        fh.close()
        out = open(self.current_output, 'w')
        for reference, records in results.items():
            out.write(reference + "\t" + str(records) + "\n")

        out.close()
        os.remove(tmp_quant)
        ebwt_file = self.local_dir + "*.ebwt"
        for filename in glob.glob(ebwt_file) :
            os.remove( filename )
        if(hasattr(self.pipeline_args, "paired_end")):
            os.remove(self.local_dir + self.pn + "_combined_for_quant.fastq")

    def move_db_files(self):
        '''
            Moves databse files contained in self.db_files to the remote 
            computational resource, after renaming them to their full remote
            path.
        '''
        c = netutils.ssh_connect(self.ssh_credentials)
        make_remote_dir(settings.remote_dir + self.run_name, c)
        for fp in self.db_files:
            if(os.path.exists(fp)):
                fn = os.path.split(fp)[1]
                c.put(fp, self.remote_dir + fn)
        c.close()

    def prep_blast(self):
        '''
            Prepares DB and query files for a quantification job using BLAST.
            
            1) adds self.local_db_path to db_files
            2) constructs a BLAST index for the file pointed to by 
                local_db_path
        '''
        self.db_files = [self.local_db_path]
        if(self.db_type == "NT"):
            exts = ['.nin', '.nsq', '.nhr']
        else:
            exts = ['.pin', '.psq', '.phr']

        for ext in exts:
            self.db_files.append(self.local_db_path + ext)
            self.hpc_job_files.append(self.local_db_path + ext)

        make_blast_index(self.local_db_path, self.db_type)

    def prep_bowtie(self):
        '''
            Prepares database and query files for quantification using bowtie.
            
            1) Makes a bowtie index of the file pointed to by 
                self.local_db_path
            2) Gets a list of the files created during step 1
            3) Creates a remote directory for the job 
            4) Appends the list of files from step 2 to hpc_job_files, so that
                they are moved to the remote computational resource
        '''
        #make and move the bowtie index
        make_bowtie_index(self.local_db_path, self.local_db_path)
        #get the list of all bowtie indicies
        proc = subprocess.Popen("ls " + self.local_dir + ' | grep -e ".ebwt$"',
                                shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        (out, err) = proc.communicate()
        out = out.split("\n")
        #make the remote directory
        c = netutils.ssh_connect(self.ssh_credentials)
        if self.special_run:
            make_remote_dir(settings.remote_dir + self.job_name, c)
        else:
            make_remote_dir(settings.remote_dir + self.run_name, c)
        c.close()
        for line in out:
            if(line.startswith("\n") or not ".ebwt" in line):
                continue
            line = line.rstrip("\n")
            local = self.local_dir + line
            self.db_files.append(local)
            self.hpc_job_files[local] = local

        if(self.paired_end):
            #need to combine the two files to a single, combined fastq file
            combined_fn = self.local_dir + self.pn + "_combined_for_quant.fastq"

            p = subprocess.Popen("cat " + self.input_files['query'] + " " +\
                                  self.input_files['query2'] + " > " +\
                                   combined_fn, shell=True)
            p.wait()
            self.hpc_qsub_files['query'] = combined_fn

    def start(self):
        self.prepare_for_aln[self.quant_type]()
        self.resources['modules'] = self.aligner_modules[self.quant_type]
        
        self.executable = self.pipeline_args.aligner
        HPCJob.start(self)


class Qsub:
    '''
        Class providing functionality to start and monitor qsub jobs on 
        remote computing resources utilizing a PBS system.
        
        Attributes:
            qsub (str):
                The full text of a the qsub script to be executed on the 
                remote cluster
                
            commands (str):
                String representation of commands to be executed via the 
                qsub script
                
            name (str):
                The unique name of this qsub object
                
            mem_req (int):
                The current memory request used when creating the job on the
                remote cluster
                
            walltime (str):
                The current walltime request used when creating the job
                
                Format: hh:mm:ss
                
            email (str):
                The email address to pass to the cluster for reporting job
                status
                
            q (str):
                The name of the q that the remote job should be submitted to
                
            nodes (int):
                The number of total nodes requested for the job
                
            ppn (int):
                The processor per node request for the job
                
            job_id (str):
                The identifier for this job on the remote cluster
            
            inputs ([str]):
                A list of input files for the job
                
            local_dir (str):
                Path to the local working directory for this job
               
            remote_dir (str):
                Path to the remote working directory for the job
                
            stats (str):
                Current status of the job
                
            finished_count (int):
                The number of times this job has been 'seen' as complete by 
                this object
                
            current_output (str):
                The name of the currently-used output file
    '''

    def __init__(self, n, d, inputs, local_dir, out_file = None):
        '''
            n (str):
                The name to use for this qsub
                
            d (str):
                Path to the remote working directory for this qusb
                
            input ([placeholder:str]):
                A dictionary of placeholder -> paths to the input files 
                for this qsub
                
            local_dir (str):
                Path to the local working diretory for this qsub
            
            out_file (str):
                Optional path to the resired local output file for this qsub
        '''
        self.qsub = ""
        self.commands = ""
        self.name = n
        self.mem_req = 0
        self.walltime = ""
        self.email = ""
        self.q = ""
        self.nodes = 0
        self.ppn = 0
        self.job_id = ""
        self.inputs = inputs
        self.local_dir = local_dir
        if(self.local_dir == ""):
            self.local_dir = "./"
        else:
            self.local_dir = self.local_dir + "/"
        self.remote_dir = d
        self.status = "idle"
        self.finished_count = 0;
        if(out_file is None): 
            self.current_output = self.name + ".stdout"
        else:
            self.current_output = out_file

    def _get_connection(self, h, u, p):
        '''
            Returns a new SSHConnection created using the provided host (h),
            username (u),and password (p).
        '''
        return netutils.SSHConnection(h, username=u, password=p)

    def append_qsub_command(self, command):
        '''
            Appends command to the string of commands to be run by the qsub 
            stiored in self.commands.
        '''
        self.commands += command + "\n"

    def append_to_pbs(self, directive):
        '''
            Appends directive to the string stored in self.qsub.
        '''
        self.qsub += directive + "\n"

    def cleanup(self):
        '''
            Removes each file stored in self.inputs.
        '''
        os.remove(self.local_dir + self.name + ".qsub")
        for infile in self.inputs.values():
            try:
                os.remove(infile)
                d, f = os.path.split(infile)
            except:
                pass

    def create_pbs_header(self, walltime, q, mem, nodes, ppn, email):
        '''
            Sets self.qsub using the supplied arguments.
            
            walltime (str):
                The walltime this qsub should request for it's job
                
            q (str):
                The queue this qsub should submit its job to
                
            mem (int):
                The memory request for this job, in megabytes. 
                
            nodes (int):
                The total number of nodes this qsub should request for the job
                
            ppn (int):
                The processors per node to use for this job
                
            email (str):
                The email address that messages about job status should be sent
                to
        '''
        self.qsub = "#! /bin/bash\n";
        self.qsub += "#PBS -r n\n";
        self.qsub += "#PBS -N " + self.name + "\n";
        self.qsub += "#PBS -o " + self.name + ".stdout\n";
        self.qsub += "#PBS -e " + self.name + ".stderr\n";
        if(q == 'billed'): self.qsub += "#PBS -W group_list=billed\n"
        else: self.qsub += "#PBS -q " + q + "\n";
        self.qsub += "#PBS -m a\n";
        self.qsub += "#PBS -M " + email + "\n";
        self.qsub += "#PBS -l walltime=" + walltime + "\n";
        self.qsub += "#PBS -l pmem=" + str(mem) + "\n";
        self.qsub += "#PBS -l nodes=" + str(nodes) + ":ppn=" + str(ppn) + "\n";
        self.qsub += "\n"
        self.mem_req = int(mem[:-2])
        self.walltime = walltime
        self.email = email
        self.q = q
        self.nodes = nodes
        self.ppn = ppn

    def retrieve_output(self, c):
        '''
            Retrieves the remote output of the qsub - looks for the file at:
                self.remote_dir + "/" + self.current_output
                
            Returns the local path of the file it retrieved.
        '''
        #gets the output created by this Qsub job, returns the local path to the file
        #sys.stdout.write("Retrieving: " + self.remote_dir + "/" + self.current_output + "\n")
        try:
            rfile = self.remote_dir + "/" + self.current_output
            lfile = self.local_dir + self.current_output
            c.get(rfile, lfile)
        except IOError as e:
            sys.stdout.write("\nError " + str(e.errno) +  "  "  + e.strerror + rfile + "\n")
            sys.stdout.write("Unable to retrieve these rusults;\njob will need to be rerun after checking: " + self.remote_dir)
            return None
#         return self.local_dir + self.current_output
        return lfile

    def resubmit(self, c, memIncrease = 0):
        '''
            Resubmits a failed job to the cluster using the remote connection, 
            c, and increasing the memory by the amount specified in memIncrease.
        '''
        self.create_pbs_header(self.walltime, self.q, 
                               str(self.mem_req + memIncrease) + "mb",
                               self.nodes, self.ppn, self.email)
        self.submit(c)

    def set_dir(self, directory):
        '''
            Sets self.remote_dir = directory
        '''
        self.remote_dir = directory

    def set_mem_request(self, req):
        '''
            Sets self.mem_req = req
        '''
        self.mem_req = req

    def submit(self, c):
        '''
            Submits the job to the remote cluster.
            
            1) writes the qsub file to the local disk as a .qsub script
            2) moves the qsub script to the remote working directory for the
                job
            3) moves all input files to the remote working directory
            4) changes to the remote working directory
            5) executes the qsub with 'qsub <scriptname>'
            6) sets this qsub's status to running and finished_count to 0
        '''
        #write the qsub submit file
        self.write_qsub()
        #put this qsub script on the server
        retries = 10
        while (retries > 0):
              try:
                 c.put(self.local_dir + self.name + ".qsub", self.remote_dir + "/" + self.name + ".qsub")
                 break
              except:
                 sys.stderr.write(" c.put(" + self.local_dir + self.name + ".qsub", self.remote_dir + "/" + self.name + ".qsub)\n")
                 sys.stderr.write("Exception raised during qsub submission to cluster \
                     sleeping for two minutes and trying again.\n")
                 time.sleep(120)
                 retries -= 1
                 if(retries == 0):
                    raise
        #put input files on the server
        for infile in self.inputs.values():
            f = os.path.split(infile)[1]
            retries = 10
            while (retries > 0):
              try:
                 c.put(infile, self.remote_dir + "/" + f)
                 break
              except:
                 sys.stderr.write("c.put(" + infile, self.remote_dir + "/" + f + ")\n")
                 sys.stderr.write("Exception raised during input file submission to cluster \
                     sleeping for two minutes and trying again.\n")
                 time.sleep(120)
                 retries -= 1
                 if(retries == 0):
                    raise

        #start the job
        command = "cd " + self.remote_dir + ";"
        command += "qsub " + self.name + ".qsub"
        retries = 10
        sleep_time = 30
        while (retries > 0):
              job_id = c.execute(command)
              if settings.debug_pipe: print "job_id from qsub submission: ", job_id
              if (job_id == [] or settings.QSUB_OK not in job_id[0]):
                  if(retries > 0):
                     sys.stderr.write("error executing qsub, retrying command: " + command + "\n")
                     retries -= 1
                     time.sleep(sleep_time)
                     sleep_time = sleep_time * 2
                  else:
                      raise ValueError("error executing qsub, abandoning it " + command + "\n")
              else:
                  break
        self.job_id = job_id[0].strip()
        self.status = "running"
        self.finished_count = 0

    def write_qsub(self):
        '''
            Writes a qsub script file named self.local_dir + self.name + '.qsub'
            
            Writes self.qsub and then self.commands to the script, in that 
            order.
        '''
        out = open(self.local_dir + self.name + ".qsub", 'w')
        out.write(self.qsub)
        out.write(self.commands)
        out.close()

class HPCProcess(PipeProcess):
    '''
        This PipeProcess provides the functionality to start and manage remote
        analyses on a high-performance computing cluster. This class was
        designed with the idea of handling multiple cluster architectures with
        a single class, although it currently only supports those clusters 
        utilizing a PBS system.
        
        Attributes:
            input_files ({flag_name:path}):
                A dictionary of input file paths, keyed on the input flag
                
            hpc_job_files ({flag_name:path}):
                A dictionary of files that need to be moved to the remote 
                cluster for job execution, but should not be split during the 
                process of query splitting
                
            hpc_qusb_files ({flag_name:path}):
                A dictionary of files that need to be split and then moved to 
                the remote computational resource. Each split file is given
                to an individual qsub file to work on
                
            job_type (str):    
                The type of this job
                
            retries (int):
                The number of retries to use when attemping a command/action
                via SSH
                
            resources ({res_name:res_amt}):
                A dictionary of resources this job will use during its 
                execution
                
            job_name (str):
                A name for this job
                
            local_dir (str):
                The local working directory for this job
                
            executable (str): 
                The executable that this object is planning to use to do its 
                work. The executable may change over the life of the object, as
                dictated by the needs of the analysis. This value is used to 
                look up command templates in the executables table, prior to 
                command execution
                
            processes ([Process]):
                A list of processes that belong to this object
                
            hpc_command (str):
                The fully-formed remote command to be executed at the HPC as
                a part of the qsub script
                
            current_output (str):
                Path to the current output file for this process
                
            mem_increment (int):
                The amount to increment the memory request when resubmitting
                a job that failed due to exceeding its curent request
                
            mail_credentials (settings.Credentials):
                A credentials object used to log in to an imap-enabled mail 
                server. This account accessed on this server should be set up
                to receive job status updates from the HPC cluster
            
            ssh_credentials (settings.Credentials):
                A credentials object used to connect to the remote cluster via
                SSH
            
            restart (Boolean):
                Whether jobs should be resubmitted upon failure 
                
            mail_connect:
                Holds open connections to the mail server
                
            remote_connect:
                Holds open ssh connection to the cluster    
    '''
    def __init__(
                 self,
                 input_files,
                 hpc_job_files,
                 hpc_qsub_files,
                 job_type,
                 resources,
                 remote_dir,
                 job_name,
                 local_dir,
                 executable,
                 process_args,
                 hpc_command,
                 current_output,
                 special_run,
                 jid,
                 output_file
                 ):
        Process.__init__(self)
        self.input_files = input_files
        self.hpc_job_files = hpc_job_files
        self.hpc_qsub_files = hpc_qsub_files
        self.job_type = job_type
        self.retries = 5
        self.resources = resources
        self.remote_dir = remote_dir
        self.job_name = job_name
        self.local_dir = local_dir
        self.executable = executable
        self.processes = []
        self.process_args = process_args
        self.hpc_command = hpc_command
        self.current_output = current_output
        self.special_run = special_run
        self.jid = jid
        self.outputFile = output_file
        self.mem_increment = 200
        self._adjust_cpu4_ppn()
        self._adjust_mem4_ppn()
        self.mail_credentials = settings.mail_cred
        self.ssh_credentials = settings.hpc_cred
        self.location = Locations.HPC
        self.retry_interval = 30
        self.restart = settings.restart_jobs
        self.mail_connect = None
        self.remote_connect = None
        self.mail_enabled = True

    def _adjust_cpu4_ppn(self):
        '''
            Adjusts the cpu request in self.resources['cpu'] by dividing the
            current value by self.resources['ppn']
        '''
        if (self.resources['cpu'] < self.resources['ppn']):
          print "WARNING: in _adjust_cpu4_ppn self.resources[cpu] < resources[ppn] will result in resources = 0"
        self.resources['cpu'] = self.resources['cpu'] / self.resources['ppn']

    def _adjust_mem4_ppn(self):
        '''
            On the HPC cluster, memory requests are given 'per node'. This
            method adjusts the memory requirement given to each qsub script
            by dividing self.resources['mem'] by self.resources['ppn']
        '''
        mem = self.resources['mem'].replace("mb", "")
        mem = int(mem)/self.resources['ppn']
        self.resources['mem'] = str(mem) + "mb"
        self.mem_increment = self.mem_increment / self.resources['ppn']

    def check(self, conn = None):
        '''
            Checks the status of this analysis process. 
            
            The method checks the status of each qsub job submitted by this
            object. If all qsub jobs have finished on the cluster, check returns
            JobState.FINISHED. Otherwise it returns JobState.RUNNING
        '''
        #create connections for the mailServer and reomte client for this check
        self.mail_connect = self._mail_connect()
        c = netutils.ssh_connect(self.ssh_credentials)
        num_finished = 0
        errors = []
        num_procs = len(self.processes)
        for index, proc in enumerate(self.processes):
            stat = self.check_individual_job(index, c)
            if(stat == 'finished'):
                proc.finished_count += 1
                #give error messages time to reach the mail server
                if(proc.finished_count * settings.mainLoopSleepInterval >=\
                    settings.waitBeforeMarkComplete):
                    num_finished += 1
                    proc.status ="processed"
            elif(stat == 'exceedMem'):
                #check if we should restart this job
                if(self.restart):
                    print "in stat=exceedMem' self.mem_increment: ", self.mem_increment
                    proc.resubmit(c, self.mem_increment)
                else:
                    num_finished += 1
                    proc.status = "error"
            elif(stat == "error"):
                errors.append(index)
            index += 1

        for i in reversed(errors):
            self.processes.pop(i)
        self.mail_connect.close()
        self.mail_connect.logout()
        c.close()

        if(len(self.processes) == num_finished):
            return JobState.FINISHED
        else:
            return JobState.RUNNING

    def _cleanup_local_files(self):
        '''
            Iterates over all processes in self.processes and calls cleanup()
            on each of them.
        '''
        #call each qsub job's cleanup
        for qsub in self.processes:
            qsub.cleanup()

    def _convert_dict_remotepaths(self, file_dict):
        ''' file_dict (dict):
                A dictionary of flag->file, or flag->file_list pairs. 
                
            Each path present in the dictionary is convereted to the remote path
            expected by the HPC cluster.

        '''
        ret = {}
        for flag, f in file_dict.items():
            if(isinstance(f, (list, tuple))):
                tmp = []
                for path in f:
                    tmp.append(self._convert_to_remotepath(path))
                ret[flag] = tmp
            elif(isinstance(f, basestring)):
                ret[flag] = self._convert_to_remotepath(f)
            else:
                raise Exception("Unsupported type given as value for file \
                                flag: " + str(flag))

        return ret

    def _convert_to_remotepath(self, f):
        '''
            Converts the local path, f, into a remote path. Returns the remote
            path.
        '''
        fn = os.path.split(f)[1]
        return self.remote_dir + fn

    def _rm_working_dir(self):
        '''
            Removes this job's working directory on the remote computational
            resource.
        '''
        connection = netutils.ssh_connect(self.ssh_credentials)
        command = "rm -r " + self.remote_dir
        connection.execute(command)
        connection.close()

    def _replace_using_dict(self, s, file_dict):
        '''
            Iterates over the key, value pairs in file_dict. Replaces each 
            instance of key in s with value. 
            
            Returns the modified s.
        '''
        for key, value in file_dict.items():
            s = s.replace("<" + key + ">", value)
        return s

    def _disable_mail(self):
        '''
            Prints an error message and sets mail_enabled = False.
        '''
        print("Error accessing mail account after " + str(self.retries) + "\
                retries, disabling mail checking for: " + str(self.job_name) \
                 + "/n")
        self.mail_enabled = False

    def _mail_connect(self):
        '''
            Returns an active connection to an IMAP-enabled mail server.
            
            Tries to create the connection maximum self.retries times. After
            self.retries attempts, disables mail checking for this object and
            returns None.
        '''
        retries = self.retries
        while(True and self.mail_enabled):
            try:
                mail_connect = imaplib.IMAP4_SSL(self.mail_credentials.host,
                                                 '993')
                mail_connect.login(self.mail_credentials.user,
                                   self.mail_credentials.passwd)
                mail_connect.select(settings.checkedMailbox)
                return mail_connect
            except Exception as e:
                if(retries > 0):
                    sys.stderr.write("Exception in _mail_connect() for " + self.job_name + "\n")
                    sys.stderr.write(e.message + "\n")
                    retries -= 1
                    sys.stderr.write("retries left: " + str(retries) + "\n")
                    time.sleep(self.retry_interval)
                else:
                    sys.stderr.write("invoking self._disable_mail\n")
                    self._disable_mail()
                    return None

        return None

    def _exceeded_memory(self, job):
        '''
            Checks if job exceeded its memory request during execution on
            the remote cluster. If it did, this method returns the amount of
            memory the job was using when it was killed. If this job did not
            exceed it's memory request, returns 0.
        '''
        retries = self.retries
        while(True and self.mail_enabled):
            try:
                typ, msg_ids = self.mail_connect.search(None, '(SUBJECT "PBS\
                 JOB ' + str(job.job_id) + '" BODY "exceeded MEM\
                  usage hard limit")')

                if(typ == "OK"):
                    if(msg_ids[0] != ''):
                        typ, data = self.mail_connect.fetch(msg_ids[0]\
                                                        ,'(BODY.PEEK[TEXT])')
                        return int(self._parse_email_body(data, "mem"))
                return 0

            except:
                if(retries > 0):
                    sys.stderr.write("Exception raised during mail access -- retrying!")
                    retries -= 1
                    #sleep, then start a new connection to mail server
                    time.sleep(self.retry_interval * 5)
                    self.mail_connect = self._mail_connect()
                    self.mail_connect.select(settings.checkedMailbox)
                else:
                    sys.stderr.write("Exception raised during mail access -- giving up no more retries left!")
                    self.disableMail()
                    return 0

        return 0

    def _exit_with_error(self, job):
        '''
            Not implemented.
        '''
        return False

    def _parse_email_body(self, response, directive):
        '''
            Parses the email text found in response, looking for information
            as determined by the directive.
            
            Currently supports only the 'mem' directive, which instructs the 
            method to search response for evidence that the job exceeded its
            initial memory request and was killed.
            
            mem directive:
                returns the memory usage of the job when it was killed, or an
                empty string if the job was not killed due to memory issues
        '''
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

    def _split_files(self, file_dict):
        '''
            Splits the files found in file_dict.values() into multiple files. 
            
            The number of files is dependent on the number of cpus (nodes)
            requested on the remote computational resource. There should be
            one split file for every node requested.
            
            The method returns the number of split files, and a new dictionary
            where each original key from the dictionary now indexes a list of 
            of the split files instead.
        '''
        split_files = {}
        num_parts = 0
        if(len(self.input_files) == 0):
            split_files['input'] = split_input(self.input_file, self.job_type,
                                               self.resources['cpu'])
            num_parts = len(split_files['input'])
        else:
            for flag, f in file_dict.items():
                if(self.resources['cpu'] == 1):
                    split_files[flag] = [f]
                else:
                    split_files[flag] = split_input(f, self.job_type,
                                                    self.resources['cpu'])

                num_parts = len(split_files[flag])
        return (num_parts, split_files)

    def check_individual_job(self, job_index, c):
        '''
            Checks the job status of the job held at self.processes[job_index],
            using the SSH connection, c.
            
            Returns the job's status.
        '''
        proc = self.processes[job_index]
        #check if this job terminated due to exeeding memory requests
        exceeded = self._exceeded_memory(proc)
        if(exceeded > 0):
            print "check_individual_job()  exceeded: ", exceeded
            proc.status = "exceedMem"
            if(self.restart):
                #set the new memory request to the amount used when the job was killed
                proc.mem_req = exceeded
            return proc.status

        #check if this job exited with an error status ######IMPLEMENT THIS##############
        elif(self._exit_with_error(proc)):
            proc.status = "error"
            return proc.status

        command =  "qstat -u " + settings.hpc_user +  "| grep -c " + proc.job_id[:16]
        retries = 10
        sleep_time = 120
        while (retries > 0):
            try:
                running = int(c.execute(command)[0])
                break
            except:
                retries -= 1
                if(retries == 0):
                    sys.stderr.write("Error retrieving status of running job " + proc.name + " - abandoning.\nCommand: " + command + "\n")
                    running = -1
                else:
                    time.sleep(sleep_time)
                    sleep_time = sleep_time * 2
                    sys.stderr.write("Error retrieving status of running job " + proc.name + " -retrying after sleeping " + str(sleep_time) + "\nCommand: " + command + "\n")

        if(running == 0):
            proc.status = "finished"
        elif(running == -1):
            proc.status = "error"
        else:
            proc.status = "running"

        return proc.status

    def complete(self):
        #retrieve the output files
        self.retrieve_output()
        self.cleanup()

    def cleanup(self):
        self._cleanup_local_files()
#        self._rm_working_dir()

    def move_job_data(self):
        '''
            Moves files in self.hpc_job_files to the remote computational 
            resource.
        '''
        c = netutils.ssh_connect(self.ssh_credentials)
        #make the remote dir, just in case
        c.execute("mkdir " + self.remote_dir)
        for data_file in self.hpc_job_files.values():
            d, fn = os.path.split(data_file)
            c.put(data_file, self.remote_dir + fn)

        c.close()

    def retrieve_output(self):
        '''
            Retrieves the remote output of all processes associated with this
            job.
        '''
        out = open(self.current_output, 'w')
        c = netutils.ssh_connect(self.ssh_credentials)
        #iterate over the processes that are a part of this job, getting their outputs
        for proc in self.processes:
            single_f = proc.retrieve_output(c)
            if(not single_f is None):
                tmp = open(single_f, 'r')
                out.write(tmp.read())
                tmp.close()
                os.remove(single_f)
        c.close()
        out.close()
        return self.current_output

    def run(self):
        '''
            Performs the majority of the work for starting and monitoring
            individual remote HPC cluster jobs. 
            
            1) moves job files that do not need to be split to the remote 
                cluster
            2) converts the local paths for input files and job files to
                remote paths
            3) splits relevant input files (those in self.hpc_qsub_files)
            4) creates a dictionary of all remote paths for replacement
                into a command template
            5) creates a remote working directory for the job
            6) iterates over each split input file
                - gets a command template for each input file
                - replaces remote input file paths into template
                - replaces process arguments into template
            7) creates a Qsub object
            8) Qsub writes a qsub script and moves it to remote cluster
            9) qsub is submitted to the remote cluster and the job is started
            10) monitors status of all qsub jobs and waits for completion
        '''
        print("HPC Process Starting job: " + self.job_name)
        session = netutils.make_db_session()
        #move job-level (non-split) files to the remote cluster
        self.move_job_data()

        ###move this code block to a function at some point
        remote_input_files = self._convert_dict_remotepaths(self.input_files)
        remote_job_files = self._convert_dict_remotepaths(self.hpc_job_files)

        (num_qsubs, qsub_files) = self._split_files(self.hpc_qsub_files)
        remote_qsub_files = self._convert_dict_remotepaths(qsub_files)

        remote_input_files.update(remote_job_files)
        remote_input_files.update(remote_qsub_files)

        qsub_flags = remote_qsub_files.keys()

        #for multiple-HPC system support, have this look up in a dictionary by location for credentials to use
        c = netutils.SSHConnection(self.ssh_credentials.host,
                                   username=self.ssh_credentials.user,
                                   password=self.ssh_credentials.passwd)

        working_dir = self.local_dir

        #make the remote directory for this job
        make_remote_dir(self.remote_dir, c)

        exec_plus_argstr = self.executable + " " + self.process_args.arg_string
        discard, outfile = os.path.split(self.current_output)
        for i in range(0, num_qsubs):
            ''' now we have:
                    self.hpc_job_files: local paths of all job-level files
                    remote_job_files: remote_paths of all job-level files
                    self.hpc_qsub_files: local paths of all qsub-level files
                    remote_qsub_files: remote_paths of all qsub-level files

                    need to:
                    combine hpc and qsub remotes and use for _prepare_cmd
                    create a qsub
                    submit the qsub
            '''
            job_no = i + 1
            job_name = self.job_name + str(job_no)
            qsub_out =job_name + \
                FileExtensions.for_program_output(self.job_type)

            this_input = {
                          flag: fpath for flag,
                          fpath in remote_input_files.items() \
                          if flag not in qsub_flags
                          }
            this_input.update({flag: remote_input_files[flag][i] \
                               for flag in qsub_flags})

            self.hpc_command = self._get_cmd_template(self.executable,
                                                      self.job_type, session)
            #            self.hpc_command = self._prepare_cmd(self.executable, self.job_type, input_files=this_input, output_file=self.remote_dir + qsub_out)

            if(self.hpc_command is None):
                self.hpc_command = exec_plus_argstr
                self.hpc_command = self._replace_using_dict(self.hpc_command,
                                                            this_input)
                self.hpc_command = self.hpc_command.replace("<output>",
                                                            self.remote_dir +\
                                                             qsub_out)
            else:
                param_dict = self._proc_options_dict(self.job_type, session)
                param_dict.update(this_input)
                param_dict['output'] = self.remote_dir + qsub_out
                param_dict['num_threads'] = self.resources['ppn']
                self.hpc_command = self._prepare_cmd(self.hpc_command,
                                                     param_dict)
            #get the files needed by this qsub
            these_qsub_files = {flag: qsub_files[flag][i]
                                for flag in qsub_flags}

            #get un-replaced command args to check for the presence of custom output field
            command_args = self._get_cmd_template(self.executable,
                                                  self.job_type, session)
            if(command_args is None):
                command_args = exec_plus_argstr
            if("<output>" in command_args or "<output_hpc>" in command_args):
                qsub = Qsub(job_name, self.remote_dir, these_qsub_files,
                            working_dir, qsub_out)
            else:
                qsub = Qsub(job_name, self.remote_dir, these_qsub_files,
                            working_dir)
            qsub.create_pbs_header(self.resources['wall'],
                                   self.resources['queue'],
                                   self.resources['mem'],
                                   self.resources['nodesPerSubJob'],
                                   self.resources['ppn'],
                                   settings.MAIL_ACCOUNT + "@" + settings.MAIL_PROVIDER + "," + settings.MAIL_ACCOUNT2 + "@" + settings.MAIL_PROVIDER)
            qsub.append_qsub_command(self.resources['modules'])
            touch_file = self.remote_dir + "/start." + job_name
            tcommand = "touch " + touch_file
            qsub.append_qsub_command(tcommand)
            qsub.append_qsub_command(self.hpc_command)
            touch_file = self.remote_dir + "/done." + job_name
            tcommand = "touch " + touch_file
            qsub.append_qsub_command(tcommand)
#           c = ssh connections
            qsub.submit(c)
            self.processes.append(qsub)

        c.close()

        while (True):
            stat = self.check()
            if stat == JobState.FINISHED:
                break
            time.sleep(60)
        self.complete()
        print("HPC process " , self.job_name,   "done.")
