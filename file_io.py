'''
Created on Jun 27, 2012

@author: Mathew Citarella

This module contains methods and classes that Autonomics uses for reading, writing, and manipulating files.

'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Crypto.Cipher import AES
import hashlib, math, os, random, re, sys, struct


#http://eli.thegreenplace.net/2010/06/25/aes-encryption-of-files-in-python-with-pycrypto/
def decrypt_file(key, in_filename, out_filename=None, chunksize=24*1024):
    """ Decrypts a file using AES (CBC mode) with the
        given key. Parameters are similar to encrypt_file,
        with one difference: out_filename, if not supplied
        will be in_filename without its last extension
        (i.e. if in_filename is 'aaa.zip.enc' then
        out_filename will be 'aaa.zip')
    """
    if not out_filename:
        out_filename = os.path.splitext(in_filename)[0]

    with open(in_filename, 'rb') as infile:
        origsize = struct.unpack('<Q', infile.read(struct.calcsize('Q')))[0]
        iv = infile.read(16)
        decryptor = AES.new(key, AES.MODE_CBC, iv)

        with open(out_filename, 'wb') as outfile:
            while True:
                chunk = infile.read(chunksize)
                if len(chunk) == 0:
                    break
                outfile.write(decryptor.decrypt(chunk))

            outfile.truncate(origsize)


#http://eli.thegreenplace.net/2010/06/25/aes-encryption-of-files-in-python-with-pycrypto/
def encrypt_file(key, in_filename, out_filename=None, chunksize=64*1024):
    """ Encrypts a file using AES (CBC mode) with the
        given key.

        key:
            The encryption key - a string that must be
            either 16, 24 or 32 bytes long. Longer keys
            are more secure.

        in_filename:
            Name of the input file

        out_filename:
            If None, '<in_filename>.enc' will be used.

        chunksize:
            Sets the size of the chunk which the function
            uses to read and encrypt the file. Larger chunk
            sizes can be faster for some files and machines.
            chunksize must be divisible by 16.
    """
    if not out_filename:
        out_filename = in_filename + '.enc'

    iv = ''.join(chr(random.randint(0, 0xFF)) for i in range(16))
    encryptor = AES.new(key, AES.MODE_CBC, iv)
    filesize = os.path.getsize(in_filename)

    with open(in_filename, 'rb') as infile:
        with open(out_filename, 'wb') as outfile:
            outfile.write(struct.pack('<Q', filesize))
            outfile.write(iv)

            while True:
                chunk = infile.read(chunksize)
                if len(chunk) == 0:
                    break
                elif len(chunk) % 16 != 0:
                    chunk += ' ' * (16 - len(chunk) % 16)

                outfile.write(encryptor.encrypt(chunk))


def make_record(line, t):
    ''' line (str):
            Line of text from a file to convert to a record object. Assumes the text is of the type t.

        t (str):
            Type of record to create, currently only supports "SAM" type record

        Record Fields:
            SAM - ['query', 'flags', 'reference', 'ref_start', 'score', 'cigar', 'mate_reference', 'mate_ref_start', 'insert_size', 'query_seq', 'query_qual']
            
        Returns a Record object with attributes set depending on the type of the record being parsed, t. Currently takes a tab-delimited line and returns a Record representing a SAM object. Additional single-line record formats will be supported in the future.
    '''

    if(t == "sam" or t == "SAM"):
        fields = ['query', 'flags', 'reference', 'ref_start', 'score', 'cigar', 'mate_reference', 'mate_ref_start', 'insert_size', 'query_seq', 'query_qual']
        el = line.split("\t")
        record = Record()
        for i in range(0, len(fields)):
            setattr(record, fields[i], el[i])
        return record


#add encryption to this later - http://eli.thegreenplace.net/2010/06/25/aes-encryption-of-files-in-python-with-pycrypto/
def read_credentials(f):
    ''' f (str):
            Path to file containing the credentials.

        Reads credentials stored in f to a dictionary with the following keys: user, passwd, dbuser, dbpasswd.

    '''
    fh = open(f, 'r')
    ret = {"user": None, "passwd": None, "dbuser": None, "dbpasswd": None}
    for line in fh:
        line = line.rstrip("\n")
        if(":" in line):
            el = line.split(":")
            ret[el[0].strip()] = el[1].strip()
    return ret


def sha1_file(afile):
    ''' afile (str):
            Path to file you wish to checksum using the SHA1 hashing algorithm.

       Hashes the contents of afile and returns a string representation of the 
       hexdigest.

    '''
    fh = open(afile, 'r')
    hasher = hashlib.sha1()
    for line in fh:
        hasher.update(line)

    fh.close()
    return hasher.hexdigest()


def translate_seq_file(filePath, fmt):
    ''' filePath (str):
            Path to file containing sequences you wish to translate.

        fmt (str):
            Format of the file stored at filePath. Support formats are 'fastq' 
            and 'fasta'.

        Translates the sequences in a given file from DNA or RNA to protein in 
        all six reading frames. Writes output to a file named 
        filePath_translated.fa.
    '''

    f = open(filePath, 'r')
    translated = open(filePath + "_translated.fa", "w")
    for record in SeqIO.parse(f, fmt):
        #write each translated sequence to the new file
        for i in range(-3, 4):
            try:
                SeqIO.write(translate_seq(record, i), translated, "fasta")
            except Exception as e:
                print("Unable to translate sequence: " + record.id  + " - " + e.message)

    f.close()
    translated.close()
    return filePath + "_translated.fa"


def translate_seq(seq, frame):
    ''' seq (str):
            DNA or RNA sequence you would like to translate.

        frame (int):
            Frame in which to translate the provided sequence. 
            Range: [-3, -2, -1, 1, 2, 3]

    '''
    if (not frame is None):
        sequence = seq.seq
        if frame < 0: sequence = sequence.reverse_complement()
        sequence = sequence[abs(frame)-1:].translate()
        if len(sequence) == 0:
            sequence = Seq('X')

        return SeqRecord(sequence, id = seq.id + "_" + str(frame))

    return None


class Alignment:
    '''
        Base class for storing alignment details.

        Supports the following attributes:
            query_id
            query_name
            query_seq
            reference_id
            significance
            score
            query_hit_len
            reference_hit_len
            query_length
            reference_length
            aln_length
            aln_text
            database_length
            identity
            start['query']
            end['query']
            start['reference']
            end['reference']
            hsps[]
            direction
    '''

    query_id = ""
    query_name = ""
    query_seq = ""
    reference_id = ""
    significance = -1.0
    score = -1.0
    query_hit_len = None
    reference_hit_len = None
    query_length = 0
    reference_length = 0
    aln_length = 0
    aln_text = ""
    database_length = 0
    identity = 0.0
    start = {}
    end = {}
    hsps = []
    direction = ""

    def __init__(self):
        '''
        Constructor
        '''

    def great_alignment(self):
        '''
            Check to see if this alignment meets the criteria for a 'great' 
            alignment. Currently this means capturing either 90% of the query 
            or reference sequence, as well as having 90% identity over the 
            length of the match.

            To do: add support for the method to take arguments specifying 
            cutoffs.
        '''

        if(self.query_hit_len is None or self.query_length == 0):
            query_len_frac = 1
        else:
            query_len_frac = float(self.query_hit_len)/float(self.query_length)

        if(self.reference_hit_len is None or self.query_length == 0):
            reference_len_frac = 1
        else:
            reference_len_frac = float(self.reference_hit_len)/float(self.reference_length)

        if(self.identity > 90 and 
           (query_len_frac > .80 or reference_len_frac > .80) and 
           self.significance < 1e-04):
            return True


        return False


class BlastTabAlignment(Alignment):
    '''
        Subclass of Alignment that represents tab-based BLAST alignment output.
         Constructor takes a line of BLAST-tab output and initializes the 
         appropriate attributes for the object.
    '''

    key_columns = [0, 1]

    def __init__(self, l):
        Alignment.__init__(self)
        self.query_id = l[0]
        self.reference_id = l[1]
        self.identity = l[2]
        self.hit_length = l[3]
        self.start['query'] = float(l[6])
        self.end['query'] = float(l[7])
        self.start['reference'] = float(l[8])
        self.end['reference'] = float(l[9])
        self.significance = float(l[10])
        self.query_length = float(l[12])
        self.reference_length = float(l[13])
        self.query_hit_len = math.fabs(self.start['query'] - self.end['query']) + 1
        self.reference_hit_len = math.fabs(self.start['reference'] - \
                                           self.end['reference']) + 1


class BlatAlignment(Alignment):
    '''
        Subclass of Alignment that represents PSL-formatted output from the 
        BLAT alignment software.

        Adds the attributes:

            blockSizes:

                List of alignment block sizes produced by BLAT between the 
                query and the reference

            qstarts:

                List of start positions of the alignment blocks in the query
                 sequence.

            tStarts:

                List of start positions of the aligmment blocks in the 
                reference sequence.

            key_columns:

                List of columns from which a unique text key can be 
                constructed. Keys should be constructed such that they are 
                unique across all records present in the alignment file.
    '''

    blockSizes = []
    qStarts = []
    tStarts = []
    key_columns = [9, 13]

    def __init__(self, elements):
        Alignment.__init__(self)
        self.query_id = elements[9]
        self.reference_id = elements[13]
        self.query_length = float(elements[10])
        self.reference_length = float(elements[14])
        self.start['query'] = float(elements[11])
        self.end['query'] = float(elements[12])
        self.start['reference'] = float(elements[15])
        self.end['reference'] = float(elements[16])
        self.direction = elements[8]
        self.identity = float(1 - (float(elements[1])/float(elements[0]))) * 100
        self.query_hit_len = math.fabs(self.start['query'] - \
                                       self.end['query']) + 1
        self.reference_hit_len = math.fabs(self.start['reference'] - \
                                           self.end['reference']) + 1


class WUBlastAlignment(Alignment):
    '''
        Subclass of Alignment that represents tab-delimited output from the 
        WU-BLAST alignment software.
    '''
    key_columns = [0, 1]

    def __init__(self, elements):
        Alignment.__init__(self)
        self.query_id = elements[0]
        self.reference_id = elements[1]
        self.signifncance = float(elements[2])
        self.identity = float(elements[10])
        self.start['query'] = float(elements[17])
        self.end['query'] = float(elements[18])
        self.start['reference'] = float(elements[20])
        self.end['reference'] = float(elements[21])


class AlignmentReader:
    '''
        This class contains the functionality to parse a variety of alignment 
        files and internally store a dictionary of annotations.

        --Supported Alignment Formats [RecordClass]--
            blast-tab [BlastTabAlignment]
            hpc-blast [BlastAlignment]
            wu-blast [WUBlastAlignment]
            sam [Alignment]
            psl [BlatAlignment]

        --Data Structures--

            annotations: Dictionary of key->Alignment object pairs.
            unannotated: List of keys for unannotated sequences.

    '''

    file_name = ""
    fmt = ""
    annotations = {}
    unannotated = []
    databases = {'1':1, '2':2, 'nr':1, 'sp':2}

    def __init__(self, file_path, fmt, database = '1', verbose = False):
        ''' file_path:
                Path to the file to be read.

            fmt:
                Format of the file to be read. Can be one of: ['blast-tab', 
                'wu-blast', 'hpc-blast', 'psl', 'sam'].

            database:
                Usually not needed, 1 specifies an alignment against NR, 2 an
                 alignment against SwissProt.

            verbose:
                Print detailed progress information.

        '''
        self.file_name = file_path
        self.fmt = fmt
        self.v = verbose
        self.database = self.databases[str(database)]

    def read(self):
        '''
            Calls internal parse methods, depending on the format attribute of
            this AlignmentReader.

            Result of this call is population of the annotations and 
            unannotated attributes of the AlignmentReader object.

        '''
        if(self.fmt == 'sam'):
            self._parseSAM()
        elif(self.fmt == 'wu-blast'):
            self._parse_tab_file(WUBlastAlignment)
        elif(self.fmt == 'blast-tab'):
            self._parse_tab_file(BlastTabAlignment)
        elif('blast' in self.fmt):
            self._parse_blast()
        elif(self.fmt == 'psl'):
            self._parse_tab_file(BlatAlignment)

    def _parse_blast(self):
        a = open(self.file_name, 'r')
        state = "START"
        query_reg = re.compile("^Query=\s+")
        if('mpi' in self.fmt): 
            length_reg = re.compile("\((\d+\,)*\d+\s+letters\)")
        else: 
            length_reg = re.compile("^Length=")
        if('mpi' in self.fmt): 
            hit_length_reg = re.compile("\s+Length\s=")
        else: 
            hit_length_reg = length_reg
        new_hit_reg = re.compile("^>")
        no_hits_reg = re.compile("^\*{5} No hits found")
        score_reg = re.compile("Score\s\=")
        iden_reg = re.compile("(\(\d{2,3}\%\))")
        subj_strand = re.compile("(Sbjct)\:*\s+\d+")
        query_strand = re.compile("(Query)\:*\s+\d+")
        reg_spaces = re.compile("\s+")
        if('mpi' in self.fmt): 
            effSearchSpace = re.compile("^\s+Number of sequences in database")
        elif("hpc" in self.fmt): 
            effSearchSpace = re.compile("^Effective search space used")
        else: 
            effSearchSpace = re.compile("^Reference: Altschul")
        annotation = Alignment()
        first_score = True
        query_name = ""
        examined = 0
        after_query = False
        if(self.v):print("Parsing BLAST annotations in: " + self.file_name)
        for l in a:
            l = l.rstrip('\n')
            if(state == "START"):
                if(query_reg.search(l)):
                    examined += 1
                    if(self.v):
                        sys.stderr.write("%d query sequeneces processed.\r" \
                                         % examined)
                    state = "QUERY"
                    query_name = query_reg.sub("", l)
                    first_score = True
            elif(state == "QUERY"):
                if(length_reg.search(l)):
                    state = "HIT"
                else:
                    query_name = query_name + l

            elif(state == "HIT"):
                if(new_hit_reg.search(l)):
                    #this is a new hit for this query sequence
                    hit_id = new_hit_reg.sub("", l)
                    annotation = Alignment()
                    annotation.query_id = query_name,
                    annotation.reference_id = hit_id
                    self.annotations[query_name] = []
                    state = "HIT_DESCRIPTION"
                    #add the annotation text to the appropriate attribute
                    annotation.aln_text = l + "<br />"
                elif(no_hits_reg.search(l)):
                    #there were no hits for this query sequence
                    self.unannotated.append(query_name)
                    state = "START"

            elif(state == "HIT_DESCRIPTION"):
                #append this line to the growing annotation text
                annotation.aln_text += l + "<br/>"

                if(hit_length_reg.search(l) and first_score == True):
                    state = "WAITING_FOR_SCORE"
                else:
                    #continue parsing reference id
                    annotation.reference_id = annotation.reference_id + l

            elif(state == "WAITING_FOR_SCORE"):
                annotation.aln_text += l + "<br/>"
                if(score_reg.search(l)):
                    #get the evalue
                    elements = re.split("\=\s+", l)
                    elements[2] = re.sub(",.+", "", elements[2])
                    if(elements[2].startswith("e")):
                        elements[2] = "1" + elements[2]
                    annotation.significance = float(elements[2].strip())
                    #get the score
                    score = re.split("\(", elements[1])
                    score = re.split("\s", score[0])
                    annotation.score = float(score[0])
                    state = "GOT_SCORE"

            elif(state == "GOT_SCORE"):
                #append this line to the growing annotation text
                annotation.aln_text += l + "<br/>"
                #get the % identity
                if(first_score == True):
                    matches = iden_reg.search(l)
                    annotation.identity = float(re.sub("(\(|\)|\%)", "", 
                                                       matches.group(1)))
                state = "GOT_IDENT"

            elif(state == "GOT_IDENT"):
                #append this line to the growing annotation text
                annotation.aln_text += l + "<br/>"
                #get the direction
                if(first_score == True):
                    if(re.search("Plus\ Minus", l)):
                        annotation.direction = 'antisense'
                    else:
                        annotation.direction = 'sense'

                    first_score = False

                state = "ALIGNMENTS"

            elif(state == "ALIGNMENTS"):
                #start of the next hit
                if(new_hit_reg.search(l)):
                    self.annotations[query_name].append(annotation)
                    hit_id = new_hit_reg.sub("", l)
                    annotation = Alignment()
                    annotation.query_id = query_name
                    annotation.reference_id = hit_id
                    annotation.aln_text = l + "<br/>"
                    first_score = True
                    state = "HIT_DESCRIPTION"

                elif(effSearchSpace.search(l)):
                    self.annotations[query_name].append(annotation)
                    annotation = None
                    state = "START"

                elif(query_strand.search(l)):
                    after_query = True
                    elements = reg_spaces.split(l)
                    annotation.aln_text += "<table style=\"margin-top: 20px\">\
                        <tr><td>" + str(elements[0]) + " " + str(elements[1]) +\
                         "</td><td>" + str(elements[2]) + "</td><td>" +\
                          str(elements[3]) + "</td></tr>"
                elif(subj_strand.search(l)):
                    elements = reg_spaces.split(l)
                    annotation.aln_text += "<tr><td>" + str(elements[0]) + " "\
                     + str(elements[1]) + "</td><td>" + str(elements[2]) + \
                     "</td><td>" + str(elements[3]) + "</td></tr>" + "</table>"
                else:
                    if(after_query):
                        after_query = False
                        l = reg_spaces.sub("&nbsp;", l)
                        annotation.aln_text += "<tr><td></td><td>" + l + \
                            "</td><td></td></tr>"
                    else:
                        annotation.aln_text += l + "<br/>"

        if(self.v): print("Done.")
        return self.annotations

    def _parse_tab_file(self, record_maker):
        theseAnnotations = []
        thisQuery = None
        fh = open(self.file_name, 'r')
        for l in fh:
            record = record_maker(l.split("\t"))
            if(not thisQuery is None):
                #check if the query matches the current query
                if(record.query_id != thisQuery):
                    #store the annotations, reset the annotation list and 
                    #store the new queryID
                    self.annotations[thisQuery] = theseAnnotations
                    thisQuery = record.query_id
                    theseAnnotations = []
                theseAnnotations.append(record)
            else:
                thisQuery = record.query_id
                theseAnnotations.append(record)
        self.annotations[thisQuery] = theseAnnotations
        return self.annotations

    def _parseSAM(self):
        pass


class FileExtensions:
    '''
        This class contains data structures for supported input and output 
        extensions, and methods to access the data structures. Details of the 
        data structures are provided below.

        --Class Attributes--
            input_exts: A dictionary of 'job_type' -> 'job_input_extension' 
                pairs. Jobs look up their input extensions when determining the 
                default input file to use.
            
            output_exts: A dictionary of 'job_type' -> 'job_output_extension' 
                pairs. Jobs look up their own output extensions when 
                automatically determining how to name their output files.

    '''

    input_exts = {"panther": ".fasta", "quality_trim": ".fastq", 
                  "adapter_trim": ".fastq", "read_normalization": ".fastq", 
                  "blast": "_project.fasta", 
                  "blast_swissprot": "_project.fasta", 
                  "blast_nr": "_project.fasta", "blat" : "_project.fasta", 
                  "pfam": "_project.fasta", "go": "_blast_swissprot.txt", 
                  "kegg": "_blast_swissprot.txt"}
    output_exts = {"panther": "_panther", "assemble": "_project.fasta", 
                   "blast": "_blast.txt", 
                   "blast_swissprot": "_blast_swissprot.txt", 
                   "blast_nr": "_blast_nr.txt", "blat" : "_blat.psl", 
                   "kegg": "_KEGG.txt", "go": "_GO.txt", 
                   "quality_trim": ".fastq", "adapter_trim": ".fastq", 
                   "read_normalization": ".fasta", "pfam": "_pfam.txt", 
                   "quantification": "_quantification.txt", 
                   "go_categories": "_gocats.txt"}

    @staticmethod
    def for_input_to_prog(program):
        ''' program:
                Name of an executable for which an input extension should be 
                returned.
        '''
        return FileExtensions.input_exts[program]

    @staticmethod
    def for_input(jobType):
        ''' jobType:
                Job type for which an input extension should be looked up and 
                returned.

        '''
        return FileExtensions.input_exts[jobType]

    @staticmethod
    def for_output(jobType):
        ''' jobType:
                Same as 'for_input', but retrieves output extensions instead.

        '''
        return FileExtensions.output_exts[jobType]

    @staticmethod
    def for_program_output(program):
        ''' program:
                Same as 'for_input_to_prog', but returns an output extension 
                instead.

        '''
        return  FileExtensions.output_exts[program]


class Record:

    '''
        This class represents an empty Record object, to which arbitrary 
        attributes can be assigned depending on the record type, at run-time.
    '''

    def __init__(self):
        pass


class FileFormats:
    '''
        Placeholder for a future feature implementation.
    '''


wu_blast_typeflags = {"NT": '-n', "AA": '-p'}
