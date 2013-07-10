'''
Created on Jun 27, 2012

@author: mat
'''
import math
import sys
import re
from files.io import TabFile

class Alignment:

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

    def great_alignment(self, significance_cutoff=1e-04, identity_cutoff=90):


        if(self.query_hit_len is None or self.query_length == 0):
            query_len_frac = 1
        else:
            query_len_frac = float(self.query_hit_len)/float(self.query_length)

        if(self.reference_hit_len is None or self.query_length == 0):
            reference_len_frac = 1
        else:
            reference_len_frac = float(self.reference_hit_len)/float(self.reference_length)

        if(self.identity > identity_cutoff and (query_len_frac > .80 or reference_len_frac > .80) and self.significance < significance_cutoff):
           return True

        return False


class BlastTabAlignment(Alignment):

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
        self.reference_hit_len = math.fabs(self.start['reference'] - self.end['reference']) + 1


class BlatAlignment(Alignment):

    blockSizes = []
    qstarts = []
    tStarts = []
    key_columns = [9, 13]

    def __init__(self, elements):
        Alignment.__init__(self)
        self.query_id = elements[9]; self.reference_id = elements[13]; self.query_length = float(elements[10]); self.reference_length = float(elements[14]); self.start['query'] = float(elements[11]); self.end['query'] = float(elements[12]); self.start['reference'] = float(elements[15]); self.end['reference'] = float(elements[16]); self.direction = elements[8]; self.identity = float(1 - (float(elements[1])/float(elements[0]))) * 100
        self.query_hit_len = math.fabs(self.start['query'] - self.end['query']) + 1
        self.reference_hit_len = math.fabs(self.start['reference'] - self.end['reference']) + 1


class WUBlastAlignment(Alignment):

    key_columns = [0, 1]

    def __init__(self, elements):
        Alignment.__init__(self)
        self.query_id = elements[0]
        self.reference_id = elements[1]
        self.significance = float(elements[2])
        self.identity = float(elements[10])
        self.start['query'] = float(elements[17])
        self.end['query'] = float(elements[18])
        self.start['reference'] = float(elements[20])
        self.end['reference'] = float(elements[21])


class Reader:

    fileName = ""
    fmt = ""
    annotations = {}
    unannotated = []
    databases = {'1':1, '2':2, 'nr':1, 'sp':2}

    def __init__(self, fileName, fmt, database = '1', verbose = False):
        '''
        Constructor
        '''
        self.fileName = fileName
        self.fmt = fmt
        self.v = verbose
        self.database = self.databases[str(database)]

    def read(self):
        '''
        MODIFY THIS FUNCTION TO RETURN AN ALIGNMENTRESULT OBJECT, DO IT, SLACKER!!!!
        The object should hold the annotations, unannotated records, database, alignment program used
        '''
        '''
        Code to parse the annotation file associated with this parser
        '''
        if(self.fmt == 'sam'):
            self.parseSAM()
        elif(self.fmt == 'wu-blast'):
            self._parse_tab_file(WUBlastAlignment)
        elif(self.fmt == 'blast-tab'):
            self._parse_tab_file(BlastTabAlignment)
        elif('blast' in self.fmt):
            self._parseBLAST()
        elif(self.fmt == 'psl'):
            self._parse_tab_file(BlatAlignment)

    def _parseBLAST(self):
        a = open(self.fileName, 'r')
        state = "START"
        queryReg = re.compile("^Query=\s+")
        if('mpi' in self.fmt): lengthReg = re.compile("\((\d+\,)*\d+\s+letters\)")
        else: lengthReg = re.compile("^Length=")
        if('mpi' in self.fmt): hitLengthReg = re.compile("\s+Length\s=")
        else: hitLengthReg = lengthReg
        newHitReg = re.compile("^>")
        removeQuotes = re.compile("(\'|\")")
        noHitsReg = re.compile("^\*{5} No hits found")
        scoreReg = re.compile("Score\s\=")
        idenReg = re.compile("(\(\d{2,3}\%\))")
        subjectStrand = re.compile("(Sbjct)\:*\s+\d+")
        queryStrand = re.compile("(Query)\:*\s+\d+")
        regSpaces = re.compile("\s+")
        if('mpi' in self.fmt): effSearchSpace = re.compile("^\s+Number of sequences in database")
        elif("hpc" in self.fmt): effSearchSpace = re.compile("^Effective search space used")
        else: effSearchSpace = re.compile("^Reference: Altschul")
        annotation = Alignment()
        firstScore = True
        queryName = ""
        inHitDesc = 0
        examined = 0
        afterQuery = False
#        if(self.v):print("Parsing BLAST annotations in: " + self.fileName)
        for l in a:
            l = l.rstrip('\n')
            if(state == "START"):
                if(queryReg.search(l)):
                    examined += 1
#                    if(self.v):sys.stderr.write("%d query sequeneces processed.\r" % examined)
                    state = "QUERY"
                    queryName = queryReg.sub("", l)
                    firstScore = True
            elif(state == "QUERY"):
                if(lengthReg.search(l)):
                    state = "HIT"
                else:
                    queryName = queryName + l

            elif(state == "HIT"):
                if(newHitReg.search(l)):
                    #this is a new hit for this query sequence
                    hitId = newHitReg.sub("", l)
                    annotation = Alignment()
                    annotation.query_id = queryName,
                    annotation.reference_id = hitId
                    self.annotations[queryName] = []
                    state = "HIT_DESCRIPTION"
                    #add the annotation text to the appropriate attribute
                    annotation.aln_text = l + "<br />"
                elif(noHitsReg.search(l)):
                    #there were no hits for this query sequence
                    self.unannotated.append(queryName)
                    state = "START"

            elif(state == "HIT_DESCRIPTION"):
                #append this line to the growing annotation text
                annotation.aln_text += l + "<br/>"

                if(hitLengthReg.search(l) and firstScore == True):
                    state = "WAITING_FOR_SCORE"
                else:
                    #continue parsing reference id
                    annotation.reference_id = annotation.reference_id + l

            elif(state == "WAITING_FOR_SCORE"):
                annotation.aln_text += l + "<br/>"
                if(scoreReg.search(l)):
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
                if(firstScore == True):
                    matches = idenReg.search(l)
                    annotation.identity = float(re.sub("(\(|\)|\%)", "", matches.group(1)))
                state = "GOT_IDENT"

            elif(state == "GOT_IDENT"):
                #append this line to the growing annotation text
                annotation.aln_text += l + "<br/>"
                #get the direction
                if(firstScore == True):
                    if(re.search("Plus\ Minus", l)):
                        annotation.direction = 'antisense'
                    else:
                        annotation.direction = 'sense'

                    firstScore = False

                state = "ALIGNMENTS"

            elif(state == "ALIGNMENTS"):
                #start of the next hit
                if(newHitReg.search(l)):
                    self.annotations[queryName].append(annotation)
                    hitId = newHitReg.sub("", l)
                    annotation = Alignment()
                    annotation.query_id = queryName
                    annotation.reference_id = hitId
                    annotation.aln_text = l + "<br/>"
                    firstScore = True
                    state = "HIT_DESCRIPTION"

                elif(effSearchSpace.search(l)):
                    self.annotations[queryName].append(annotation)
                    annotation = None
                    state = "START"

                elif(queryStrand.search(l)):
                    afterQuery = True
                    elements = regSpaces.split(l)
                    annotation.aln_text += "<table style=\"margin-top: 20px\"><tr><td>" + str(elements[0]) + " " + str(elements[1]) + "</td><td>" + str(elements[2]) + "</td><td>" + str(elements[3]) + "</td></tr>"
                    #''.join([annotation.aln_text,"<table style=\"margin-top: 20px\">", "<tr><td>" + str(elements[0]) + " " + str(elements[1]) + "</td><td>" + str(elements[2]) + "</td><td>" + str(elements[3]) + "</td></tr>"])
                elif(subjectStrand.search(l)):
                    elements = regSpaces.split(l)
                    annotation.aln_text += "<tr><td>" + str(elements[0]) + " " + str(elements[1]) + "</td><td>" + str(elements[2]) + "</td><td>" + str(elements[3]) + "</td></tr>" + "</table>"
                    #''.join([annotation.aln_text, "<tr><td>" + str(elements[0]) + " " + str(elements[1]) + "</td><td>" + str(elements[2]) + "</td><td>" + str(elements[3]) + "</td></tr>", "</table>"])
                else:
                    if(afterQuery):
                        afterQuery = False
                        l = regSpaces.sub("&nbsp;", l)
                        annotation.aln_text += "<tr><td></td><td>" + l + "</td><td></td></tr>"
                    else:
                        annotation.aln_text += l + "<br/>"

#        if(self.v): print("Done.")
        return self.annotations

    def _parse_tab_file(self, record_maker):
        theseAnnotations = []
        thisQuery = None
        fh = open(self.fileName, 'r')
        for l in fh:
            record = record_maker(l.split("\t"))
            if(not thisQuery is None):
                #check if the query matches the current query
                if(record.query_id != thisQuery):
                    #store the annotations, reset the annotation list and store the new queryID
                    self.annotations[thisQuery] = theseAnnotations
                    thisQuery = record.query_id
                    theseAnnotations = []
                theseAnnotations.append(record)
            else:
                thisQuery = record.query_id
                theseAnnotations.append(record)
        self.annotations[thisQuery] = theseAnnotations
        return self.annotations

    def parseSAM(self):
        pass
