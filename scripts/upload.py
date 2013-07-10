#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      mat
#
# Created:     08/08/2012
# Copyright:   (c) mat 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from alignment.io import Reader
from Bio.Blast import NCBIStandalone
from Bio import SeqIO
from sqlalchemy import *
from sqlalchemy.sql import *
from autonomics.file_io import Record
from autonomics import settings, netutils, utility
import argparse
import datetime
import os
import MySQLdb
import re
import sqlalchemy
import subprocess
import sys
import time
import utils
from datetime import date
import warnings

perlPath=settings.PERLPATH
baseDir = settings.NEUROBASE_DATA_PATH
debug = 0
debug_kegg = 0

def allPipe(args, s, v):

    s.activateConnection()
    session = s
    #ftp directory
    baseDir = settings.NEUROBASE_DATA_PATH
    if(not args.data_dir is None):
        baseDir = args.data_dir
        
    baseName = baseDir + args.pn + "/" + args.pn

    reportStatus("Initiating project...\n", v)
    initProject(args.publicName, s)
    reportStatus("done initiating project.\n", v)

    reportStatus("Loading project sequences...\n", v)
    seqs = baseName + "_project.fasta"
    loadSequences(args.pn, args.seqType, s, seqs, v, sbStart=args.sbStart)
    reportStatus("done loading seqs.\n", v)

    projectID = getProjectID(args.pn, s, True)

    reportStatus("Linking sequence identifiers in FASTA file to NeuroBase ID...\n", v)
    linkDict = utils.link_dbid_to_fastaid(projectID, s)

    if (quant):
      reportStatus("Loading quantification data...\n", v)
      if(args.quantFile is None):
          args.quantFile = baseName + "_quantification.txt"
      loadQuantification(projectID, args.pn, args.quantFile, args.abundanceCol, s, v)
      reportStatus("done loading quantif.\n", v)
    
    reportStatus("Loading swissprot homology data...\n", v)
    p = Reader(baseName + "_blast_swissprot.txt", args.alnFmt, 2, v)
    p.read()
    loadHomology(projectID, linkDict, p, s, args.deleteH)
    linfo = utils.createTableObject('load_info', session)
    insert = linfo.insert()
    insert.execute(project_id=projectID, data='swiss')

    reportStatus("done loading swiss.\n", v)
    
    reportStatus("Storing annotation alignments...\n", v)
    loadHomologyAlignments(projectID, linkDict, p, s, v)
    reportStatus("done storing anno.\n", v)
    linfo = utils.createTableObject('load_info', session)
    insert = linfo.insert()
    insert.execute(project_id=projectID, data='anno_align')

    reportStatus("Assigning best annotations to sequences..\n", v)
    assignBestAnnotation(projectID, s, v)
    reportStatus("done assign best\n", v)

    reportStatus("Loading GO terms...\n", v)
    if(args.goFile is None):
        args.goFile = baseName + "_GO.txt"
    loadGO(projectID, linkDict, args.goFile, s, v)
    reportStatus("done loading go.\n", v)
    linfo = utils.createTableObject('load_info', session)
    insert = linfo.insert()
    insert.execute(project_id=projectID, data='go')

    reportStatus("Loading GO categories..\n", v)
    catFile = baseName + "_gocats.txt"
    loadGOCategories(projectID, args.pn, catFile, s, v)
    reportStatus(".done loading gocats\n", v)
    linfo = utils.createTableObject('load_info', session)
    insert = linfo.insert()
    insert.execute(project_id=projectID, data='gocat')
    
    reportStatus("Loading KEGG pathways...\n", v)
    if(args.keggFile is None):
        args.keggFile = baseName + "_KEGG.txt"
    loadKEGG(projectID, linkDict, args.keggFile, s, v)
    reportStatus("done loading kegg.\n", v)
    
    reportStatus("Loading pfam annotations...\n", v)
    if(args.pfamFile is None):
            args.pfamFile = baseName + "_pfam.txt"
    loadPfam(projectID, args.pfamFile, args.seqType, linkDict, s)
    reportStatus("done loading pfam.\n", v)
    linfo = utils.createTableObject('load_info', session)
    insert = linfo.insert()
    insert.execute(project_id=projectID, data='pfam')


def assignBestAnnotation(projectID, session, v):

    if(v):
        reportStatus("Assigning best annotation.\n", v)

    session.activateConnection()
    sorted_homology = netutils.get_table_object("sorted_homology", session)
    annotation_db = netutils.get_table_object('annotation_db', session)
    seq_table = netutils.get_table_object(str(projectID) + "_sequences", session)
    best_annot = netutils.get_table_object('best_annotations', session)
    datab = 'best_anno'
    session.conn.execute("DELETE FROM load_info WHERE data='" + datab + "' and project_id='" + str(projectID) + "'")

    res = seq_table.select().execute()
    for row in res.fetchall():
        conn = MySQLdb.connect(host=session.host, db=session.db,
                           user=session.user, passwd=session.passwd)
        cursor = conn.cursor(MySQLdb.cursors.DictCursor)

        sb_id = row.sb_id
        
        q = "".join(["SELECT text, evalue FROM sorted_homology JOIN annotation_db",
                     " on sorted_homology.annotation_id=annotation_db.annotation_id",
                     " AND project_id=%s AND sb_id=%s AND sort_id=1 ORDER BY ranking",
                     " LIMIT 0, 1"])
        cursor.execute(q, (projectID, sb_id))

        annot_row = cursor.fetchone()
        if(annot_row is None):
            annot_row = {'text': "No annotation",
                         'evalue': 0}
        i = best_annot.insert().values(project_id=projectID, sb_id=sb_id, annot=annot_row['text'], 
                                       eval=annot_row['evalue'])
        try:
            session.conn.execute(i)
        except sqlalchemy.exc.IntegrityError as e:
            u = best_annot.update().where(and_(best_annot.c.project_id==projectID, 
                                               best_annot.c.sb_id==sb_id)).values(annot=annot_row['text'], 
                                                                                  eval=annot_row['evalue'])
            session.conn.execute(u)
        
    linfo = utils.createTableObject('load_info', session)
    insert = linfo.insert()
    insert.execute(project_id=projectID, data='best_anno')


def collapseAssociationFile(lines, keyColumns, rankColumn, keep):

    ret = {}
    for line in lines:
        key = "-".join(line[c] for c in keyColumns)
        if(key in ret):
            if(keep == "largest"):
                if(float(line[rankColumn]) > ret[key][rankColumn]):
                    ret[key] = line
            else:
                if(float(line[rankColumn]) < ret[key][rankColumn]):
                    ret[key] = line
        else:
            ret[key] = line

    return ret

# not used
def getGOParents(go_id, session):
    go_relations = utils.createTableObject("go_relationships", session)
    results = go_relations.select(go_relations.c.child_id==go_id).execute()
    return [row.parent_id for row in results]

def getProjectID(projectName, session, v):
    t = utils.createTableObject('project_directory', session)
    reportStatus("Getting the projectID from the project name...\n", v)
    select = t.select(t.c.project_name == projectName)
    row = select.execute().fetchone()
    if(row): 
        projectID = row.projectID
        reportStatus("It's " + str(projectID) + "\n", v)
        return projectID
    else: 
        return None

def linkDescToID(desc, linkDict, attempt = 0):
    sbID = None
    while(sbID is None):
        try:
            sbID = linkDict[desc]
        except:
            attempt += 1
            if(attempt < 2):
                desc = repairDescription(desc, attempt)
            else:
                    return None

    return sbID


def initProject(publicName, session):
    #create empty project directories, project_directory entry, and sequence table
    #first, check if this project already exists within the database
    t = utils.createTableObject('project_directory', session)
    results = t.select(t.c.project_name == publicName).execute()
    if(not results.fetchone()):
        #set up the directory structure for this project
        session.activateConnection()
        results = session.conn.execute(select([func.max(t.c.projectID).label('maxID')]))
        row = results.fetchone()
        if(row is None):
            cur_id = 0
        else:
            cur_id = row.maxID
            
        if(cur_id is None):
            cur_id = 0
            
        projectID = int(cur_id) + 1
        projectPath = settings.NEUROBASE_FASTADB_PATH + str(projectID)
        if(not os.path.exists(projectPath)):
            os.mkdir(projectPath)
        #create a sequence table for this project
        seqTable = utils.create_seq_table(str(projectID) + "_sequences", session)
        seqTable.create()
        #insert an entry in the project_directory for the new project
        insert = t.insert()
        insert.execute(projectID=projectID, path=projectPath, assembly='Y', project_name=publicName)

        linfo = utils.createTableObject('load_info', session)
        insert = linfo.insert()
        insert.execute(project_id=projectID, data='init')

    else:
        sys.stderr.write("This project already exists, no need to create tables - exiting!")
        return



def loadHomology(projectID, linkDict, alnParser, session, remove = False):


    if debug:
      print "a few members of linkDict";
      num = 0
      for key, value in linkDict.iteritems():
        print key, value
        num = num +1
        if num > 5: 
          break

    cutoff = 1e-04
    #open the load files for writing
    tmpFilePrefix = os.getcwd() + "/"
    pid = os.getpid()    
    homologyLoad = tmpFilePrefix + "_homology_load.txt." + str(pid)
    annotationLoad = tmpFilePrefix + "_" + str(alnParser.database) + "_annotation_update_load.txt." + str(pid)

    hl = open(homologyLoad, 'w')
    al = open(annotationLoad, 'w')

    firstAln = True
    firstAnnot = True
    seen = {}
    noEntries = 0
    for query, annotation in alnParser.annotations.iteritems():
        query = query.strip()
        query = query.split(" ")[0]
        sbid = linkDict.get(query)
        n = 0
        if(len(annotation) < 5):
            maxHits = len(annotation)
        else:
            maxHits = 5
        while(n < maxHits):
            if(annotation[n].significance > cutoff):
                break
            #check if annotation for this hit sequence already exists
            refParts = annotation[n].reference_id.split("|")
            hitId = refParts[1]
            hitText = " ".join(refParts[2:])
            #if it doesn't exist, insert the annotation into the annotation load file
            t = utils.createTableObject('annotation_db', session)
            results = t.select(and_(t.c.annotation_id==hitId, t.c.source==alnParser.database)).execute()
            if(not results.fetchone()):#s.getData("SELECT annotation_id FROM annotation_db WHERE annotation_id =%s and source=%s", (hitId, args.annotDb))):
                if(not hitId in seen):
                    seen[hitId] = True
                    if(firstAnnot):
                        firstAnnot = False
                        al.write("\t".join([str(hitId), hitText, str(alnParser.database)]))
                    else:
                        al.write("\n" + "\t".join([str(hitId), hitText, str(alnParser.database)]))
            #check if this is the first homology entry or not, either way write to homology load file
            if(firstAln):
                hl.write("\t" + "\t".join([str(projectID), str(sbid),str(hitId), str(annotation[n].significance), str(alnParser.database)]))
                firstAln = False
            else:
                hl.write("\n\t" + "\t".join([str(projectID), str(sbid),str(hitId), str(annotation[n].significance), str(alnParser.database)]))
            n += 1
        noEntries += 1
    hl.close()
    al.close()

    #load the annotation files
    session.activateConnection()
    if settings.machine == 'oem':
      with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        session.conn.execute("LOAD DATA INFILE '" + homologyLoad + "' INTO TABLE homology")
    else:
      with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        session.conn.execute("LOAD DATA LOCAL INFILE '" + homologyLoad + "' INTO TABLE homology")
    if settings.machine == 'oem':
      session.conn.execute("LOAD DATA INFILE '" + annotationLoad + "' REPLACE INTO TABLE annotation_db")
    else:
      session.conn.execute("LOAD DATA LOCAL INFILE '" + annotationLoad + "' REPLACE INTO TABLE annotation_db")

    if not debug:
      os.remove(homologyLoad)
      os.remove(annotationLoad)

    #store sorted homology

    t1 = ""
    t2 = ""
    if settings.machine == 'oem':
      t1 = "/home/oem/python_software/autonomics/scripts/ss.done1"
      t2 = "/home/oem/python_software/autonomics/scripts/ss.done2"
    else:
      t1 = "/home/pwilliams/python_software/autonomics/scripts/ss.done1"
      t2 = "/home/pwilliams/python_software/autonomics/scripts/ss.done2"

    subprocess.Popen("perl \"" + perlPath + "storeSortedHomology.pl\" " + str(projectID) + " 1 " + str(debug) , shell=True).wait()
    proc = subprocess.Popen("perl \"" + perlPath + "storeSortedHomology.pl\" " + str(projectID) + " 2 " + str(debug), shell=True)

    proc.wait()
    
    if(not os.path.exists(t1)):
        print t1, " does not exist"
        print "storeSortedHomology.pl failed for sort 1"
        sys.exit()        

    if(not os.path.exists(t2)):
        print t2, " does not exist"
        print "storeSortedHomology.pl failed for sort 2"
        sys.exit()        

    os.remove(t1)
    os.remove(t2)

    session.conn.execute("DELETE FROM homology")

def loadHomologyAlignments(projectID, linkDict, alnReader, session, v = False):
    d = os.getcwd()
    pid = os.getpid()
    alignmentLoad = d + "/" + str(alnReader.database) + "_alignment_load.txt." + str(pid)
    alignments = open(alignmentLoad, 'w')
    reportStatus("Writing load files.\n", v)
    firstAln = True
    noEntries = 0
    for query, annotations in alnReader.annotations.iteritems():
        query = query.split(" ")[0]
        try:
            sbid = linkDict[query.strip()]
        except KeyError as e:
            print("Could not find sb_id for sequence identifier: " +
                  query.strip() + ", skipping.")
            continue
        
        n = 0
        if(len(annotations) < 5): maxHits = len(annotations)
        else: maxHits = 5
        while(n < maxHits):
            #check if annotation for this hit sequence already exists
            refParts = annotations[n].reference_id.split("|")
            hitId = refParts[1]
            #check if this is the first homology entry or not, either way write to alignment load file
            if(firstAln):
                alignments.write("\t".join([str(sbid), str(hitId), str(alnReader.database), annotations[n].aln_text]))
                firstAln = False
            else:
                alignments.write("\n" + "\t".join([str(sbid), str(hitId), str(alnReader.database), annotations[n].aln_text]))

            n += 1
        noEntries += 1
    alignments.close()
    reportStatus("\n", v)
    #store annotation alignments
    reportStatus("Loading mysql table  annotation alignments...\n", v)
    session.activateConnection()

    if settings.machine == 'oem':
      session.conn.execute("LOAD DATA INFILE '" + alignmentLoad + "' REPLACE INTO TABLE annotation_alignments")
    else:
      session.conn.execute("LOAD DATA LOCAL INFILE '" + alignmentLoad + "' REPLACE INTO TABLE annotation_alignments")

    if not debug:    os.remove(alignmentLoad)


def loadQuantification(projectID, projectName, abundanceFile, abundanceCol, s, v=False):
    if(abundanceFile == None):
        sys.stderr.write("No abundance file specified, using default.\n")
        pid = os.getpid()    
        abundanceFile = baseDir + projectName + "/" + projectName + "_quantification.txt." + str(pid)

    reportStatus("Loading quantification data...\n", v)
    seq_table = utils.createTableObject(str(projectID) + "_sequences", s)
    abundances = {}
    abundanceCol -= 1
    fh = open(abundanceFile, 'r')
    for l in fh:
        l = l.rstrip("\n")
        el = l.split("\t")
        abundances[el[0]] = el[abundanceCol]
        
    results = seq_table.select().execute()
    for row in results.fetchall():
        if(row.description in abundances):
            u = seq_table.update().where(seq_table.c.sb_id==row.sb_id).values(abundance=abundances[row.description])
            u.execute()
        
    proj_dir = utils.createTableObject("project_directory", s)
    u = proj_dir.update().where(proj_dir.c.projectID==projectID).values(has_abundance='Y')

    linfo = utils.createTableObject('load_info', session)
    insert = linfo.insert()
    insert.execute(project_id=projectID, data='quant')



def loadKEGG(projectID, linkDict, keggFile, session, v=False):
    if(v): reportStatus("Parsing KEGG file...\n")
    lines = parseCSV(keggFile, "\t")

    if(v):
        reportStatus("Collapsing KEGG file...\n")
    collapsedLines = collapseAssociationFile(lines, [0, 2, 3], 6, "smallest")
    if(v): reportStatus("Removing old KEGG entries for this project...\n")
    if(session.conn == None):
        session.activateConnection()
    session.conn.execute("DELETE FROM kegg_annotations WHERE project_id='" + str(projectID) + "'")
    datak = 'kegg'
    session.conn.execute("DELETE FROM load_info WHERE data='" + datak + "' and project_id='" + str(projectID) + "'")

    kegg = utils.createTableObject('kegg_annotations', session)

    if(v): reportStatus("Inserting KEGG annotations...\n")
    for key, elements in collapsedLines.items():
#        sbID = linkDescToID(elements[0].replace(" ",""), linkDict)
        sbID = linkDescToID(elements[0].strip(), linkDict)
        insert = kegg.insert()
        try:
            if debug and debug_kegg:
              print "insert.execute(project_id=",projectID," sb_id = ",sbID, " kegg_id= ", str(elements[2]).strip(), " path_id= ",str(elements[3]).strip(), " kegg_pathway_description= ",str(elements[4]).strip(), " annotating_acc= ", elements[1]," evalue= ",elements[6]
              print "elem[0] = ", str(elements[0])

            insert.execute(project_id=projectID, sb_id = sbID, kegg_id=str(elements[2]).strip(), path_id=str(elements[3]).strip(), kegg_pathway_description=str(elements[4]).strip(), annotating_acc=elements[1], evalue=elements[6])
        #catch sqlalchemy exception
        except sqlalchemy.exc.OperationalError:
            #do not insert this value, print warning
            print("Error inserting KEGG record for: " + str(elements[0]))


    linfo = utils.createTableObject('load_info', session)
    insert = linfo.insert()
    insert.execute(project_id=projectID, data='kegg')


def loadGO(projectID, linkDict, goFile, session, v = False):

    tmpFilePrefix = os.getcwd() + "/"
    pid = os.getpid()    
    loadFile = tmpFilePrefix + "go_load.txt." + str(pid)
    curId = ""
    reportStatus("Parsing GeneOntology CSV file...\n", v)
    f = parseCSV(goFile, "\t")

    gLoad = open(loadFile, 'w')
    reportStatus("Creating load file...\n", v)
    records = 0
    
    for elements in f:
        elen = len(elements)
        curId = linkDescToID(elements[0], linkDict)
        goId = elements[2].split(":")[1]
        if debug:
           if records < 10:
             print "==============="
             print str(goId), "\t", str(curId), "\t", elements[1], "\t", elements[elen-1]
             print "==============="
        gLoad.write("".join(["\t".join([str(goId), str(curId), elements[1], elements[elen-1]]), "\n"]))
        records += 1
    gLoad.close()

    session.activateConnection()
    print "Loading ", records, " GO records into sql database...\n"
    if settings.machine == 'oem':
      session.conn.execute("LOAD DATA INFILE '" + loadFile + "' REPLACE INTO TABLE go_annotation_new")
    else:
      session.conn.execute("LOAD DATA LOCAL INFILE '" + loadFile + "' REPLACE INTO TABLE go_annotation_new")

    if not debug:    
      print "os.remove(loadFile)" 
      os.remove(loadFile)
    reportStatus("loadGO done.\n", v)


def loadGOCategories(projectID, projectName, catFile, session, v):
    #open the go category file
    fh = open(catFile, 'r')
    tmpFilePrefix = os.getcwd() + "/"
    pid = os.getpid()    
    tmpFile = tmpFilePrefix + "gocat_load.txt." + str(pid)

    tmp = open(tmpFile, 'w')
    for line in fh:
        tmp.write(line.replace(projectName, str(projectID)))
    
    fh.close()
    tmp.close()
    reportStatus("Loading mysql table go_categories...\n", v)
    session.activateConnection()

    if settings.machine == 'oem':
      session.conn.execute("LOAD DATA INFILE '" + tmpFile + "' REPLACE INTO TABLE go_categories" )
    else:
      session.conn.execute("LOAD DATA LOCAL INFILE '" + tmpFile + "' REPLACE INTO TABLE go_categories" )

    if not debug:    os.remove(tmpFile)


def loadPfam(projectID, pfamFile, seqType, linkDict, session, v=True):
    pfamFile = rewritePfam(pfamFile, seqType, linkDict)
    tmpFilePrefix = os.getcwd() + "/"
    pid = os.getpid()    
    tmpFile = tmpFilePrefix + "pfam_load_tmp.txt." + str(pid)
    flattened = pfamFile + "_flattened.txt"
    reportStatus("Flattening Pfam annotation file...\n", v)
    proc = subprocess.Popen("python \"" + settings.python_scripts_path + "filetools.py\" --fields seq_id:1 alignment_start:2 alignment_end:3 hmm_acc:6 E-value:13:float --key-col 1 2 3 6 --filter E-value:lt:1e-04 --flatten \"" + pfamFile + "\"", shell=True)
    proc.wait()

    reportStatus("Writing load file...\n", v)
    fh = open(flattened, 'r')
    tmpOut = open(tmpFile, 'w')
    for line in fh:
        line = line.rstrip("\n")
        el = line.split("\t")
        tmpOut.write("".join(["\t".join([str(projectID), str(el[0]), el[3].split(".")[0], str(el[4]), str(el[1]), str(el[2])]), "\n"]))
        
    fh.close()
    tmpOut.close()
    if not debug:    os.remove(flattened)

    
    reportStatus("Loading pfam annotations...\n", v)
    #load the data into the pfam_annotations table
    if(session.conn == None):
        session.activateConnection()

    if settings.machine == 'oem':
      session.conn.execute("LOAD DATA INFILE '" + tmpFile + "' REPLACE INTO TABLE pfam_annotations")
    else:
      session.conn.execute("LOAD DATA LOCAL INFILE '" + tmpFile + "' REPLACE INTO TABLE pfam_annotations")

    reportStatus("Loading pfam domain counts...\n", v)
    #get category information
    pfam = utils.createTableObject('pfam_annotations', session)
    categories = utils.createTableObject('pfam_domain_counts', session)
    pfama = utils.createTableObject('pfama', session)
    #count the occurences of each pfam_acc
    results = session.conn.execute(select([func.count("*").label('pfamCatCount'), pfam.c.pfamA_acc, pfama.c.description], and_(pfam.c.pfamA_acc==pfama.c.pfamA_acc, pfam.c.project_id==projectID)).group_by(pfam.c.pfamA_acc))
    rows = results.fetchall()
    for row in rows:
        #generate a new insert statement
        i = categories.insert()
        i.execute(project_id=projectID, acc=row.pfamA_acc, description = row.description, counts = row.pfamCatCount)
    if not debug:    os.remove(tmpFile)
    if not debug:    os.remove(pfamFile)


def loadSequences(projectName, seqType, session, seqFile=None, v = False, sbStart=None):
    #check if this project has abundances yet
    if(session.conn == None):
        session.activateConnection()

    uploadDir = baseDir + projectName + "/"
    if(seqFile is None):
        seqFile = uploadDir + projectName + "_project.fasta"

    projectID = getProjectID(projectName, session, v)

    #get the next file identifier
    catalog = utils.createTableObject('sb_catalog', session)

    results = session.conn.execute(select([func.max(catalog.c.fileID).label("maxFile")]))
    row = results.fetchone()
    try:
        fileID = int(row['maxFile']) + 1
    except:
        fileID = 1

    #get the starting sb#
    if(sbStart is None):
        results = session.conn.execute(select([func.max(catalog.c.end).label("max")]))
        row = results.fetchone()
        
        try:
            curId = int(row['max']) + 1
        except:
            curId = 1
        sbStart = curId
    else:
        curId = int(sbStart)

    #open the database file for writing
    directory = utils.createTableObject('project_directory', session)
    dbPath = settings.NEUROBASE_FASTADB_PATH + str(projectID)
    if(not os.path.exists(dbPath)):
        os.mkdir(dbPath)
    dbFile = open(dbPath + "/" + seqType + "DatabaseFile.fas", 'w')

    #open a temporary load file
    loadFilePath = uploadDir + "tempSeqLoad" + str(time.strftime("%b%d%Y%H%M%S")) + ".txt"
    loadFile = open(loadFilePath, 'w')
    numSeqs = 0
    
    seqhandle = open(seqFile, 'r')
    
    #write the sequences to the load file
    print "writing seqs to load file: ", loadFilePath
    today = date.today()
    for rec in SeqIO.parse(seqhandle, "fasta"):
        if(seqType == "NT"):
            loadFile.write("\t".join([str(curId), str(fileID), str(rec.seq), "", str(rec.id), seqType, str(len(str(rec.seq))),str(today),"1"]))
            loadFile.write("\n")
        else:
            loadFile.write("\t".join([str(curId), str(fileID), "", str(rec.seq), str(rec.id), seqType, str(len(str(rec.seq))),str(today),"1"]))
            loadFile.write("\n")
            
        rec.id = "sb|" + str(curId) + "|"
        dbFile.write(">" + rec.id + "\n" + str(rec.seq) + "\n")
        curId += 1
        numSeqs += 1
    print "DONE writing seqs to load file"

    seqhandle.close()
    loadFile.close()
    dbFile.close()

    #load the sequences to the table
    retries = 5
    while(retries > 0):
        try:
            if settings.machine == 'oem':
              session.conn.execute("TRUNCATE TABLE " + str(projectID) + "_sequences")
              session.conn.execute("LOAD DATA INFILE '" + loadFilePath + "' INTO TABLE " + str(projectID) + "_sequences")
            else:
              session.conn.execute("delete from " + str(projectID) + "_sequences")
              session.conn.execute("LOAD DATA LOCAL INFILE '" + loadFilePath + "' INTO TABLE " + str(projectID) + "_sequences")
            break
        
        except sqlalchemy.exc.OperationalError as e:
            time.sleep(60)
            print("MySQL Server went away, reconnecting and trying again.")
            session = netutils.make_db_session("moroz_lab")
    #remove the temporary file
    if not debug:    os.remove(loadFilePath)

    print "DONE loading seqs to nn_sequences"

    print "formating DB for blast"

    #format the database file for BLAST                                       
    db_file = dbPath + "/" + seqType + "DatabaseFile.fas"
    #first try with formatdb                                                  
    if(seqType == "NT"): typeFlg = "-p F"
    else: typeFlg = "-p T"
    p = subprocess.Popen('formatdb -i ' + db_file + ' ' + typeFlg,shell=True)
    p.wait()
    
    if(p.returncode != 0):
        #try with makeblastdb
        print "formating with makeblastdb instead"
        if(seqType == "NT"):
            typeFlg = "-dbtype nucl"
        else:
            typeFlg = "-dbtype prot"

        p = subprocess.Popen('makeblastdb -in ' + db_file + " " + typeFlg, shell=True)
        p.wait()
        utility.die_on_error(p.returncode)

    #insert records into the sb_catalog
    i = catalog.insert()
    i.execute(begin=sbStart, end=curId - 1, fileID=fileID, projectID=projectID)

    #update the number of sequences in the project directory
    if(seqType == "NT"): 
       session.conn.execute(directory.update().where(directory.c.projectID == projectID).values(num_NT_seqs=numSeqs))
    else: 
       session.conn.execute(directory.update().where(directory.c.projectID == projectID).values(num_AA_seqs=numSeqs))

    linfo = utils.createTableObject('load_info', session)
    insert = linfo.insert()
    insert.execute(project_id=projectID, data='seqs')

    #close the session


def parseCSV(file, d = ",", q = "\""):
    f = open(file, 'r')
    parsed = []
    for l in f:
        l = l.strip()
        token = ""
        delim = d
        quote = q
        quoteOpen = False
        elements = []
        for c in l:
            if(c == delim):
                if(not quoteOpen):
                    elements.append(token)
                    token = ""
            elif(c == quote):
                quoteOpen = not quoteOpen
            else:
                token += c
        #append the last element if things are cool
        if(not quoteOpen):
            elements.append(token)
        parsed.append(elements)
    return parsed
    f.close()


def repairDescription(description, attempt):

    if(attempt == 1):
        return description.split(" ")[0]
    else:
        return None

def rewritePfam(pfamFile, seqType, linkDict):
    f = open(pfamFile, 'r')
    new = open(pfamFile + "_linked.out", 'w')
    spaceSplit = re.compile("\s+")
    translateSub = re.compile("_-*\d$")
    for l in f:
        l = l.strip()
        if(not l.startswith("#")):
            elements = spaceSplit.split(l)
            if(len(elements) > 1):
                if(seqType == 'NT'):
                    #this means the query sequences were translated and have the frame hanging off the end
                    elements[0] = translateSub.sub("", elements[0])
                elements[0] = linkDict[elements[0]]
                new.write("\t".join(map(str, elements)) + "\n")
    f.close()
    new.close()
    return pfamFile + "_linked.out"


def reportStatus(txt, v = True):
    if(v): sys.stdout.write(txt); sys.stdout.flush()


def setDBs(db, args):
    if(db == 'nr'):
        args.annotDb=1
    else:
        args.annotDb=2

# not used:
def storeSortedHomology(projectID, sortType):
    subprocess.Popen("perl \"" + perlPath + "storeSortedHomology.pl\" " + str(projectID) + " " + str(sortType), shell=True).wait()

def main():
    print sys.argv
    parser = argparse.ArgumentParser(description = "This script contains various methods of loading data into NeuroBase.")
    parser.add_argument('-a', '--append', dest='append', default=False, const=True, action='store_const')
    #parser.add_argument('-w', '--write', dest='write', default=False, const=True, action='store_const')
    #parser.add_argument('-c', '--create-tables', dest='createTables', default=False, const=True, action='store_const')
    #parser.add_argument("--annotDb", dest="annotDb", default="nr")
    parser.add_argument("--alignments", dest="alignments", default=False, const=True, action='store_const')
    parser.add_argument("--annotation", dest="annotation", default=False, const=True, action='store_const')
    parser.add_argument("--data-dir", dest='data_dir', default=None)
    parser.add_argument("--init-project", dest="initProject", default=False, const=True, action='store_const')
    parser.add_argument("--nr", dest='nr', default=False, const=True, action='store_const')
    parser.add_argument("--swissprot", dest='swissprot', default=False, const=True, action='store_const')
    parser.add_argument("--alignment-format", dest='alnFmt', default='hpc-blast')
    parser.add_argument("--delete-homology", dest="deleteH", default=False, const=True, action='store_const')
    parser.add_argument('--GO', dest='GO', default=False, const=True, action='store_const')
    parser.add_argument("--GOCategories", dest='GOCategories', default=False, const=True, action='store_const')
    parser.add_argument('--go-file', dest='goFile', default=None)
    parser.add_argument("--KEGG", dest='KEGG', default=False, const=True, action='store_const')
    parser.add_argument("--kegg-file", dest='keggFile', default=None)
    parser.add_argument('--pfam', dest='pfam', default=False, const=True, action='store_const')
    parser.add_argument('--pfam-file', dest='pfamFile', default=None)
    parser.add_argument("--quantification", dest='quant', default=False, const=True, action='store_const')
    parser.add_argument('--quantification-file', dest='quantFile', default=None)
    parser.add_argument('--abundance-column', dest='abundanceCol', default=2)
    parser.add_argument("--seqs", dest="loadSeqs", default=False, const=True, action="store_const")
    parser.add_argument("--seq-type", dest='seqType', default='NT')
    parser.add_argument("--sbid-start-value", dest="sbStart", default=None)
    parser.add_argument("--best-annotation", dest='bestAnnotation', default=False, const=True, action='store_const')
    parser.add_argument('-p', '--project-name', dest='pn', required=True, help='Name of the project you\'re working with.')
    parser.add_argument('--public-name', dest='publicName', default=None)
    #parser.add_argument('-f', '--data-file', dest='data', required=True, help='File containing the data to load to NeuroBase.')
    parser.add_argument('-s', '--seq-file', dest='seqFile', default=None)
    parser.add_argument('--user', dest='user', required=True, help='Username for mysql database containing NeuroBase data')
    parser.add_argument('--password', dest='passwd', required=True, help='Password for mysql database containing NeuroBase data')
    parser.add_argument('-v', dest='verbose', default=False, const=True, action='store_const')
    parser.add_argument('--debug', dest='debug', default=False, const=True, action='store_const')
    parser.add_argument('--no-quant', dest='no_quant', default=False, const=True, action='store_const')
    args = parser.parse_args()

    global debug
    global debug_kegg
    global quant

    v = args.verbose
    debug = args.debug
    no_quant = args.no_quant

    if debug: debug = 1
    else: debug = 0

    if no_quant: quant = 0
    else: quant = 1

    print "debug = ", debug
    print "quant = ", quant
    debug_kegg = 0
    print "debug_kegg = ", debug_kegg

    if (quant):
      print "running quant"
    else:      print "NOT running quant"

    reportStatus("Starting a new MySQLseq session with NeuroBase\n", v)
    s = utils.DBSession(args.user, args.passwd, 'localhost')
    session = s

    if(not args.data_dir is None):
        baseDir = args.data_dir
    else:
        baseDir = settings.NEUROBASE_DATA_PATH

    if(not baseDir.endswith("/")):
            baseDir += "/"
            
    baseName = baseDir + args.pn + "/" + args.pn
    if(args.publicName is None): 
        args.publicName = args.pn
    
    projectID = getProjectID(args.pn, s, v)

    if(args.annotation):

        projectID = getProjectID(args.pn, s, v)
        #create linkage of assembly ID to neurobase ID
        reportStatus("Linking sequence identifiers in FASTA file to NeuroBase ID\n", v)
        linkDict = utils.link_dbid_to_fastaid(projectID, s)

        p = None

        name = baseName + "_blast_swissprot.txt"
        if(args.swissprot):
            reportStatus("Loading swissprot homology data...\n",)
            p = Reader(baseName + "_blast_swissprot.txt", args.alnFmt, 2, v)
        elif(args.nr):
            name = baseName + "_blast_nr.txt"
            reportStatus("Loading nr homology data.\n", v)
            p = Reader(baseName + "_blast_nr.txt", args.alnFmt, 1, v)
        p.read()
        loadHomology(projectID, linkDict, p, s, args.deleteH)
        if(args.swissprot):
          linfo = utils.createTableObject('load_info', session)
          insert = linfo.insert()
          insert.execute(project_id=projectID, data='swiss')
        elif(args.nr):
          linfo = utils.createTableObject('load_info', session)
          insert = linfo.insert()
          insert.execute(project_id=projectID, data='nr')
        assignBestAnnotation(projectID, s, v)

    elif(args.alignments):
        projectID = getProjectID(args.pn, s, v)
        reportStatus("Linking sequence identifiers in FASTA file to Neurobase ID\n", v)
        linkDict = utils.link_dbid_to_fastaid(projectID, s)
        if(args.swissprot):
            r = Reader(baseName + "_blast_swissprot.txt", args.alnFmt, 2, v)

        elif(args.nr):
            r = Reader(baseName + "_blast_nr.txt", args.alnFmt, 1, v)

        r.read()
        reportStatus("Loading annotation alignments.\n", v)
        loadHomologyAlignments(projectID, linkDict, r, s, v)
        linfo = utils.createTableObject('load_info', session)
        insert = linfo.insert()
        insert.execute(project_id=projectID, data='anno_align')

    elif(args.bestAnnotation):
        assignBestAnnotation(projectID, s, v)

    elif(args.initProject):
        initProject(args.publicName, s)

    elif(args.quant):
        projectID = getProjectID(args.pn, s, v)
        if(args.quantFile is None):
           args.quantFile = baseName + "_quantification.txt"
        loadQuantification(projectID, args.pn, args.quantFile, args.abundanceCol, s, v)

    elif(args.KEGG):
        projectID = getProjectID(args.pn, s, v)
        #create linkage of assembly ID to neurobase ID
        linkDict = utils.link_dbid_to_fastaid(projectID, s)
        if(args.keggFile is None):
            args.keggFile = baseName + "_KEGG.txt"
        loadKEGG(projectID, linkDict, args.keggFile, s, v)

    elif(args.GO):
        projectID = getProjectID(args.pn, s, v)
        #create linkage of assembly ID to neurobase ID
        linkDict = utils.link_dbid_to_fastaid(projectID, s)
        if(args.goFile is None):
            args.goFile = baseName + "_GO.txt"
        loadGO(projectID, linkDict, args.goFile, s, v)
        linfo = utils.createTableObject('load_info', session)
        insert = linfo.insert()
        insert.execute(project_id=projectID, data='go')
        
    elif(args.GOCategories):
        projectID = getProjectID(args.pn, s, v)
        catFile = baseName + "_gocats.txt"
        loadGOCategories(projectID, args.pn, catFile, s, v)
        linfo = utils.createTableObject('load_info', session)
        insert = linfo.insert()
        insert.execute(project_id=projectID, data='gocat')

    elif(args.pfam):
        projectID = getProjectID(args.pn, s, v)
        linkDict = utils.link_dbid_to_fastaid(projectID, s)
        if(args.pfamFile is None):
            args.pfamFile = baseName + "_pfam.txt"
            
        loadPfam(projectID, args.pfamFile, args.seqType, linkDict, s)
        linfo = utils.createTableObject('load_info', session)
        insert = linfo.insert()
        insert.execute(project_id=projectID, data='pfam')

    elif(args.loadSeqs):
        reportStatus("Loading project sequences...\n")
        loadSequences(args.pn, args.seqType, s, args.seqFile, v, sbStart=args.sbStart)
        linfo = utils.createTableObject('load_info', session)
        insert = linfo.insert()
        insert.execute(project_id=projectID, data='seqs')

    else:
        allPipe(args, s, v)

    print "Completed upload!!\n"

if __name__ == '__main__':
    main()
