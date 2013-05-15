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
from zeroclick.file_io import Record
from zeroclick import settings, netutils, utility
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

perlPath="/home/oem/PerlScripts/Database/"
pythonPath = "/home/oem/moroz-python/"
baseDir = settings.NEUROBASE_DATA_PATH

def allPipe(args, s, v):

    s.activateConnection()
    
    #ftp directory
    baseDir = settings.NEUROBASE_DATA_PATH
    if(not args.data_dir is None):
        baseDir = args.data_dir
        
    baseName = baseDir + args.pn + "/" + args.pn

    reportStatus("Initiating project...", v)
    initProject(args.publicName, s)
    reportStatus("done.\n", v)

    reportStatus("Loading project sequences...", v)
    seqs = baseName + "_project.fasta"
    loadSequences(args.pn, args.seqType, s, seqs, v, sbStart=args.sbStart)
    reportStatus("done.\n", v)

    projectID = getProjectID(args.pn, s, True)

    reportStatus("Linking sequence identifiers in FASTA file to NeuroBase ID...", v)
    linkDict = utils.link_dbid_to_fastaid(projectID, s)
    reportStatus("done.\n", v)

    reportStatus("Loading quantification data...", v)
    if(args.quantFile is None):
        args.quantFile = baseName + "_quantification.txt"
    loadQuantification(projectID, args.pn, args.quantFile, args.abundanceCol, s, v)
    reportStatus("done.\n", v)
    
    reportStatus("Loading swissprot homology data...", v)
    p = Reader(baseName + "_blast_swissprot.txt", args.alnFmt, 2, v)
    p.read()
    loadHomology(projectID, linkDict, p, s, v, args.deleteH)
    reportStatus("done.\n", v)
    
    reportStatus("Storing annotation alignments...", v)
    loadHomologyAlignments(projectID, linkDict, p, s, v)
    reportStatus("done.\n", v)

    reportStatus("Assigning best annotations to sequences..", v)
    assignBestAnnotation(projectID, s, v)
    reportStatus("done\n", v)

    reportStatus("Loading GO terms...", v)
    if(args.goFile is None):
        args.goFile = baseName + "_GO.txt"
    loadGO(projectID, linkDict, args.goFile, s, v)
    reportStatus("done.\n", v)

    reportStatus("Loading GO categories..", v)
    catFile = baseName + "_gocats.txt"
    loadGOCategories(projectID, args.pn, catFile, s, v)
    reportStatus(".done\n", v)
    
    reportStatus("Loading KEGG pathways...", v)
    if(args.keggFile is None):
        args.keggFile = baseName + "_KEGG.txt"
    loadKEGG(projectID, linkDict, args.keggFile, s, v)
    reportStatus("done.\n", v)
    
    reportStatus("Loading pfam annotations...", v)
    if(args.pfamFile is None):
            args.pfamFile = baseName + "_pfam.txt"
    loadPfam(projectID, args.pfamFile, args.seqType, linkDict, s)
    reportStatus("done.\n", v)


def assignBestAnnotation(projectID, session, v):
    if(v):
        reportStatus("Assigning best annotation.\n", v)

    session.activateConnection()
    sorted_homology = netutils.get_table_object("sorted_homology", session)
    annotation_db = netutils.get_table_object('annotation_db', session)
    seq_table = netutils.get_table_object(str(projectID) + "_sequences", session)
    best_annot = netutils.get_table_object('best_annotations', session)
    
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
        #fields = [annotation_db.c.text, sorted_homology.c.evalue]
        #q = select(fields, and_(sorted_homology.c.sort_id==1, 
#                                        sorted_homology.c.project_id==projectID, 
#                                        sorted_homology.c.sb_id==sb_id, 
#                                        sorted_homology.c.annotation_id==annotation_db.c.annotation_id)).order_by(sorted_homology.c.ranking).offset(0).limit(1)
        #res = q.execute()
        #annot_row = res.fetchone()
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
        
    '''
    session.activateConnection()
    skipped = 0
    for seqID, annotations in alnReader.annotations.iteritems():
        seqID = seqID.split(" ")[0]
        seqID = seqID.strip()
        if(not seqID in linkDict):
            skipped += 1
            continue
        sb_id = linkDict[seqID]
        best = annotations[0]
        values = "(\"" + str(projectID) + "\", \"" + str(sb_id) + "\", \"" + best.reference_id + "\", \"" + str(best.significance) + "\")"
        update = "ON DUPLICATE KEY UPDATE annot = (IF(eval > " + str(best.significance) + ", \"" + str(best.reference_id) + "\", annot)), eval = (IF(eval > " + str(best.significance) + ", " + str(best.significance) + ", eval))"
        query = "INSERT INTO best_annotations (project_id, sb_id, annot, eval) values " + values + " " + update
        session.conn.execute(query)
    print(str(skipped) + " sequences skipped while assigning best annotations.")
    '''

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
        sys.stderr.write("Couldn't find the project name you specified, please check your submission if the project is already created.")
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

    else:
        sys.stderr.write("This project already exists, no need to create tables - exiting!")
        return



def loadHomology(projectID, linkDict, alnParser, session,
                 cutoff=1e-04, v = False, remove = False):
    #check if this project has abundance information loaded
    if(session.conn == None):
        session.activateConnection()
    directory = utils.createTableObject("project_directory", session)
    #result = session.conn.execute(select([directory.c.has_abundance], directory.c.projectID==projectID))
    #if(result.fetchone()["has_abundance"] != "Y"):
    #    sys.stderr.write("This project does not yet have sequence abundance associated with it, which is required for storing homology data. Please load the abundance data and try again.\n")
    #    sys.exit()
    #open the load files for writing
    tmpFilePrefix = os.getcwd() + "/"
    homologyLoad = tmpFilePrefix + "_homology_load.txt"
    annotationLoad = tmpFilePrefix + "_" + str(alnParser.database) + "_annotation_update_load.txt"
    hl = open(homologyLoad, 'w')
    al = open(annotationLoad, 'w')
    #remove existing annotations for this project to avoid duplications
    '''if(remove):
        removeExistingHomology(projectID, session)
        ^^ NOT IMPLEMENTED
    '''
    firstAln = True
    firstAnnot = True
    seen = {}
    reportStatus("Writing annotation and homology entries to load files\n")
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
        reportStatus("%d lines written to load files.\r" % noEntries)
    reportStatus("\n")
    hl.close()
    al.close()

    #load the annotation files
    session.activateConnection()
    reportStatus("Importing load files to database...\n", v)
    reportStatus("homology data...\n", v)
    session.conn.execute("LOAD DATA INFILE '" + homologyLoad + "' INTO TABLE homology")
    reportStatus("annotation data...\n", v)
    session.conn.execute("LOAD DATA INFILE '" + annotationLoad + "' REPLACE INTO TABLE annotation_db")
    os.remove(homologyLoad)
    os.remove(annotationLoad)

    #store sorted homology
    reportStatus("sorting homology entries...\n", v)
    subprocess.Popen("perl \"" + perlPath + "storeSortedHomology.pl\" " + str(projectID) + " 1", shell=True).wait()
    proc = subprocess.Popen("perl \"" + perlPath + "storeSortedHomology.pl\" " + str(projectID) + " 2", shell=True)
    proc.wait()
    reportStatus("done.\n", v)

    #empty the homology table
    session.conn.execute("TRUNCATE TABLE homology")
    #if(v): reportStatus("done.\n")

def loadHomologyAlignments(projectID, linkDict, alnReader, session, v = False):
    d = os.getcwd()
    alignmentLoad = d + "/" + str(alnReader.database) + "_alignment_load.txt"
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
        reportStatus("%d queries written to load files.\r" % noEntries, v)
    alignments.close()
    reportStatus("\n", v)
    #store annotation alignments
    reportStatus("Loading annotation alignments...", v)
    session.activateConnection()
    session.conn.execute("LOAD DATA INFILE '" + alignmentLoad + "' REPLACE INTO TABLE annotation_alignments")
    os.remove(alignmentLoad)
    reportStatus("done.\n", v)


def loadQuantification(projectID, projectName, abundanceFile, abundanceCol, s, v=False):
    if(abundanceFile == None):
        sys.stderr.write("No abundance file specified, using default.\n")
        abundanceFile = baseDir + projectName + "/" + projectName + "_quantification.txt"

    reportStatus("Loading quantification data...", v)
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

    reportStatus("done.\n", v)    


def loadKEGG(projectID, linkDict, keggFile, session, v=False):
    if(v): reportStatus("Parsing KEGG file...")
    lines = parseCSV(keggFile, "\t")
    if(v): reportStatus("done\n")
    if(v):
        reportStatus("Collapsing KEGG file...")
    collapsedLines = collapseAssociationFile(lines, [0, 2, 3], 6, "smallest")
    if(v):
        reportStatus("done.\n")
    if(v): reportStatus("Removing old KEGG entries for this project...")
    if(session.conn == None):
        session.activateConnection()
    session.conn.execute("DELETE FROM kegg_annotations WHERE project_id='" + str(projectID) + "'")
    kegg = utils.createTableObject('kegg_annotations', session)
    if(v): reportStatus("done\n")
    if(v): reportStatus("Inserting KEGG annotations...")
    for key, elements in collapsedLines.items():
        sbID = linkDescToID(elements[0], linkDict)
        insert = kegg.insert()
        try:
            insert.execute(project_id=projectID, sb_id = sbID, kegg_id=str(elements[2]).strip(), path_id=str(elements[3]).strip(), kegg_pathway_description=str(elements[4]).strip(), annotating_acc=elements[1], evalue=elements[6])
        #catch sqlalchemy exception
        except sqlalchemy.exc.OperationalError:
            #do not insert this value, print warning
            print("Error inserting KEGG record for: " + str(elements[0]))

    if(v): reportStatus("done\n")


def loadGO(projectID, linkDict, goFile, session, v = False):
    #seqs = utils.createTableObject(str(projectID) + "_sequences", session)
    #goAnnot = utils.createTableObject('go_annotation_new', session)
    #reportStatus("Removing existing GeneOntology annotations...", v)
    #delete the existing annotations for this project
    #results = select([seqs.c.sb_id]).execute()
    #for row in results.fetchall():
    #    goAnnot.delete(goAnnot.c.sb_id==row.sb_id).execute()
    #reportStatus("done.\n", v)
    loadFile = "go_load.txt"
    curId = ""
    reportStatus("Parsing GeneOntology CSV file...", v)
    f = parseCSV(goFile, "\t")
    reportStatus("done.\n", v)
    #reportStatus("Collapsing Gene Ontology file...", v)
    #f = collapseAssociationFile(f, [0, 2], 7, "smallest")
    #reportStatus("done.\n", v)
    gLoad = open(loadFile, 'w')
    reportStatus("Creating load file...\r", v)
    records = 0
    for elements in f:
        curId = linkDescToID(elements[0], linkDict)
        goId = elements[2].split(":")[1]
        gLoad.write("".join(["\t".join([str(goId), str(curId), elements[1], elements[3]]), "\n"]))
        records += 1
        reportStatus("Creating load file..." + str(records) + " records processed\r", v)
        #i.execute(go_id=goId, sb_id=curId, project_id=projectID, go_higher=curHigher, annotator_id=elements[5], annotation_ev=elements[7], date='CURDATE()')

    gLoad.close()
    reportStatus("\nDone.\n", v)
    session.activateConnection()
    reportStatus("Loading raw GO data...", v)
    session.conn.execute("LOAD DATA INFILE 'go_load.txt' REPLACE INTO TABLE go_annotation_new")
    os.remove(loadFile)
    reportStatus("done.\n", v)


def loadGOCategories(projectID, projectName, catFile, session, v):
    #open the go category file
    fh = open(catFile, 'r')
    tmp = open("tmp.txt", 'w')
    for line in fh:
        tmp.write(line.replace(projectName, str(projectID)))
    
    fh.close()
    tmp.close()
    reportStatus("Loading GO category data...\n", v)
    session.activateConnection()
    session.conn.execute("LOAD DATA INFILE 'tmp.txt' REPLACE INTO TABLE go_categories" )
    os.remove("tmp.txt")
    reportStatus("done.\n", v)



def loadPfam(projectID, pfamFile, seqType, linkDict, session, v=True):
    pfamFile = rewritePfam(pfamFile, seqType, linkDict)
    tmp = "pfam_load_tmp.txt"
    flattened = pfamFile + "_flattened.txt"
    reportStatus("Flattening Pfam annotation file...", v)
    proc = subprocess.Popen("python \"" + pythonPath + "/zeroclick/scripts/filetools.py\" --fields seq_id:1 alignment_start:2 alignment_end:3 hmm_acc:6 E-value:13:float --key-col 1 2 3 6 --filter E-value:lt:1e-04 --flatten \"" + pfamFile + "\"", shell=True)
    proc.wait()
    reportStatus("done.\n", v)
    #os.chdir("C:\\Documents and Settings\\Administrator\\My Documents\\Dropbox\\Work\\JAva\\Bioinformatics\\build\\classes")
    #proc = subprocess.Popen("java " + javaClassPathPfam + " -Xmx1500m annotation.scripts.AnnotatePfamDomains " + str(projectID) + " " + pfamFile + " 1e-04 database > " + tmp, shell=True)
    #proc.wait()
    #write the records to tmp
    reportStatus("Writing load file...", v)
    fh = open(flattened, 'r')
    tmpOut = open(tmp, 'w')
    for line in fh:
        line = line.rstrip("\n")
        el = line.split("\t")
        tmpOut.write("".join(["\t".join([str(projectID), str(el[0]), el[3].split(".")[0], str(el[4]), str(el[1]), str(el[2])]), "\n"]))
        
    fh.close()
    tmpOut.close()
    os.remove(flattened)
    reportStatus("done.\n")
    
    reportStatus("Loading pfam annotations...", v)
    #load the data into the pfam_annotations table
    if(session.conn == None):
        session.activateConnection()
    session.conn.execute("LOAD DATA INFILE '" + tmp + "' REPLACE INTO TABLE pfam_annotations")
    reportStatus("done.\n", v)
    
    reportStatus("Loading pfam domain counts...", v)
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
    os.remove(tmp)
    os.remove(pfamFile)
    reportStatus("done.\n", v)

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
    for rec in SeqIO.parse(seqhandle, "fasta"):
        if(seqType == "NT"):
            loadFile.write("\t".join([str(curId), str(fileID), str(rec.seq), "", str(rec.id), seqType, str(len(str(rec.seq))),"CURDATE()","1"]))
            loadFile.write("\n")
        else:
            loadFile.write("\t".join([str(curId), str(fileID), "", str(rec.seq), str(rec.id), seqType, str(len(str(rec.seq))),"CURDATE()","1"]))
            loadFile.write("\n")
            
        rec.id = "sb|" + str(curId) + "|"
        dbFile.write(">" + rec.id + "\n" + str(rec.seq) + "\n")
        curId += 1
        numSeqs += 1
        reportStatus(str(numSeqs) + " written to load file.\r", v)
            

    seqhandle.close()
    loadFile.close()
    dbFile.close()

    #load the sequences to the table
    retries = 5
    while(retries > 0):
        try:
            session.conn.execute("TRUNCATE TABLE " + str(projectID) + "_sequences")
            session.conn.execute("LOAD DATA INFILE '" + loadFilePath + "' INTO TABLE " + str(projectID) + "_sequences")
            break
        
        except sqlalchemy.exc.OperationalError as e:
            time.sleep(60)
            print("MySQL Server went away, reconnecting and trying again.")
            session = netutils.make_db_session("moroz_lab")
    #remove the temporary file
    os.remove(loadFilePath)

    #format the database file for BLAST                                       
    db_file = dbPath + "/" + seqType + "DatabaseFile.fas"
    #first try with formatdb                                                  
    if(seqType == "NT"): typeFlg = "-p F"
    else: typeFlg = "-p T"
    p = subprocess.Popen('formatdb -i ' + db_file + ' ' + typeFlg,
                         shell=True)
    p.wait()
    if(p.returncode != 0):
        #try with makeblastdb                                                 
        if(seqType == "NT"):
            typeFlg = "-dbtype nucl"
        else:
            typeFlg = "-dbtype prot"

        p = subprocess.Popen('makeblastdb -in ' + db_file + " " + typeFlg,
                             shell=True)
        p.wait()
        utility.die_on_error(p.returncode)

    #insert records into the sb_catalog
    i = catalog.insert()
    i.execute(begin=sbStart, end=curId - 1, fileID=fileID, projectID=projectID)

    #update the number of sequences in the project directory
    if(seqType == "NT"): session.conn.execute(directory.update().where(directory.c.projectID == projectID).values(num_NT_seqs=numSeqs))
    else: session.conn.execute(directory.update().where(directory.c.projectID == projectID).values(num_AA_seqs=numSeqs))

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


def storeSortedHomology(projectID, sortType):
    subprocess.Popen("perl \"" + perlPath + "storeSortedHomology.pl\" " + str(projectID) + " " + str(sortType), shell=True).wait()

def main():
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
    args = parser.parse_args()

    v = args.verbose
    reportStatus("Starting a new MySQLseq session with NeuroBase\n", v)
    s = utils.DBSession(args.user, args.passwd, 'localhost')

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

        if(args.swissprot):
            reportStatus("Loading swissprot homology data...\n",)
            p = Reader(baseName + "_blast_swissprot.txt", args.alnFmt, 2, v)

        elif(args.nr):
            reportStatus("Loading nr homology data.\n", v)
            p = Reader(baseName + "_blast_nr.txt", args.alnFmt, 1, v)


        p.read()
        loadHomology(projectID, linkDict, p, s, v, args.deleteH)
        assignBestAnnotation(projectID, s, v)
        reportStatus("Done.\n")
        
    elif(args.alignments):
        projectID = getProjectID(args.pn, s, v)
        reportStatus("Linking sequence identifiers in FASTA file to Neurobase ID\n", v)
        linkDict = utils.link_dbid_to_fastaid(projectID, s)
        reportStatus("Done.\n", v)
        if(args.swissprot):
            r = Reader(baseName + "_blast_swissprot.txt", args.alnFmt, 2, v)

        elif(args.nr):
            r = Reader(baseName + "_blast_nr.txt", args.alnFmt, 1, v)

        r.read()
        reportStatus("Loading annotation alignments.\n", v)
        loadHomologyAlignments(projectID, linkDict, r, s, v)


    elif(args.bestAnnotation):
    
        '''if(args.swissprot):
            reportStatus("Loading swissprot best homology data...\n",)
            p = Reader(baseName + "_blast_swissprot.txt", args.alnFmt, 2, v)

        elif(args.nr):
            reportStatus("Loading nr best homology data.\n", v)
            p = Reader(baseName + "_blast_nr.txt", args.alnFmt, 1, v)

        p.read()
        '''
        assignBestAnnotation(projectID, s, v)
        reportStatus("Done.\n")

    elif(args.initProject):
        initProject(args.publicName, s)

    elif(args.quant):
        projectID = getProjectID(args.pn, s, v)
        loadQuantification(projectID, args.pn, args.quantFile, args.abundanceCol, s, v)

    elif(args.KEGG):
        projectID = getProjectID(args.pn, s, v)
        #create linkage of assembly ID to neurobase ID
        if(v):print("Linking sequence identifiers in FASTA file to NeuroBase ID")
        linkDict = utils.link_dbid_to_fastaid(projectID, s)
        if(v):print("Done.")
        if(args.keggFile is None):
            args.keggFile = baseName + "_KEGG.txt"
        loadKEGG(projectID, linkDict, args.keggFile, s, v)

    elif(args.GO):
        projectID = getProjectID(args.pn, s, v)
        #create linkage of assembly ID to neurobase ID
        if(v):print("Linking sequence identifiers in FASTA file to NeuroBase ID")
        linkDict = utils.link_dbid_to_fastaid(projectID, s)
        if(args.goFile is None):
            args.goFile = baseName + "_GO.txt"
        loadGO(projectID, linkDict, args.goFile, s, v)
        
    elif(args.GOCategories):
        projectID = getProjectID(args.pn, s, v)
        catFile = baseName + "_gocats.txt"
        loadGOCategories(projectID, args.pn, catFile, s, v)

    elif(args.pfam):
        projectID = getProjectID(args.pn, s, v)
        linkDict = utils.link_dbid_to_fastaid(projectID, s)
        if(args.pfamFile is None):
            args.pfamFile = baseName + "_pfam.txt"
            
        loadPfam(projectID, args.pfamFile, args.seqType, linkDict, s)

    elif(args.loadSeqs):
        reportStatus("Loading project sequences...")
        loadSequences(args.pn, args.seqType, s, args.seqFile, v, sbStart=args.sbStart)
        reportStatus("done.")
    else:
        allPipe(args, s, v)

if __name__ == '__main__':
    main()
