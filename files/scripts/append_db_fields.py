import argparse
import MySQLdb
import Databases

parser = argparse.ArgumentParser(description="This script appends database information in tab-format to the end of an existing file")
parser.add_argument('-f', '--input-file', dest='inputFile', required=True)
parser.add_argument('-db', '--database-name', dest='database', required=True)
parser.add_argument('-c', '--identifier-column', dest='column', required=True)
parser.add_argument('-d', '--id-is-description', dest='idDescription', default=False, const=True, action='store_const')
parser.add_argument('-p', '--project-id', dest='projectID', default='-1')
parser.add_argument('--homology', dest='homology', default=False, const=True, action='store_const')
parser.add_argument('--abundance', dest='abundance', default=False, const=True, action='store_const')
parser.add_argument('--gene-ontology', dest='GO', default=False, const=True, action='store_const')

args = parser.parse_args()

#connect to the specified database
if(args.database == 'neurobase' or args.database == 'NeuroBase'):
    db = MySQLdb.connect(host='localhost', user='root', passwd='meow12', db='moroz_lab')

#open the file and generate a list of the identifiers
ids = []
textFile = open(args.inputFile, 'r')
for line in textFile:
    elements = line.split("\t")
    ids[elements[args.column - 1]] = elements

#check if we're getting homology, abundance, or gene ontology data
if(args.homology):
    for seqID in ids:
        fields = "*"
        lookupID = seqID
        if(args.idDescription):
            query = "SELECT sb_id FROM " + args.projectID + "_sequences WHERE description LIKE '%" + seqID + "'"
            cursor = db.cursor()
            cursor.execute(query)
            lookupID = cursor.fetchone()[0]
        
        table = "annotation_db, sorted_homology"
        where = "project_id='" + args.projectID + "' AND annotation_db.annotation_id = sorted_homology.annotation_id AND sorted_homology.sb_id='" + lookupID + "'"
         
        homologyResults = Databases.NeuroBaseUtils.RetrieveData(db, fields, table, where)
        
        


