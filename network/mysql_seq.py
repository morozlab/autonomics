#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Administrator
#
# Created:     09/08/2012
# Copyright:   (c) Administrator 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import MySQLdb
from sqlalchemy import *


class session:

    def __init__(self, user, passwd, output = "retriever_output.txt", ):
        self.user = user
        self.passwd = passwd
        self.output = output
        self.engine = create_engine('mysql+mysqldb://' + self.user + ':' + self.passwd + '@localhost/moroz_lab')

    def _connect(self):
        return  MySQLdb.connect("localhost", self.user, self.passwd, "moroz_lab")

    def execute(self, query, params = None):
        '''
        if(query.startswith("SELECT") or query.startswith("select")):
            if(params): return self.getData(query, params)
            else: return self.getData(query)
        else:
            c = self._connect().cursor(MySQLdb.cursors.DictCursor)
            if(params): c.execute(query, params)
            else: c.execute(query)
            c.close()
            return None
        '''

    def getData(self, query, params = None):
        c = self._connect().cursor(MySQLdb.cursors.DictCursor)
        print(query)
        if(params): c.execute(query, params)
        else: c.execute(query)
        c.close()
        l = [r for r in c.fetchall()]
        print(l)
        return l


    def linkDatabaseIDToFastaID(self, projectID):
        returnDict = {}
        query = "SELECT sb_id, description FROM " + str(projectID) + "_sequences"
        c = self._connect()
        cursor = c.cursor(MySQLdb.cursors.DictCursor)
        cursor.execute(query)
        row = cursor.fetchone()
        while(row != None):
            returnDict[row['description']] = row['sb_id']
            row = cursor.fetchone()
        cursor.close()
        return returnDict

    def setOutput(self, outputName):
        self.output = outputName

    '''
    def getIDsFromFile(self, file, fmt, projectID):
        c = self._connect()
        if(fmt == 'fasta'):
            ids = []
            for record in SeqIO.parse(file, fmt):
                #get the neurobase ID for this sequence
                query = "SELECT sb_id FROM " + str(self.projectID) + "_sequences WHERE description "
    '''
def main():
    pass

if __name__ == '__main__':
    main()
