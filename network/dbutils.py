'''
Created on Aug 21, 2012

@author: maggiemcg
'''

from sqlalchemy import *
from sqlalchemy import Column
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

def createSeqTable(tableName, session):
    return Table(tableName, session.meta, Column('sb_id', Integer(20), primary_key=True), Column('file_id', Integer), Column('NT_sequence', Text), Column("AA_sequence", Text), Column("description", Text), Column("type", String(2)), Column("length", Integer), Column("date", Date), Column("abundance", Integer(10)), Column("best_annot", Text), Column("best_annot_eval", Float), mysql_engine='InnoDB')

def createTableObject(tableName, session):
    return Table(tableName, session.meta, autoload=True, autoload_with=session.engine)


def linkDatabaseIDToFastaID(projectID, session):
        returnDict = {}
        t = globals()['createTableObject'](str(projectID) + '_sequences', session)
        results = t.select().execute()
        for row in results.fetchall():
            desc = row.description
            desc = desc.strip()
            returnDict[desc] = row.sb_id
        return returnDict


class DBSession:

    def __init__(self, h, d, u, p, driver="mysql+oursql"):
        self.engine = create_engine(driver + '://' + u + ':' + p + '@' + h + '/' + d)
        self.meta = MetaData(bind=self.engine)
        self.conn = None

    def activateConnection(self):
        self.conn = self.engine.connect()

    def close(self):
        self.conn.close()



