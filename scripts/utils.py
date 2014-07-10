'''
Created on Aug 21, 2012

@author: maggiemcg
'''

from sqlalchemy import *
from sqlalchemy import Column
from sqlalchemy.ext.declarative import declarative_base
from autonomics import settings

Base = declarative_base()


def create_seq_table(tableName, session):
    return Table(tableName, session.meta, Column('sb_id', Integer, primary_key=True), Column('file_id', Integer), Column('NT_sequence', Text), Column("AA_sequence", Text), Column("description", Text), Column("type", String(2)), Column("length", Integer), Column("date", Date), Column("abundance", Integer), mysql_engine='InnoDB')


def createTableObject(tableName, session):
    return Table(tableName, session.meta, autoload=True, autoload_with=session.engine)


def link_dbid_to_fastaid(projectID, session):
        returnDict = {}
        t = globals()['createTableObject'](str(projectID) + '_sequences', session)
        results = t.select().execute()
        for row in results.fetchall():
            desc = row.description
            desc = desc.strip()
            returnDict[desc] = row.sb_id
        return returnDict


class DBSession:

    def __init__(self, u, p, h, db):
        self.host = h
        self.user = u
        self.passwd = p
        self.db = db
#        print "create_engine(",'mysql+mysqldb://' + u + ':' + p + '@' + h + '/' + self.db
        self.engine = create_engine('mysql+mysqldb://' + u + ':' + p + '@' + h + '/' + self.db)
        self.meta = MetaData(bind=self.engine)
        self.conn = None

    def activateConnection(self):
        self.conn = self.engine.connect()

    def close(self):
        self.conn.close()



