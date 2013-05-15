
'''
Created on Aug 21, 2012

@author: Mathew Citarella

This is a utility module containing classes and methods for Autonomics to deal with data and objects over networks.

'''

import os
import subprocess
import paramiko
import time
from sqlalchemy.exc import IntegrityError
from sqlalchemy import *
from sqlalchemy import Column
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.sql import and_
from zeroclick import settings, utility
import sys

Base = declarative_base()

def ssh_connect(ssh_credentials, retries=None):
    '''
        ssh_credentials (settings.Credentials):
            An object holding the connection details for the to-be-created
            SSHConnection

        retries (int):
            How many times should this method attempt to create an SSHConnection
            object before raising an exception

        Creates and returns a SSHConnection object using the provided
        credentials. If an exception is raised during the connection
        process, will retry up to n=retries times, before raising
        the exception.

    '''
    if(retries is None):
        retries = settings.SSH_RETRIES
    try:
        return SSHConnection(
                             ssh_credentials.host,
                             username=ssh_credentials.user,
                             password=ssh_credentials.passwd
                             )
    except:
        sys.stderr.write("".join["Error creating SSHConnection to ",
                                 str(ssh_credentials.host),
                                 ", sleeping for two minutes and trying \ "
                                 "again.\n"])
        time.sleep(120)
        if(retries > 0):
            return ssh_connect(ssh_credentials, retries - 1)
        else:
            raise


def check_and_insert(s, insert):
    ''' s (Select):
            SQLAlchemy Select object created using either table.select() or
            select() from sqlalchemy.sql

        insert (Insert):
            SQLAlchemy Insert object created using either table.insert() or
            insert() from sqlalchemy.sql

        Uses the select method to check for a database entry. If the entry
        does not exist, it calls i.execute() and returns the lastrowid from
        the insert results.

    '''
    check = s.execute().fetchone()
    if(check is None):
        result = insert.execute()
        return result.lastrowid

    return None


def check_rsync_error(msg):
    '''Returns True if msg contains error msg, False otherwise.'''

    if("rsync error" in msg):
        return True
    return False


def create_seq_table(tableName, session):
    ''' tableName (str):
            String name of the table you wish to create.

        session (DBSession):
            Active session with some relational database instance.

        Creates the table specified by tableName in the MySQL database
        associated with session and returns a Table object.

        Uses this table definition for tableName:
            Column('sb_id', Integer(20), primary_key=True),
            Column('file_id', Integer), Column('NT_sequence', Text),
            Column("AA_sequence", Text), Column("description", Text),
            Column("type", String(2)), Column("length", Integer),
            Column("date", Date), Column("abundance", Integer(10)),
            Column("best_annot", Text), Column("best_annot_eval", Float),
            mysql_engine='InnoDB')

    '''
    return Table(
                 tableName, session.meta,
                 Column('sb_id', Integer(20), primary_key=True),
                 Column('file_id', Integer), Column('NT_sequence', Text),
                 Column("AA_sequence", Text), Column("description", Text),
                 Column("type", String(2)), Column("length", Integer),
                 Column("date", Date), Column("abundance", Integer(10)),
                 Column("best_annot", Text), Column("best_annot_eval", Float),
                 mysql_engine='InnoDB'
                 )


def delete_project(pid, session):
    '''
        pid (int):
            Identifier of the project to delete

        session (netutils.DBSession):
            A session object with an active connection to the Autonomics DB

        Deletes the project specified by the given project_id from the database.

    '''
    pa = get_table_object("project_adapters", session)
    pn_mapping = get_table_object("pn_mapping", session)
    jn_mapping = get_table_object("jn_mapping", session)
    config = get_table_object("configuration", session)
    rn_2_pid = get_table_object("runname_to_pid", session)
    args = get_table_object("args", session)
    trans = session.conn.begin()
    try:
        jn_results = session.conn.execute(jn_mapping.select(jn_mapping.c.project_id == pid))
        for jn_row in jn_results.fetchall():
            # remove the args for each project
            session.conn.execute(args.delete(args.c.job_id == jn_row.job_id))

            # remove adapter entries
            session.conn.execute(pa.delete(pa.c.project_id == pid))

            # remove all the jn_mapping data for this project
            session.conn.execute(jn_mapping.delete(jn_mapping.c.project_id == pid))

            # remove all the configuration data for this project
            session.conn.execute(config.delete(config.c.project_id == pid))

            # remove the pn_mapping for this project
            session.conn.execute(pn_mapping.delete(pn_mapping.c.project_id == pid))

            # change the rn_2_pid information for this project
            session.conn.execute(rn_2_pid.delete(rn_2_pid.c.project_id==pid))
            trans.commit()
            return True
    except:
        trans.rollback()
        return False


def get_adapter_rows(pid, session, retries=settings.SSH_RETRIES):
    ''' pid (int):
            Integer project_id of the project adapters should be retrieved for.

        session (DBSession):
            Active session with some relational database instance.

        retries (int):
            Integer number of times this method should try to query the database for adapters before propogating the underlying error.

        Method looks up pid in known_adapters table, and then gets adapter sequences for all adapters assigned to the project by querying the known_adapters table. Returns a list of string representations of adapter sequences.
    '''
    try:
        ret = []
        ka = getTableObject("known_adapters", session)
        pa = getTableObject('project_adapters', session)
        results = select([ka.c.adapter_sequence, ka.c.end],
                         and_(ka.c.adapter_id==pa.c.adapter_id, pa.c.project_id==pid)).execute()
        return [row for row in results.fetchall()]

    except Exception as e:
        print("Error getting adaptors, sleeping and re-opening connection")
        time.sleep(10)
        session = make_db_session()
        get_adapters(pid, session, retries - 1)


def get_file(remote_path, local_path, credentials,
             xfer_mthd='rsync', tries=5):
    if(xfer_mthd == 'rsync'):
        remote_full = credentials.user + "@" + credentials.host + ":" + remote_path
        rsync(remote_full, local_path)

    elif(xfer_mthd == 'scp'):
        remote_full = credentials.user + "@" + credentials.host + ":" + remote_path
        scp(remote_full, local_path)

    elif(xfer_mthd == 'sftp'):
        if(credentials['port'] is None):
            ssh = SSHConnection(credentials.host,
                                username=credentials.user,
                                password=credentials.passwd)
        else:
            ssh = SSHConnection(credentials.host,
                                username=credentials.user,
                                password=credentials.passwd,
                                port=credentials.port)

        ssh.get(remote_path, local_path)
        ssh.close()

    if(os.path.exists(local_path)):
        return True

    return False


def get_field_list(table):
    ''' table_name (Table):
            Table object representing a reflection of a db table.

        Returns a list of field names for the table provided.

    '''

#    column_tbl = get_table_object("information_schema.columns", session)
#    res = column_tbl.select(column_tbl.c.table_name==table_name).execute()
    ret = []
    for column in table.c:
        ret.append(column.name)

    return ret


def get_jid(pid, job_type, session):
    ''' pid (int):
            Integer project_id of the project job is associated with.

        job_type (str):
            Type of the job whose job_id you are trying to find.

        session (DBSession):
            Active session with some relational database instance.

        Looks up the job_id of the job specified by pid and job_type in the jn_mapping table. Returns job_id as an integer if the job is found, otherwise returns None.

    '''
    jn_mapping = get_table_object("jn_mapping", session)
    res = jn_mapping.select(and_(jn_mapping.c.project_id==pid, jn_mapping.c.job_type==job_type)).execute()
    row = res.fetchone()
    if(row is None):
        return None
    return row.job_id


def get_jobtype_inputs(job_type, session):

    j_inputs = get_table_object("job_type_inputs", session)
    session.activateConnection()
    res = session.conn.execute(j_inputs.select(j_inputs.c.job_type==job_type))
    no_split = []
    inputs = {}
    for row in res.fetchall():
        inputs[row.flag] = row.extension_for_default
        if(row.split_on_remote == 'N'):
            no_split.append(row.flag)

    print(inputs, no_split)
    return (inputs, no_split)


def get_max_sbid(session):
    '''
        session (netutils.DBSession):
            A session object with an active connection to the Autonomics database.

        Returns the maximum sb_id assigned to any sequence across this database instance.
    '''

    sb_catalog = get_table_object('sb_catalog', session)
    q = select([func.max(sb_catalog.c.end).label('max')])
    res = session.conn.execute(q)
    max_row = res.fetchone()
    return max_row.max


def get_pid(project_name, session):
    ''' project_name (str):
            Name of the project to lookup.

        session (DBSession):
            Active session with some relational database instance.

        Looks up project_name in pn_mapping, returns project_id if project_name returns a result, returns None otherwise.

    '''
    pn_mapping = get_table_object("pn_mapping", session)
    res = pn_mapping.select(pn_mapping.c.project_name==project_name).execute()
    row = res.fetchone()
    if(row is None):
        return None

    return row.project_id


def get_table_object(tableName, session):
    ''' tableName (str):
            Name of the table to reflect with SQLAlchemy

        session (DBSession):
            Active session with some relational database instance

        Returns a SQLAlchemy Table object by reflecting the table using the tableName and the session metadata.

    '''
    return Table(tableName, session.meta, autoload=True, autoload_with=session.engine)


def get_rowfield_values(row_obj, field_list, as_dict=False):
    '''
        row_obj (RowProxy):
            A RowProxy object returned from a sqlalchemy result set, either as
            the result of a call to results.fetchone() or results.fetchall()

        field_list ([str]):
            List of fieldnames for which to return values

        as_dict (Boolean):
            Instead of returning a list of values, return a dict of
            field->value pairs.

            Defaults to False

        Returns a list of values from the row object for the specified list of fields.
    '''
    if(as_dict):
        ret = {field_list[i] : getattr(row_obj, field_list[i]) for i in field_list if hasattr(row_obj, field_list[i])}
    else:
        ret = [getattr(row_obj, field_list[i]) for i in field_list if hasattr(row_obj, field_list[i])]
    return ret


def link_dbid_to_fastaid(project_id, session):
    ''' project_id (int):
            project_id of the project to link.

        session (DBSession):
            Active session with some relational database instance.

        Looks up the sequence table for the given project_id in Neurobase, returns a dictionary with key, value => description, sb_id from the sequence table.
    '''
    return_dict = {}
    t = globals()['createTableObject'](str(project_id) + '_sequences', session)
    results = t.select().execute()
    for row in results.fetchall():
        desc = row.description
        desc = desc.strip()
        return_dict[desc] = row.sb_id
    return return_dict


def make_db_session(db_name=settings.ZC_DB_NAME):
    ''' db_name (str):
            Name of the database to open a connection with.

        Returns a DBSession object by connecting with the system db credentials to db_name.

    '''
    return DBSession(settings.db_cred.host, u=settings.db_cred.user, p=settings.db_cred.passwd, d=db_name)


def put_file(local_path, remote_path, credentials, xfer_mthd='rsync', tries=5):
    '''
        local_path (str):
            Path to the source file

        remote_path (str):
            Path to the remote file

        credentials (settings.Credentials):
            Credentials storing the means of connection to the remote
            computer

        xfer_method (str)='rsync':
            The method used to transfer the file

            Defaults to 'rsync'

            Supported values: 'rsyc', 'sftp', 'scp'

        tries (int)=5:
            How many times to re-attempt the transfer when an exception is
            raised during the transfer process. If attempts = tries the
            exception is raised

        Transfers file local_path to remote_path, using the specified
        transfer method. Note that using rsync or scp as the transfer
        method requires that public key authentication be setup for the
        user account attempting the transfer.
    '''
    if(xfer_mthd == 'rsync'):
        remote_full = credentials['user'] + "@" + credentials['host'] + ":" + remote_path
        rsync(local_path, remote_full)

    elif(xfer_mthd == 'scp'):
        remote_full = credentials['user'] + "@" + credentials['host'] + ":" + remote_path
        scp(local_path, remote_full)

    elif(xfer_mthd == 'sftp'):
        while(True):
            try:
                if(credentials['port'] is None):
                    ssh = SSHConnection(credentials.host, username=credentials.user, password=credentials.passwd)
                else:
                    ssh = SSHConnection(credentials.host, username=credentials.user, password=credentials.passwd, port=credentials.port)
                ssh.put(local_path, remote_path)
                break
            except Exception as e:
                print("Error in put_file of: " + local_path)
                print(e.message)
                if(tries == 0):
                    raise
                print("Sleeping and then trying again.")
                tries -= 1
                time.sleep(settings.PROCESS_MONITOR_INTERVAL)

        ssh.close()

    if(os.path.exists(local_path)):
        return True

    return False


def remote_path_exists(conn, path):
    ''' conn (SSHConnection):
            Connection to some remote machine.

        path (str):
            Full path to the directory.

        Returns True if the remote path exists on the machine accessed via conn, False otherwise.

    '''
    results = utility.convert_if_int(conn.execute("[ -e " + path + " ] && echo 1 || echo 0")[0].strip())
    if(results == 1):
        return True
    return False


def row2dict(db_row):
    ''' Converts a RowProxy element to a more lightweight native dict object.
    '''
    ret = {}
    for column in db_row.keys():
        ret[column] = getattr(db_row, column)

    return ret

def rsync(src, dest):
    ''' Rsyncs 'src' to 'dest'. Arguments to rsync are: -acvzP.
    '''
    p = subprocess.Popen(["rsync", "-acvzP", src, dest], stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if(p.returncode != 0):
        raise IOError(stderr)


def scp(src, dest):
    ''' Scp's 'src' to 'dest'. Arguments to scp are: -Crv.
    '''
    p = subprocess.Popen(["scp", "-Crv", src, dest], stderr=subprocess.PIPE)
    out, err = p.communicate()
    if(p.returncode != 0):
        raise IOError(err)


def scp_error(std_err):
    ''' Checks std_err for presence of scp error msg. Returns True if found,
        False otherwise.
    '''
    if("No such file" in std_err):
        return True
    if("ssh:" in std_err):
        return True
    if("Permission denied" in std_err):
        return True

    return False


def update_on_exists(insert, update):
    ''' insert:
            SQLAlchemy Insert object.

        update
            SQLAlchemy Update object.

        Attempts to execute the provided insert. If an IntegrityError is raised, perform the update provided.

    '''
    try:
        insert.execute()
    except IntegrityError:
        update.execute()


"""Friendly Python SSH2 interface."""

class SSHConnection(object):
    """Connects and logs into the specified hostname.
    Arguments that are not given are guessed from the environment."""

    def __init__(self,
                 host,
                 username = None,
                 private_key = None,
                 password = None,
                 port = 22,
                 ):

        if not username:
            username = os.environ['LOGNAME']

        port = int(port)
        self.host = host
        self.user = username
        self.passwd = password
        self.private_key = private_key

        self._sftp_live = False
        self._sftp = None

        # Log to a temporary file.
        #templog = tempfile.mkstemp('.txt', 'ssh-')[1]
        #paramiko.util.log_to_file(templog)

        # Begin the SSH transport.
        self._transport = paramiko.Transport((host, port))
        self._tranport_live = True

        # Authenticate the transport.
        if password:
            # Using Password.
            self._transport.connect(username = username, password = password)

        else:
            # Use Private Key.
            if not private_key:
                # Try to use default key.
                if os.path.exists(os.path.expanduser('~/.ssh/id_rsa')):
                    private_key = '~/.ssh/id_rsa'
                elif os.path.exists(os.path.expanduser('~/.ssh/id_dsa')):
                    private_key = '~/.ssh/id_dsa'
                else:
                    raise TypeError, "You have not specified a password or key."

            private_key_file = os.path.expanduser(private_key)
            rsa_key = paramiko.RSAKey.from_private_key_file(private_key_file)
            self._transport.connect(username = username, pkey = rsa_key)

    def _sftp_connect(self):
        """Establish the SFTP connection."""
        if not self._sftp_live:
            self._sftp = paramiko.SFTPClient.from_transport(self._transport)
            self._sftp_live = True

    def get(self, remotepath, localpath = None):
        """Copies a file between the remote host and the local host."""
        if not remotepath:
            remotepath = os.path.split(localpath)[1]
        self._sftp_connect()
        self._sftp.get(remotepath, localpath)

    def put(self, localpath, remotepath = None):
        """Copies a file between the local host and the remote host."""
        if not remotepath:
            remotepath = os.path.split(localpath)[1]

        self._sftp_connect()
        self._sftp.put(localpath, remotepath)

    def execute(self, command):
        """Execute the given commands on a remote machine."""
        retries = 5
        while(retries > 0):
            try:
                channel = self._transport.open_session()
                break
            except:
                retries -= 1
                time.sleep(1)
                self = SSHConnection(self.host, username=self.user, password=self.passwd)

        channel.exec_command(command)
        output = channel.makefile('rb', -1).readlines()
        if output:
            return output
        else:
            return channel.makefile_stderr('rb', -1).readlines()

    def close(self):
        """Closes the connection and cleans up."""
        # Close SFTP Connection.
        if self._sftp_live:
            self._sftp.close()
            self._sftp_live = False
        # Close the SSH Transport.
        if self._tranport_live:
            self._transport.close()
            self._tranport_live = False

    def __del__(self):
        """Attempt to clean up if not explicitly closed."""
        self.close()


class DBSession:
    ''' DBSession - Class representing a SQLAlchemy connection to an arbitrary relational database.

        Important Attributes:
            self.engine - stores result of create_engine()
            self.meta - stores result of MetaData(bind=self.engine)
            self.conn - initially set to None, until activateConnection is invoked
    '''

    retries = 5

    def __init__(self, h, d, u, p, driver="mysql+oursql", pool_size=1):
        ''' h (str):
                Host.
            d (str):
                Database name to connect to on Host.

            u (str):
                Username to use during connection attempts.

            p (str):
                Password to use during connection attempts.

            driver (str):
                SQLAlchemy-formatted driver string for the relational database.

            Connects to the database using the supplied credentials and sets the meta attribute.
        '''

        self.host = h
        self.driver = driver
        self.db = d
        self.user = u
        self.passwd = p

        self.engine = create_engine(driver + '://' + u + ':' + p + '@' + h + '/' + d, pool_recycle=3600, pool_size=pool_size)
        self.meta = MetaData(bind=self.engine)
        self.conn = self.engine.connect()
#        print("Connection created to: " + self.host)
        #    def activateConnection(self):
        '''
            Calls self.engine.connect() and stores the return value in self.conn, effectively activating conn for use.
        '''
        # self.conn = self.engine.connect()

    def close(self):
        '''
            If the object's connection isn't None, closes the connection.
        '''
        if(not self.conn is None):
            self.conn.close()

    def safe_exec(self, query):
        ''' Deprecated. To be removed in a future release.
        '''
        my_retries = DBSession.retries
        while(True):
            try:
                self.conn.execute(query)
                break
            except Exception as e:
                print("Exception in safe_exec.")
                print(e.message)
                my_retries -= 1
                if(my_retries < 0):
                    raise
                time.sleep(1)
                self = DBSession(self.host, self.db, self.user, self.passwd, self.driver)
                self.activateConnection()
