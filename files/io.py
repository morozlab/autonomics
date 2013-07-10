import os.path
import re
from sys import exit


type_casts = {"int": int, "str": str, "float": float}


def append_to_filename(f, to_add):
    base, ext = os.path.splitext(f)
    return base + to_add + ext


def build_fields(field_list):

    if(field_list is None):
       return None
    fields = []
    for i, field in enumerate(field_list):
        el = field.split(":")
        column_index = i
        cast = str
        if(len(el) == 2):
            #check if this is a cast or a column index
            if(el[1] in type_casts):
                cast = type_casts[el[1]]
            else:
                column_index = int(el[1]) - 1
        elif(len(el) == 3):
            cast = type_casts[el[2]]
            column_index = int(el[1]) - 1
        fields.append((el[0], column_index, cast))

    return fields


def convert_if_int(s):
    try:
        return int(s)
    except ValueError:
        return s


def convert_if_number(s):
    try:
        return float(s)
    except ValueError:
        return s


def list_to_str(l, sep=","):
    ret = []
    for el in l:
        ret.append(str(el))
        ret.append(sep)
    ret = ret[:len(ret) - 1]
    return "".join(ret)


class DataFile:

    '''
    Changes needed:
        - if records is empty, should iterrecords and apply filters to each record
        - records passing filtering should be written to a new filtered file
        - new filtered file should be set as file_name_filtered(time), old self.handle closed and reopened on the new file name

    '''

    def __init__(self, file_name):
        self.file_name = file_name
        self.handle = open(self.file_name, 'r')
        self.records = {}

    def __iter__(self):
        return self

    def next(self):
        line = self.handle.readline()
        if(line == ''):
            self.reset()
            raise StopIteration

        return self.make_record(line)

    def append_if_exists(self, key, value):
        if(key in self.records):
            self.records[key].append(value)
        else:
            self.records[key] = [value]

    def apply_filters(self, filters):
        if(filters is None):
            return
        for key, records in self.records.items():
            tmp = []
            for record in records:
                append = True
                for (field, op, value) in filters:
                    if(not op(getattr(record, field), value)):
                        append = False
                        break
                if(append):
                    tmp.append(record)
            if(len(tmp) > 0):
                self.records[key] = tmp
            else:
                del self.records[key]

    def close(self):
        self.handle.close()

    def flatten(self, num_records=None):
        for key, rows in self.records.items():
            new_rows = [rows[0]]
            if(not num_records is None):
                for row in rows[1:]:
                    if(len(new_rows) == num_records):
                        break
                    new_rows.append(row)

            self.records[key] = new_rows

    def iterrecords(self):
        self.reset()
        for line in self.handle:
            record = self._make_record(line)
            yield record

    def read(self):
        return self.handle.read()

    def readlines(self):
        return self.handle.readlines()

    def reset(self):
        self.handle.close()
        self.handle = open(self.file_name, 'r')


class Record:

    def __init__(self):
        key = None

    def __getitem__(self, name):
        return getattr(self, name)

    def __setitem__(self, name, value):
        setattr(self, name, value)


class TabFile(DataFile):

    def __init__(self, file_name, fields=None, key_columns=[0], sep="\t", load_records=False, comments="#"):
        DataFile.__init__(self, file_name)
        self.comments = comments
        self.format = "tab"
        self.field_names = set()
        self.sep = sep
        if(self.sep == " "):
            self.split_re = re.compile("\s+")
        else:
            self.split_re = re.compile(self.sep)
        if(fields is None):
            fields = self._make_fields()

        for (name, column, cast) in fields:
            self.field_names.add(name)

        self.field_list = fields
        self.key_columns = key_columns
        if(self.key_columns is None):
            self.key_columns = [0]
        if(load_records):
            self.load_records()
        # self._make_records()

    def __getitem__(self, key):
        return self.records[key]

    def __setitem__(self, key, value):
        self.records[key] = value

    def _make_fields(self):
        ret = []
        for line in self.handle:
            elements = self.split_re.split(line)
            for i, el in enumerate(elements):
                ret.append(("col_" + str(i), i, str))
            break

        return ret

    def _make_record(self, line):
        if(line.startswith("\n")):
            return None
        if(line.startswith(self.comments)):
            return None
        line = line.rstrip("\n")
        elements = self.split_re.split(line)
        r = Record()
        r.source = self.file_name
        for (name, column, cast) in self.field_list:
            setattr(r, name, cast(elements[column]))

        r.key = self._generate_key(line)
        return r

    def load_records(self):
        self.reset()
        for line in self.handle:
            row = self._make_record(line)
            if(row is None):
                continue
            self.append_if_exists(row.key, row)

    def _generate_key(self, line):
        el = self.split_re.split(line)
        return "-".join([el[column] for column in self.key_columns]).strip()

    def get_ids(self):
        if(self.records_loaded()):
            return self.records.keys()

        keys = set()
        for record in self:
            keys.add(record.key)

        return keys

    def get_record(self, key):
        return self.records[key]

    def merge_with(self, other_file):
        #iterate over the records in this file
        for key, rows in self.records.items():
            if(key in other_file.records):
                shared = True
            else:
                shared = False
            for row in rows:
                for (name, column, cast) in other_file.field_list:
                    attr_name = str(name) + "_" + os.path.split(other_file.file_name)[1]
                    if(shared):
                        other_row = other_file.records[key][0]
                        setattr(row, attr_name, getattr(other_row, name))
                    else:
                        setattr(row, attr_name, None)

        #add the new field information to this file
        for (name, column, cast) in other_file.field_list:
            attr_name = str(name) + "_" + os.path.split(other_file.file_name)[1]
            self.field_names.add(attr_name)
            self.field_list.append((attr_name, column, cast))

    def records_loaded(self):
        if(len(self.records) > 0):
            return True

        return False

    def reorder(self, column_order):
        new = []
        for column in column_order:
            new.append(self.field_list[column - 1])

        self.field_list = new

    def write_flattened(self, file_name):
        out = open(file_name, 'w')
        s = []
        out.write("key" + self.sep)
        out.write(list_to_str([name for (name, column, cast) in self.field_list]))
        out.write("\n")
        for key, records in self.records.items():
            out.write(key + self.sep)
            for record in records:
                out.write(list_to_str([getattr(record, name) for (name, column, cast) in self.field_list if(column not in self.key_columns)]))
                out.write(";")
            out.write("\n")

        out.close()

    def write_to_new(self, file_name):
        out = open(file_name, 'w')
        for key, records in self.records.items():
            for record in records:
                for i, (name, column, cast) in enumerate(self.field_list):
                    if(i == 0):
                        out.write(getattr(record, name))
                    else:
                        out.write(self.sep + str(getattr(record, name)))
                out.write("\n")
        out.close()


class SamFile(TabFile):

    fields = ["query:1", "flag:2:int", "reference:3", "position:4:int", "mapq:5:int", "cigar:6", "query_seq:10", "query_qual:11"]

    def __init__(self, file_path, load_records=False, comments="#"):
        #self, file_name, fields=None, key_columns=[0], sep="\t", load_records=False, comments="#"
        TabFile.__init__(self, file_path, fields=build_fields(self.fields))

    def iterrecords(self):
        self.reset()
        for line in self.handle:
            if(not line.startswith("@")):
                record = self._make_record(line)
                yield record
