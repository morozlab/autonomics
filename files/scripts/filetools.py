'''
Created on Oct 18, 2012

@author: Administrator
'''
import argparse
from Bio import SeqIO
from files.io import TabFile, build_fields, convert_if_int, convert_if_number
import operator
import os
import sys

type_casts = {"int": int, "str": str, "float": float}

def add_if_exists(dict_of_sets, key, record):
    if(key in dict_of_sets):
        dict_of_sets[key].add(record)
    else:
        dict_of_sets[key] = set([record])


def append_if_exists(dict_of_lists, key, record):
    if(key in dict_of_lists):
        dict_of_lists[key].append(record)
    else:
        dict_of_lists[key] = [record]


def build_filter(filter_string):
    el = filter_string.split(":")
    return (el[0]), comparison_macros[el[1]], convert_if_number(el[2])


def convert_to_indices(columns):
    return [(convert_if_int(column) - 1) for column in columns]


def generate_key(line, key_columns, sep="\t"):
    el = line.split(sep)
    return "-".join([el[column] for column in key_columns]).strip()


def make_column_list(columns, sep=","):
    return map(convert_if_int, columns.split(sep))


def make_tab_files(files, fields, key_columns, sep='\t'):
    return [TabFile(f, fields, key_columns, sep) for f in files]


comparison_macros = {"lt": operator.lt, "gt": operator.gt, "eq": operator.eq, "ne": operator.ne, "ge": operator.ge, "le": operator.le}

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('files', nargs='+', help='Input files.')
    parser.add_argument('-f', '--format', dest='format', help='File format.', default='tab')
    parser.add_argument('--fields', dest='fields', nargs="+", help='Fields present in the tab-delimited file, specified in the form "field:column:cast" - cast defaults to str()', default=None)
    parser.add_argument('--fields2', dest='fields2', nargs="+", default=None, help='Fields for a second file in two-file operations.')
    parser.add_argument('--key-col', dest='key_columns', nargs="+", default="0", help='List of columns on which the key should be built for every record in the file.')
    parser.add_argument('--key-col2', dest='key_columns2', nargs='+', default=None, help='List of columns for a second file (during merges, other two file-operations) to generate a key on.')
    parser.add_argument('--key-fields', dest='key_fields', nargs='+', default=None, help='List of field names on which the key should be built for every record in the file.')
    parser.add_argument('--count-unique-ids', dest='count_unique_ids', default=False, const=True, action='store_const', help='Determine the number of unique identifiers present across all files.')
    parser.add_argument('--flatten', dest='flatten', help='Flatten a tab-delimited file on the specified key.', default=False, const=True, action='store_const')
    parser.add_argument('--flatten-depth', dest='flatten_depth', default=1, help="Maximum number of records for each key after a flatten operation.")
    parser.add_argument('--merge', dest='merge', default=False, const=True, action='store_const')
    parser.add_argument('--remove-duplicate-on', dest='remove_duplicate_on', default=False, const=True, action='store_const')
    parser.add_argument('--reorder', dest='reorder', default=None, nargs='+', help='Reorder a tab-delimited file. Argument is a list specifing the new column order, e.g. 1 4 2 5 3.')
    parser.add_argument('--rewrite-collapsed', dest='rewrite_collapsed', default=False, const=True, action='store_const')
    parser.add_argument('--separator', dest='separator', default='\t')
    parser.add_argument('--subset-of-columns', dest='subset_of_columns', default=False, const=True, action='store_const')
    parser.add_argument('--column-list', dest='column_list', default=None, help="Comma-separated list of columns to keep/work on")
    parser.add_argument('--intersection', dest='intersection', default = False, const = True, action='store_const', help='Create a file containing the intersection of the supplied files, based on keys constructed for each record using the columns in column_list.')
    parser.add_argument('--filter', dest='filter_list', nargs="+", default=None, help='Filters in the form "column:op:value". Op can be one of the following:\n>\n<\n=\n!=')
    parser.add_argument('-v', '--verbose', dest='verbose', default=False, const=True, action='store_const')

    args = parser.parse_args()

    if(args.format == 'tab'):
        #convert the supplied columns to list indices
        args.key_columns = convert_to_indices(args.key_columns)
        if(len(args.files) == 1):
            files = [TabFile(args.files[0], build_fields(args.fields), args.key_columns, sep=args.separator)]
        elif(len(args.files) == 2):
            #convert the second set of columns
            args.key_columns2 = convert_to_indices(args.key_columns2)
            files = [TabFile(args.files[0], build_fields(args.fields), args.key_columns, sep=args.separator), TabFile(args.files[1], build_fields(args.fields2), args.key_columns2, sep=args.separator)]
        else:
            files = make_tab_files(args.files, build_fields(args.fields), args.key_columns, sep=args.separator)

    if(not args.filter_list is None):
        filters = map(build_filter, args.filter_list)
        for f in files:
            f.load_records()
            f.apply_filters(filters)

    if(args.count_unique_ids):
        keys = set()
        for f in files:
            for record in f.iterrecords():
                keys.add(record.key)
                if(args.verbose):
                    print(len(keys))

        print(len(keys) + " unique identifiers across all files.")

    if(args.subset_of_columns):
        if(args.format == 'tab'):
            columns = args.column_list.split(",")
            columns = map(convert_if_int, columns)
            for f in args.files:
                out = open(f + ".subset", 'w')
                infile = open(f, 'r')
                for line in infile:
                    out.write("\t".join([field for i, field in enumerate(line.strip().split("\t")) if i + 1 in columns]) + "\n")
                out.close()
                infile.close()

    if(args.remove_duplicate_on):
        keycolumns = make_column_list(args.column_list)
        for f in args.files:
            if(args.format == 'tab'):
                with open(f, 'r') as file:
                    to_print = {}
                    for line in file:
                        line = line.strip()
                        if(len(line) == 0):
                            continue
                        elements = line.split(args.separator)
                        key = "".join([elements[column - 1] for column in keycolumns])
                        if(not key in to_print):
                            to_print[key] = line
                    for line in to_print.values():
                        print(line)
                file.close()

    if(args.rewrite_collapsed):
        for fn in files:
            fn.load_records()
            fn.write_flattened(fn.file_name + "_flattened.txt")

    if(args.intersection):
        key_columns = make_column_list(args.column_list)
        files_with_key = {}
        for f in args.files:
            with open(f, 'r') as cur_file:
                for line in cur_file:
                    line = line.rstrip("\n")
                    key = generate_key(line, key_columns)
                    add_if_exists(files_with_key, generate_key(line, key_columns), f)


        intersections = {}
        for key, filenames in files_with_key.items():
            #new key is combination of files that have this record
            file_key = "-".join(filenames)
            #add this record to the list of records shared by this file combination
            add_if_exists(intersections, file_key, key)

        del files_with_key

        #print all of the intersection information
        for file_combo, record_keys in intersections.items():
            print(file_combo)
            print("----------")
            for key in record_keys:
                print(key)

    if(args.merge):
        f = files[0]
        f2 = files[1]
        f.load_records()
        f2.load_records()
        f.merge_with(f2)
        for key, records in f.records.items():
            for record in records:
                first = True
                str_rep = ""
                for(name, column, cast) in f.field_list:
                    value = str(getattr(record, name))
                    if(first):
                        str_rep += value
                        first = False
                    else:
                        str_rep += "\t" + value
                    #print("\t".join(str(getattr(record, name)) for (name, column, cast) in f.field_list))
                print(str_rep)

    if(args.flatten):
        for f in files:
            f.load_records()
            f.flatten(convert_if_int(args.flatten_depth))
            f.write_to_new(f.file_name + "_flattened.txt")

    if(not args.reorder is None):
        for f in files:
            args.reorder = map(convert_if_int, args.reorder)
            f.reorder(args.reorder)
            f.load_records()
            f.write_to_new(f.file_name + "_reordered.txt")

if __name__ == '__main__':
    main()
