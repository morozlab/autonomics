'''
Created on Aug 9, 2012

@author: administrator
'''

from files.filetools import generate_key
import sys

def dict_from_file(file_name, key_columns=[0], sep='\t'):
    fh = open(file_name, 'r')
    return_dict = {}
    for l in fh:
        l = l.rstrip("\n")
        key = generate_key(l,key_columns)
        return_dict[key] = l.split(sep)
    return return_dict


def dictToFile(file_name, dict):
    fh = open(file_name, 'w')


def flatten(records, key_columns):

    ret = {}


def merge(file1, file2, key_columns):

    ret = []

    f1 = dict_from_file(file1, key_columns)
    f2 = dict_from_file(file2, key_columns)

    for key, value in f1.items():
        ret.append(value.extend(f2[key]))
        print(ret)
        sys.exit(-1)

    return ret
