'''
Created on Oct 18, 2012

@author: Administrator

This module is designed to be run as a script, and includes functionality to make a KEGG term loadfile which can be loaded to a MySQL database.

usage: make_kegg_loadfile.py [-h]
                             genes_to_ko genes_to_uniprot ko_to_pathway
                             map_title

positional arguments:
  genes_to_ko       File associating KO terms to gene accessions.
  genes_to_uniprot  File associating uniprot accessions with gene accessions.
  ko_to_pathway     File associating KO terms to KEGG pathway ids.
  map_title         File associating KEGG pathway ids to pathway descriptions

optional arguments:
  -h, --help        show this help message and exit

'''

import argparse
import sys

def get_field_val(field, sep=":"):
    if(sep in field):
        return field.split(sep)[1]
    else:
        return field

def hashify_file(file, h, key_column, value_column, single_value = True):
    f = open(file, 'r')
    for line in f:
        elements = line.strip().split("\t")
        key = get_field_val(elements[key_column])
        value = get_field_val(elements[value_column])
        if(single_value):
            if(not key in h):
                h[key] = value
        else:
            if(key in h):
                h[key].append(value)
            else:
                h[key] = [value]

    f.close()


parser = argparse.ArgumentParser()

parser.add_argument('genes_to_ko', help='File associating KO terms to gene accessions.')
parser.add_argument('genes_to_uniprot', help='File associating uniprot accessions with gene accessions.')
parser.add_argument('ko_to_pathway', help='File associating KO terms to KEGG pathway ids.')
parser.add_argument('map_title', help='File associating KEGG pathway ids to pathway descriptions')

args = parser.parse_args()

genes_to_ko = {}
genes_to_uniprot = {}
ko_to_pathways = {}
pathway_to_desc = {}

hashify_file(args.genes_to_ko, genes_to_ko, 0, 1)
hashify_file(args.genes_to_uniprot, genes_to_uniprot, 0, 1)

f = open(args.ko_to_pathway, 'r')
for line in f:
    elements = line.strip().split("\t")
    key = get_field_val(elements[0])
    value = get_field_val(elements[1]).strip("ko")
    if(key in ko_to_pathways):
        ko_to_pathways[key].append(value)
    else:
        ko_to_pathways[key] = [value]


hashify_file(args.map_title, pathway_to_desc, 0, 1)

#generate the uniprot_to_ko hash
uniprot_to_ko = {genes_to_uniprot[gene] : ko for gene, ko in genes_to_ko.items() if gene in genes_to_uniprot}

#iterate over it, associating the other info
for acc, ko in uniprot_to_ko.items():
    if(ko in ko_to_pathways):
        for path_id in ko_to_pathways[ko]:
            if(path_id in pathway_to_desc):
                print(acc + "\t" + ko + "\t" + path_id + "\t" + pathway_to_desc[path_id])
            else:
                print(acc + "\t" + ko + "\t" + path_id + "\tNA")
