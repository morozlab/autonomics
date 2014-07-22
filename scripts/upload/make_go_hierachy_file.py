'''
Created on Oct 25, 2012

@author: Administrator

This script parses an OBO file from the Gene Ontology and produces a tab-delimited file containing relationships between parent and child objects.

usage: make_go_hierachy_file.py [-h] --obo OBO_FILE [--output OUTPUT_FILE]

optional arguments:
  -h, --help            show this help message and exit
  --obo OBO_FILE
  --output OUTPUT_FILE

'''

import argparse
import os
import re
import sys


def writeRelationship(outfile, child_id, line, rel):
    
    elements = line.split(":")
    parent_id = elements[2].split("!")[0].strip()
    if(rel == Relationships.IS_A):
        relationship = "is_a"
        
    
    else:
        rel_elements = elements[1].strip().split(" ")
        relationship = rel_elements[0]  
        
    outfile.write(child_id + "\t" + parent_id + "\t" + relationship + "\n")

class OBOStates:
    
    ID = 1
    RELATION = 2
    
class Relationships:
    
    IS_A = 1
    PART_OF = 2
    REGULATES = 3
    POSITIVE_REGULATE = 4
    NEGATIVE_REGULATE = 5


class Relationship:
    
    def __init__(self, child, parent):
        self.child_category = child
        self.parent_category = parent


class IsA(Relationship):
    
    def __init__(self, child, parent):
        Relationship.__init__(self, child, parent)
        self.type = Relationships.IS_A
        
        
class PartOf(Relationship):
    
    def __init__(self, child, parent):
        Relationship.__init__(self, child, parent)
        self.type = Relationships.PART_OF


def main():
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--obo", dest="obo_file", required=True)
    parser.add_argument("--output", dest="output_file", default=None)
    
    args = parser.parse_args()
    
    if(args.output_file is None):
        args.output_file = args.obo_file + "_hierarchy.txt"
        
    obo = open(args.obo_file, 'r')
    out = open(args.output_file, 'w')
    
    state = OBOStates.ID    
    child_id = ""
    for line in obo:
        line = line.rstrip("\n")
        #relationships
        #child is a parent (is_a: GO:#####
        #child part of parent (relationship: part_of GO:####
        if("[Term]" in line):
            state = OBOStates.ID
            continue
        if(state == OBOStates.ID):
            if(line.startswith("id:")):
                child_id = line.split(":")[2]
                state = OBOStates.RELATION
        elif(state == OBOStates.RELATION):
            try:
                if(line.startswith("is_a:")):
                    writeRelationship(out, child_id, line, Relationships.IS_A)
                elif(line.startswith("relationship:")):
                    if("part_of GO" in line):
                        writeRelationship(out, child_id, line, Relationships.PART_OF)
                    elif(" regulates GO" in line):
                        writeRelationship(out, child_id, line, Relationships.REGULATES)
                    elif("positively_regulates GO" in line):
                        writeRelationship(out, child_id, line, Relationships.POSITIVE_REGULATE)
                    elif("negatively_regulates GO" in line):
                        writeRelationship(out, child_id, line, Relationships.NEGATIVE_REGULATE)
            except IndexError:
                print("Wonky relationship format, skipping.")
    
    obo.close()
    out.close()
                    
        

    
    
if __name__ == "__main__":
    main()