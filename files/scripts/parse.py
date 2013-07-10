import argparse
import re
import sys

def main():
    parser = argparse.ArgumentParser(description = "This module is used to parse different types of files!")
    parser.add_argument("-csv2tab", dest="csv2Tab", default=False, const=True, action='store_const')
    parser.add_argument("-cd", "--csv-delimiter", dest="csvDelimiter", default=",")
    parser.add_argument("-q", "--text-qualifier", dest="textQualifier", default='"')
    parser.add_argument("inputFile")  
    
    args = parser.parse_args()
    
    if(args.csv2Tab):
        CSVToTabDelimited(args.inputFile, args.textQualifier, args.csvDelimiter)  

def CSVToTabDelimited(fileName, textQualifier, delimiter):
    
    fileToRead = open(fileName, 'r')
    
    for line in fileToRead:
        line = line.strip()
        token = ""
        state = "parsing"
        for char in line:
            if(state == "parsing"):
                if(char == ","):
                    count = token.count('"')
                    if(count >= 2 or count == 0):
                        sys.stdout.write(re.sub("\"", "", token) + "\t")
                        token = ""
                        state = "after_comma"
                    else:
                        token += char
                        
                else:
                    token += char
                
            elif(state == "after_comma"):
                if(char != " "):
                    token += char
                    state = "parsing"
             
                    
        sys.stdout.write(re.sub("\"", "", token) + "\n")
        
main()