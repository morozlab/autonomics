from xlrd import open_workbook
import argparse
import sys

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file", dest="inputFile", required=True)

args = parser.parse_args()

wb = open_workbook(args.inputFile)
s = wb.sheets()[0]
for row in range(s.nrows):
	values = []
	for col in range(s.ncols):
		values.append(s.cell(row,col).value)
	sys.write("\t".join(values))
