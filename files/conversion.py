import sys
from Bio import SeqIO

def tab_to_fasta(tab_file, seq_col, id_cols, sep='\t'):

    print(id_cols)
    fh = open(tab_file, 'r')
    for l in fh:
        l = l.rstrip("\n")
        els = l.split(sep)
        identifier = " ".join([els[id_col] for id_col in id_cols])
        yield "".join([">", identifier, "\n", els[seq_col], "\n"])

    fh.close()
