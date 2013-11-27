from Bio import SeqIO
import argparse
import os

parser = argparse.ArgumentParser(description='parse output from prediction pipeline and produce fasta of predicted peptides')
parser.add_argument('projName', help='Unique identifier (name) for this project. Should be consistantly used across the pipeline.')
parser.add_argument('--sp_cutoff',type=float, default = 0.5)
parser.add_argument('--use_tmhmm',type=int,default=1)
parser.add_argument('--use_phobtm',type=int,default=1)
parser.add_argument('--use_phobsp',type=int,default=1)
parser.add_argument('--default',action='store_true',help='use the original +/- output from the pipeline',default=False)

args = parser.parse_args()
sp = 1
tmhmm = 2
phobtm = 3
phobsp = 4

pn = args.projName

ppath = os.environ["PIPEPATH"]
apath = os.environ["AUTONOMICS_PATH"]

pd = ppath + '/' + pn
fa = '%s_proteins.fasta' % pn
proteins = pd + "/" + fa
outname = '_peptide_prediction'
print "proteins: ", proteins

outfile =  ppath + '/' + pn + '/' + pn + "_peptides_" + str(args.sp_cutoff) + '.fasta'
pred = pd + '/' + pn + outname
print "pred: ", pred
print "outfile: ", outfile

with open(pred,'r') as f:
    f.next()
    secreted = set(preds[0] for preds in (l.strip().split(',') for l in f) if \
                float(preds[sp])>float(args.sp_cutoff) and \
                ((not args.use_tmhmm) or int(preds[tmhmm])==0 or \
                 (args.use_phobtm or args.use_phobsp) and \
                 ((not args.use_tmhmm) or int(preds[tmhmm])>0) and \
                 ((not args.use_phobtm) or int(preds[phobtm])==0) and  \
                 ((not args.use_phobsp) or int(preds[phobsp])==1) \
                 ))

print str(len(secreted)) + ' secreted peptides predicted'

with open(outfile,'w+') as h:
    for p in SeqIO.parse(proteins,'fasta'):
        if p.name in secreted:
            h.write('>' + p.name)
            h.write('\n')
            h.write(str(p.seq))
            h.write('\n')
    
    
    
