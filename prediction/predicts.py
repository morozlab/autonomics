from string import replace
import argparse
import os
from copy import copy
import subprocess
from subprocess import call
from Bio import SeqIO
# from autonomics import settings
import time
import sys

parser = argparse.ArgumentParser(description='Signaling Peptide prediction script for the zero-click pipeline')
parser.add_argument('projName', help='Unique identifier (name) for this project. Should be consistantly used across the pipeline.')
parser.add_argument('--sp_cutoff', dest="sp_cutoff", default = 0.5)
parser.add_argument('--max_threads',type=int, dest='max_threads', default=1)

parser.add_argument('--use_tmhmm',type=int,help='0 or 1',default=1)
parser.add_argument('--use_phobtm',type=int,help='0 or 1',default=1)
parser.add_argument('--use_phobsp',type=int,help='0 or 1',default=1)

args = parser.parse_args()
pn = args.projName

############################################
processor = 'x86_64'

#ppath = settings.proj_dir
ppath = os.environ["PROJECT_PATH"]
apath = os.environ["AUTONOMICS_PATH"]

pred_dir = ppath + '/' + pn + '/prediction/'
pd = ppath + '/' + pn
np_dir = apath + '/prediction/neuroscripts/' 
fa = '%s_proteins.fasta' % pn
# outname = '_peptide_prediction_' + str(args.sp_cutoff)
outname = '_peptide_prediction'
fa_outname = '_predicted_' + str(args.sp_cutoff) + '.fasta'

############################################

start = time.time()

# print "outname: ", outname
# print "pd: ", pd
# print "fa: ", fa
# print "pred_dir: ", pred_dir
#print "np_dir: ", np_dir
print "\n================================================================================================================================="
print "using: ", args.max_threads, " threads; and a confidence score threshold of: ", str(args.sp_cutoff), " for: ", pn
print "=================================================================================================================================\n"
sys.stdout.flush()

max_contigs = 2000 #4000 is too high, deterimined experimentally. Mysterious HOW errors caused if too high
def format_for_predictors():
    """Format protein file for prediction scripts. Return name-number mapping and total sequences"""
    seqs = dict(enumerate(SeqIO.parse(fa,'fasta'))) #protein seqs indexed by numbers for recovery after prediction
    t = 0 #File number
    for i in range(len(seqs)):
        if i >= max_contigs * t:
            flnm = pred_dir + pn + '-%s.prots' % t #Split into multiple files (too many seqs in one file causes crash in signalp/targetp)
            t+=1
        v = copy(seqs[i])
        v.name = str(i) #Rename proteins to reduce length (proprietary script truncates long names)
        write_fasta_sameline(flnm,v,'a+')
    return (seqs,t)
def parse_predictions(n):
    """parses prediction output files and combines them into a dictionary indexed by contig name"""
    with open(pn + '-%s_signalp4.out' % n,'r') as f: #Open signalp v4 results
        file = [l.split() for l in f.readlines()[2:]]
        results = dict([(l[0],[float(l[8])]) for l in file]) #prediction of secretory signal peptide, initialize dictionary
    with open(pn + '-%s_tmhmm.out' % n,'r') as f: #Open TMHMM results
        file = [l.split() for l in f if 'Number' in l]
    for l in file:
        results[l[1]].append(int(l[-1])) #Accumulate TMHMM results for each name
    with open(pn + '-%s_phob.out' % n,'r') as f: #Open Phobius results
        file = [l.split() for l in list(f)[1:]]
    for l in file:
        results[l[0]].extend([int(l[1]),int(l[2]=='Y')]) #Accumulate Phobius results for each name
    return results
######################
#####################
def main():
    os.chdir(pd)
    os.mkdir(pred_dir)
    cypher,t = format_for_predictors() #Format files and create name mapping
    os.chdir(pred_dir)
    threads = []
    def block(): #Wait until threads are free
        while True:
            running = sum(x.poll()==None for x in threads)
            if running >= args.max_threads:
                time.sleep(10)
            else: break

    for flnm in (pn + '-%s.prots' % n for n in range(t)):
        block()
        threads.append(subprocess.Popen(np_dir + "signalp-4.0/signalp -t euk -f short " + flnm + " > " + flnm.replace('.prots','_signalp4.out'), shell=True)) #Run SignalP v4
        print "cat " + flnm + " | " + np_dir + "tmhmm-2.0c/bin/decodeanhmm.Linux_" + processor + " -f " + np_dir + "tmhmm-2.0c/lib/TMHMM2.0.options -modelfile " + np_dir + "tmhmm-2.0c/lib/TMHMM2.0.model | " + np_dir + "tmhmm-2.0c/bin/tmhmmformat.pl > " + flnm.replace('.prots','_tmhmm.out')
        sys.stdout.flush()
        block()
        threads.append(subprocess.Popen("cat " + flnm + " | " + np_dir + "tmhmm-2.0c/bin/decodeanhmm.Linux_" + processor + " -f " + np_dir + "tmhmm-2.0c/lib/TMHMM2.0.options -modelfile " + np_dir + "tmhmm-2.0c/lib/TMHMM2.0.model | " + np_dir + "tmhmm-2.0c/bin/tmhmmformat.pl > " + flnm.replace('.prots','_tmhmm.out'),shell=True)) #Run TMHMM
        print "Running PHOBIUS: " + flnm
        sys.stdout.flush()
        block()
        threads.append(subprocess.Popen(np_dir + "phobius/phobius.pl " + flnm + " -short > " + flnm.replace('.prots','_phob.out'),shell=True)) #Run Phobius
    for thrd in threads: thrd.wait() #Wait for prediction processes to finish
    with open(pd + '/' + pn + outname,'w+') as h:
        h.write('name,sp4,tmhmm,phobtm,phobsp,is_sp  sp_cutoff: ' + args.sp_cutoff + '\n')
        sp = 0
        tmhmm = 1
        phobtm = 2
        phobsp = 3
        for i in range(t):
            results = parse_predictions(i)
            for n in results: #Write compiled prediction results, restoring original name
                preds = results[n]
                isSP = float(preds[sp])>float(args.sp_cutoff) and \
                ((not args.use_tmhmm) or int(preds[tmhmm])==0 or \
                 (args.use_phobtm or args.use_phobsp) and \
                 ((not args.use_tmhmm) or int(preds[tmhmm])>0) and \
                 ((not args.use_phobtm) or int(preds[phobtm])==0) and  \
                 ((not args.use_phobsp) or int(preds[phobsp])==1) \
                 )


                preds.append(isSP and '+' or '-')
                h.write(cypher[int(n)].name + ',' + ','.join([str(x) for x in preds]) + '\n') #Write to file, and restore original name
                if isSP: write_fasta_sameline(pd + '/' + pn + fa_outname,cypher[int(n)],'a+')
    call(['rm','-r',pred_dir]) #Clean up the prediction directory

    elapsed = (time.time() - start) / 60.0
    print "\n================================================================================================================================="
    print "Total time to run peptide prediction for ", pn, ": ", elapsed, " minutes"
    print "=================================================================================================================================\n"
    sys.stdout.flush()

####################################################
def write_fasta_sameline(flnm,seq,mode='w+'):
    with open(flnm,mode) as f:
        f.write('>' + seq.name)
        f.write('\n')
        f.write(replace(replace(str(seq.seq),'J','X'),'B','X'))
        f.write('\n')
####################
main()
