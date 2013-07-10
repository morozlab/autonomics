import sys
import os
from subprocess import call
from Bio import SeqIO

pred_dir = '/srv/data2/pipeline/prediction/%s/'
pd = '/srv/data2/pipeline/%s/'
np_dir = '/srv/data2/pipeline/prediction/neuroscripts/'
fa = '%s_final_assembly_tr.fasta' #Fasta protein file
max_contigs = 4000

def format_for_predictors(pn):
	seqs = dict(enumerate(SeqIO.parse(fa % pn,'fasta'))) #Numbers, indexed by protein seqs for recovery after prediction
	t = 0 #File number
	for i in range(len(seqs)):
		if i>= max_contigs * t:
			f = open(pred_dir % pn + pn + '-%s' % t + '.prots','w+') #Split into multiple files (too many in one causes crash in proprietary scripts)
			t+=1
		v = seqs[i]
		v.id = v.description = str(i) #Rename proteins to reduce length (proprietary script truncates long names)
		SeqIO.write(v,f,'fasta')
	return (seqs,t)
######################
#####################
def main():
	pjn = 't' #sys.argv[0] #Project name
	os.chdir(pd % pjn)
	os.mkdir(pred_dir % pjn)
	cypher,t = format_for_predictors(pjn) #Format files and create name mapping
	os.chdir(pred_dir % pjn)
	for flnm in (pjn + '-%s.prots' % n for n in range(t)):
		call(np_dir + "signalp-3.0/signalp -t euk -f short -trunc 70 " + flnm + " > " + flnm.replace('.prots','_signalp.out'), shell=True) #Run SignalP
		call(np_dir + "targetp-1.1/targetp -N " + flnm + " > " + flnm.replace('.prots','_targetp.out'),shell=True) #Run TargetP
		call("cat " + flnm + " | decodeanhmm -f " + np_dir + "tmhmm-2.0c/lib/TMHMM2.0.options -modelfile " + np_dir + "tmhmm-2.0c/lib/TMHMM2.0.model | " + np_dir + "tmhmm-2.0c/bin/tmhmmformat.pl > " + flnm.replace('.prots','_tmhmm.out'),shell=True) #Run TMHMM
		###
		call(np_dir + "phobius/phobius.pl " + flnm + " -short > " + flnm.replace('.prots','_phob.out'),shell=True) #Run Phobius
