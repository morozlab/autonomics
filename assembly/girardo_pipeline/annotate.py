from itertools import groupby
from Bio import SeqIO
import re
import sys
import os
from subprocess import call

annotatePath = "/srv/data2/pipeline/annotation/"
def pd(pn): return "/srv/data2/pipeline/" + pn +"/" #Project directory
def blast_out(pn,db): return annotatePath + pn + "_blast_results_" + db + ".txt" #Output file for blast
###################################
def blast_with_db(pn,db):
	call("mpirun -np 64 /usr/local/bin/mpiblast -m 0 -p blastx -d " + db + " -i " + pd(pn) + pn + "_final_assembly.fasta -o " + blast_out(pn,db)) #Run mpiBlast on the fasta file pn with the database db
###################################

def annot8_with_db(pn,db):
	os.chdir(annotatePath)
	call("ln -s " + blast_out(pn,db) + " blast_out/" + pn + "_" + db + "_KEGG.out")
	call("ln -s " + blast_out(pn,db) + " blast_out/" + pn + "_" + db + "_GO.out") #Format blast output for  annot8r
	call("/home/girardo/bin/annot8r.pl " + pn + "_" + db) #Run annot8r
	call("mv output/" + "KEGG.csv" + " " + pd(pn) + pn + "_" + db + "_" + "KEGG.csv")
	call("mv output/" + "GO.csv" + " " + pd(pn) + pn + "_" + db + "_" + "GO.csv") #Move annot8r output to project directory
	call("rm blast_out/*")
	call("rm output/*") #Clean up annot8r results

###################################
def get_ORF_translations(S):
	"""Returns a dictionary of all open reading frames of the given biopython fasta sequence"""
	trans = {}
	for strand, nuc in [(1,S.seq), (-1, S.seq.reverse_complement())]:
		for frame in range(3):
			trans[strand * (frame+1)] = nuc[frame:].translate(to_stop=True)
	return trans
def get_best_frame_translation(hits,read):
	"""Input: contig:frame dictionary for blast hits, contig. Output: 'Best' frame translation"""
	frames = get_ORF_translations(read) #Dictionary of 6 ORFs
	highest = 0
	if hits.has_key(read.id): #If the read has hits...
		return frames[hits[read.id]] #Use the most common hit frame
	for key in frames: #Otherwise, use the longest ORF
		n = len(str(frames[key]))
		if n > highest: 
			acc,highest = key,n
	return frames[acc]
def setup():
	f = SeqIO.parse('tst','fasta') 
	with open("sepiamix_blast_results_nr.txt") as q:
		g = get_blast_frames(q)	
	return (f,g)	
		
			
		
			
def translate(pn):
	"""Translates fasta file using the best frames from nr and swissprot blast results or longest ORF"""
	with open(blast_out(pn,"nr")) as f:
		blasts = get_blast_frames(f)
	with open(blast_out(pn,"swissprot")) as f:
		blasts.update(get_blast_frames(f)) #Get a list of 'known' frames
	reads = SeqIO.parse(pd(pn) + pn + "_final_assembly.fasta",'fasta') #Get the fasta
	def aux_fun(): #turns the NT sequence into an AA sequence
		for seq in reads:
			seq.seq = get_best_frame_translation(blasts,seq)
			yield seq
	return aux_fun()
	

def get_blast_frames(f):
	"""Input:Handle for blast results. Returns a dictionary of (contig name, most common frame for hits)"""
	rx1 = re.compile("Query=")
	rx2 = re.compile(" Frame")
	rx12 = re.compile("(Query=)|( Frame)")
	lines = (l[:-1] for l in f if rx12.match(l)) #Get the lines
	frame_pairs = {}
	for l in  lines:
		name = l[7:]
		frames = []
		for l in lines: #Collect all frame data
			if rx2.match(l): #Associated with the contig
				frames.append(float(l[8:])) 
			else: break
		if frames: #If the contig has hits...
			frame_pairs[name] = most_common(frames) #add (name,best frame)pair to the main list
	return frame_pairs
	
def most_common(L):
	"""Finds most common element in a list"""
	return max(groupby(sorted(L)),key=lambda(x,v):(len(list(v)),-L.index(x)))[0]
	
	
def pfam(pn):
	"""Runs pfam on pn"""
###################################
def main():
	proj_name = sys.argv[0] #project name
	os.chdir(annotatePath)
	blast_with_db(pn,"nr") #
	blast_with_db(pn,"swissprot") #Run blast
	annot8_with_db(pn,"nr") #
	annot8_with_db(pn,"swissprot") #Run annot8r
	with open(pd(proj_name) + proj_name + '_final_assembly_tr.fasta') as f:
		SeqIO.write(translate(pn),f,'fasta') #Translate fasta file
	os.chdir("/srv/data/pfam/PfamScan") #Run Pfam on translation
	subprocess.call('./pfam_scan.pl -fasta ' + pd(pn) + pn + '_final_assembly_tr.fasta -dir /srv/data/pfam -outfile ' + pd(pn) + pn + '_pfam.out')
	subprocess.call('mv ' + blast_out(proj_name,'nr') + ' ' + pd(proj_name)
	subprocess.call('mv ' + blast_out(proj_name,'swissprot') + ' ' + pd(proj_name) #Move blast results to project directory
