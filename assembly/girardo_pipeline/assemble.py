import sys
import os
import subprocess

mira = "/home/girardo/bin/mira"
cap3 = "/home/girardo/bin/cap3"
assembly_path(projName) = "/srv/data2/pipeline/assembly/" + projName + "/"
proj_dir(projName) = "/srv/data2/pipeline/" + projName + "/"
constats_header = "# name              \tlength\tav.qual\t#-reads\tmx.cov.\tav.cov\tGC%\tCnIUPAC\tCnFunny\tCnN\tCnX\tCnGap\tCnNoCov"

pn = sys.argv[1] #set project name
os.mkdir(assembly_path(pn)) #Make the assembly directory
os.chdir(assembly_path(pn)) #And switch to it
subprocess.call("ln -s " + proj_dir(pn) + pn + "_in.solexa.fastq " + assembly_path(pn) + pn + ".fastq") #Create a local link to the fastq file
subprocess.call(mira + " --project=" + pn + " --job=denovo,est,accurate,solexa >&log_assembly.txt") #Run mira on fastq
mira_data = post_process(pn) #Process mira output for cap3 input
subprocess.call(cap3 + " " + mira_data + " " + mira_data + ".qual > " + proj_dir(pn) + "/" + pn + "_cap3.log" #Run cap3 on processed output
subprocess.call("cat " + mira_data + ".cap.contigs " + mira_data + ".cap.singlets > " + proj_dir(pn) + pn + "_final_assembly.fasta"
contstats_path = assembly_path(pn) + pn + "_assembly/" + pn + "_d_info/" + pn + "_info_contigstats.txt"
if not os.path.exists(contstats_path): #If there's no contigstats file, make a dummy. This happens if mira finds no contigs
	with open(contstats_path,'w') as f: f.write(constats_header)
subprocess.call("rm -r " + assembly_path(pn)) #Delete the assembly directory

def post_process(pjn):
	os.chdir(assembly_path(pjn) + pjn + "_assembly/" + pjn + "_d_results")
	with open("../" + pjn + "_d_info/" + pjn + "_info_debrislist.txt","r") as f:
		singlet_ids = f.readlines() #Read the id of singlets
	with open("../../" + pjn + "_in.solexa.fastq","rU") as f:
		singlets = [seq_read for seq_read in SeqIO.parse(f,"fastq") if seq_read.id in singlet_ids] #Get the singlets from the original fastq file
	debris = "seqs_not_in_" + pjn + "_in.solexa.fastq.fas"
	mira_data = pjn + "_out.unpadded.fasta"
	with open(debris,"w") as f: SeqIO.write(singlets,f,"fasta") #Write the singlets in fasta
	with open(mira_data,"a+") as f: SeqIO.write(singlets,f,"fasta") #Append the singlets to the output
	with open(debris + ".qual","w") as f: SeqIO.write(singlets,f,"qual") #Write the singlets in qual
	with open(mira_data + ".qual","a+") as f: SeqIO.write(singlets,f,"qual") #Append the singlets in qual to the output
	return mira_data
