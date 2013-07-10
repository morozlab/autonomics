import argparse
import os
import sys
import threading
from assembly.consensus.pseudo_contig import *
from assembly.consensus.contig_extension import *
from assembly.consensus.base_weights import *
from alignment.io import BlatAlignment
from Bio import SeqIO

parse = argparse.ArgumentParser(description="This python script merges to FASTA files aligned to each other using BLAT")
parse.add_argument("queryFile", help="Query FASTA file used in the BLAT alignment.")
parse.add_argument("dbFile", help="Database FASTA file used in the BLAT alignment.")
parse.add_argument("alignmentFile", help="BLAT-formatted alignment file.")
parse.add_argument("-t", "--num-threads", dest="numThreads", help="Number of threads to create when creating contigs.")

args = parse.parse_args()

def adjustWeights(positions, index, sequence):
    #adjust the weights of bases for the positions this overlap occupies
    for base in sequence:
        weights = positions[index]
        if(base == 'A'):
            weights.A += 1
        elif(base == 'T'):
            weights.T += 1
        elif(base == 'C'):
            weights.C += 1
        elif(base == 'G'):
            weights.G += 1
        else:
            weights.N += 1
        positions[index] = weights
        index += 1
        
def callConsensus(baseWeights):
    
    if(baseWeights.A == baseWeights.T and baseWeights.T == baseWeights.C and baseWeights.C == baseWeights.G):
        #if impossible to call a base, consensus = backbone, which defaults to N
        baseWeights.consensusBase = baseWeights.backbone
    else:    
        calling = baseWeights.A
        baseWeights.consensusBase = 'A'
        if(baseWeights.T > calling):
            calling = baseWeights.T
            baseWeights.consensusBase = 'T'
        if(baseWeights.C > calling):
            calling = baseWeights.C
            baseWeights.consensusBase = 'C'
        if(baseWeights.G > calling):
            calling = baseWeights.G
            baseWeights.consensusBase = 'G'
        if(baseWeights.N > calling):
            baseWeights.consensusBase = 'N'
    
    sys.stdout.write(baseWeights.consensusBase)    

def parseOverlaps(lines):
    for line in lines:
        elements = line.split("\t")
        reference = elements[13]
        query = elements[9]
        if(not (reference in pseudoContigs)):
            pseudoContigs[reference] = PseudoContig(reference, dbSeqs[reference])

        #get the alignment length and  % identity 
        sizes = elements[18].split(",")
        total = 0
        for size in sizes:
            if(size != ''):
                total += int(size)
    
        #check if this query is already aligned to something
        if(query in bestAlignment):
            #determine if this alignment is better than the current best for the query sequence
            if(bestAlignment[query].length < total):
                #if it is, replace the alignment with this alignment
                bestAlignment[query] = BlatAlignment(query, reference, total, elements[11], elements[12], elements[10], elements[15], elements[16], elements[14], int(elements[0])/total)
    
        else:
            bestAlignment[query] = BlatAlignment(query, reference, total, elements[11], elements[12], elements[10], elements[15], elements[16], elements[14], int(elements[0])/total) 
     
        
def setBackbone(positions, index, sequence):
    adjustWeights(positions, index, sequence)
    for base in sequence:
        weights = positions[index]
        weights.backbone = base
        index += 1
        
    
    

#store the sequences for each file in a dictionary
querySeqs = {}
dbSeqs = {}

for seq_record in SeqIO.parse(args.queryFile, "fasta"):
    querySeqs[seq_record.id] = str(seq_record.seq.upper())
    
for seq_record in SeqIO.parse(args.dbFile, "fasta"):
    dbSeqs[seq_record.id] = str(seq_record.seq.upper())
    
#iterate over the alignment file and assign query sequences to pseudo-contigs
file = open(args.alignmentFile, 'r')

pseudoContigs = {}
bestAlignment = {}
badAlignments = {}
for line in file:
        pass

file.close()

triggering = open("extension_pairs.txt", 'w')
right = open("right_overhangs.txt", 'w')
left = open("left_overhangs.txt", 'w')
#iterate over the alignments, assigning them as extensions to the contigs where possible
for key in bestAlignment.iterkeys():
    alignment = bestAlignment[key]
    
    #temporarily store the pseudoContig for this alignment
    tmp = pseudoContigs[alignment.dbID]
    updated = 0
    leftOverhang = 0
    rightOverhang = 0
    
    #check if alignment has left overhang
    if(alignment.qstart > alignment.dstart):
        leftOverhang = alignment.qstart - alignment.dstart
        left.write(alignment.queryID + "\t" + alignment.dbID + "\t" + str(alignment.qstart) + "\t" + str(alignment.dstart) + "\n")
        #store the extension
        tmp.extensions.append(ContigExtension(leftOverhang, "left", alignment.queryID, querySeqs[alignment.queryID], -1))
        #check if this extension is longer than the longest extension
        if(leftOverhang > tmp.leftOverhangLen):
            tmp.leftOverhangLen = leftOverhang
        
        updated = 1
            
    #check if alignment has right overhang
    if((alignment.dstart + alignment.qlength - alignment.qstart) > alignment.dlength):
        #rightOverhang = length of query - endpoint of query in alignment - distance to end of datbase seq from alignment
        rightOverhang = alignment.qlength - alignment.qend - (alignment.dlength - alignment.dend)
        right.write(alignment.queryID + "\t" + alignment.dbID + "\t" + str(alignment.qlength) + "\t" + str(alignment.qend) + "\t" + str(alignment.dlength) + "\t" + str(alignment.dend) + "\n")
        #store the extension
        tmp.extensions.append(ContigExtension(rightOverhang, "right", alignment.queryID, querySeqs[alignment.queryID], (alignment.dstart - alignment.qstart)))
        #check if this extension is longer than the longest extension
        if(rightOverhang > tmp.rightOverhangLen):
            tmp.rightOverhangLen = rightOverhang

        updated = 1
        
    if(updated == 1):
        triggering.write(alignment.queryID + "\t" + alignment.dbID + "\t" + str(leftOverhang) + "\t" + str(rightOverhang) + "\n")
        pseudoContigs[alignment.dbID] = tmp
    else:
        #check if this alignment is crappy and we shouldn't include it merged with another contig
        if(alignment.length/alignment.qlength < .80):
            badAlignments[alignment.queryID] = 1
        
        
triggering.close()
left.close()
right.close()

 
#iterate over the pseudoContigs, calling consensus sequences
extensionLog = open("extension_log.txt", 'w')   
pseudoIndex = 1
for key in pseudoContigs.iterkeys():
    
    pc = pseudoContigs[key]
    
    #initalize the base weights
    totalLen = pc.leftOverhangLen + pc.baseLen + pc.rightOverhangLen
    index = 0
    while(index < totalLen):
        pc.positions.append(BaseWeights())
        index += 1
    
    #set the weights for the original sequence
    index = pc.leftOverhangLen
    setBackbone(pc.positions, index, pc.baseSeq)
    
    sys.stdout.write(">" + pc.id + "_merged_" + str(pseudoIndex) + "\n")
    
    overhangIndex = 0
    #iterate over the extensions and set the base weights for each position
    for extension in pc.extensions:
        if(extension.direction == "left"):
            overhangIndex = pc.leftOverhangLen - extension.overhang
            
        else:
            overhangIndex = totalLen - (pc.rightOverhangLen - extension.overhang) - len(extension.sequence)
        extensionLog.write(str(overhangIndex) + "\t" + extension.extensionID + "\t" + extension.direction + "\t" + extension.sequence + "\n")
        #adjust the weights of bases for the positions this overlap occupies
        adjustWeights(pc.positions, overhangIndex, extension.sequence)
            
    
    #call the consensus and print the sequence
    for position in pc.positions:
        callConsensus(position)
    sys.stdout.write("\n")
    
    #increment the index and free the memory for this pseudoContig
    pseudoIndex +=1 
    pseudoContigs[key] = PseudoContig('','')
    
extensionLog.close()
        
#print the unused db seqs
for key in dbSeqs.iterkeys():
    if(not (key in pseudoContigs)):
        sys.stdout.write(">" + key + "\n")
        sys.stdout.write(dbSeqs[key] + "\n")

#print the unused query seqs if they're > 5kb
for key in querySeqs.iterkeys():
    if((not (key in bestAlignment)) or (key in badAlignments)):
        if(len(querySeqs[key]) > 2000):
            sys.stdout.write(">" + key + "\n")
            sys.stdout.write(querySeqs[key] + "\n")


   