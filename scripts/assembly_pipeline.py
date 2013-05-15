import argparse
import shutil
import subprocess
import os
from Bio import SeqIO
from sqlalchemy.sql import and_
from zeroclick import netutils
from zeroclick.settings import SCRIPTPATH
from zeroclick.utility import replace_nonstandard_bases

#collect arguments for the assembly pipeline
parser = argparse.ArgumentParser(description='Submission script for two-stage transcriptome assembly pipeline.')
parser.add_argument('-nml', '--newbler_minimum_length', default='40', dest='newblerMinimumLength', help='Minimum length of raw reads to be included in assembly input.')
parser.add_argument('-nmi', '--newbler_minimum_identity', default="90", dest='newblerMinimumIdentity')
parser.add_argument('-nmo', '--newbler_minimum_overlap', default="40", dest='newblerMinimumOverlap')
parser.add_argument('-margs', '--mira_argument_list', default='', dest='miraArgumentList')
parser.add_argument('-cpu', '--num_processors', default="1", dest='numCPU')
parser.add_argument('-o', '--output_prefix', dest='outputPrefix', default='transcriptome_assembly', help='Prefix for all assembly output files')
parser.add_argument('-a', '--adapter-file', dest='adapters', default='none', help='Optional file of adapters to trim off')
parser.add_argument('-v', '--vector-file', dest='contaminants', default='none', help='Optional file of contaminants to screen')
#parser.add_argument('-b', '--barcode-file', dest='barcodeFile', help='File containing multiplexing barcodes.')
parser.add_argument('-mira', '--mira-only', dest='miraOnly', default=False, const=True, action='store_const')
parser.add_argument('-nq', '--no-quality', dest='noQual', default=False, const=True, action='store_const')
parser.add_argument('fileType', help='Type of the files being used for this assembly. Currently support FASTQ and SFF files.')
parser.add_argument('directory', help='Directory containing the files for this assembly.')
args = parser.parse_args()
#if args.pn == 'none': args.pn = [x for x in args.directory.split('/') if x!=''][-1]

#function definitions
def getFileList(extension, dir):
    fileList = [os.path.normcase(f) for f in os.listdir(dir)]
    fileList = [os.path.join(args.directory, f) for f in fileList if(os.path.splitext(f)[1] == extension)]
    return fileList

def main():
    if(args.miraOnly):
        assemble(args.directory, "")
    else:
        #pre-process the input files
        if(args.fileType == 'fastq'):
            assemblyFile = args.directory + "/" + args.outputPrefix + "_input_reads.fastq"

            os.chdir(args.directory)
            if(os.path.isfile(assemblyFile)):
                os.remove(assemblyFile)

            fileList = [os.path.normcase(f) for f in os.listdir(args.directory)]
            fileList = [os.path.join(args.directory, f) for f in fileList if(os.path.splitext(f)[1] == '.fastq')]

            combinedFile = open(assemblyFile, 'w')
            cat(fileList, combinedFile)
            combinedFile.close()

            #process barcodes, if they exist
            #if(args.barcodeFile):
            #    handleBarcodes(args.barcodeFile, args.fileType, assemblyFile)

#            else:
                #start the assembly
            assemble(args.directory, assemblyFile)

        elif(args.fileType == 'fasta'):
            assemblyFile = args.directory + "/" + args.outputPrefix + "_input_reads.fasta"
            qualFile = args.directory + "/" + args.outputPrefix + "_input_reads.qual"
            os.chdir(args.directory)
            if(os.path.isfile(assemblyFile)):
                os.remove(assemblyFile)

            #combine sequence file
            fileList = getFileList(args.directory, '.fasta')
            combinedFile = open(assemblyFile, 'w')
            cat(fileList, combinedFile)
            combinedFile.close()

            #combine quality file
            fileList = getFileList(args.directory, '.qual')
            combinedFile = open(qualFile, 'w')
            cat(fileList, combinedFile)
            combinedFile.close()

            assemble(args.directory, assemblyFile, qualFile)


        elif(args.fileType == 'sff'):
            #if this directory contains more than 1 sff file, merge them together
            assemblyFile = args.directory + "/" + args.outputPrefix + "_input_reads.sff"

            os.chdir(args.directory)

            if(os.path.isfile(assemblyFile)):
                os.remove(assemblyFile)

            #get the pertinent files for assembly in the specified directory
            fileList = [os.path.normcase(f) for f in os.listdir(args.directory)]
            fileList = [os.path.join(args.directory, f) for f in fileList if os.path.splitext(f)[1] == '.sff']

            sffList = []
            if(len(fileList)):

                for f in fileList:
                    sffList.append(f)

                #create the combined SFF file for the assembly
                argList = ['/opt/454/bin/sfffile', '-o', assemblyFile]
                argList[len(argList):] = sffList
                process = subprocess.Popen(argList)
                process.wait()

                assemble(args.directory, assemblyFile)


#basically the cat shell command
def cat(fileList, outFile):
    fileList = list(fileList) #Just in case
    last = len(fileList) - 1
    for i,f in enumerate(fileList):
        readFile = open(f, 'r')
        line = readFile.readline()
        while line:
            outFile.write(line)
            line = readFile.readline()
        readFile.close()
        if i!=last: outFile.write("\n")


def assemble(directory, inputFile, qualFile = None):

    #combine the assembly products with the unused reads
    miraAssemblyID = args.outputPrefix

    if(not args.miraOnly):
        os.chdir(directory)
        newblerCmd = ["/opt/454/bin/runAssembly", "-cdna", "-minlen", args.newblerMinimumLength,"-cpu", args.numCPU, "-m", "-ml", args.newblerMinimumOverlap, "-mi", args.newblerMinimumIdentity, "-o", directory + "/NewblerAssembly"]

        if(args.contaminants != 'none'):
            newblerCmd.extend(["-vs ", args.contaminants])
        if(args.adapters != 'none'):
            newblerCmd.extend(["-v", args.adapters])

        newblerCmd.extend(["-force", inputFile])
        #start the newbler assembly on the merged sff file
        process = subprocess.Popen(newblerCmd)
        process.wait()

        #parse the output of the newbler assembly
        #first, get the set of sequences that weren't assembled
        readStatus = open(directory + "/NewblerAssembly/454ReadStatus.txt", 'r')
        line = readStatus.readline()
        listPath = directory + "/unassembledReads.list"
        listFile = open(listPath, 'w')
        while line:
            split = line.split()
            if(split[1] == "Singleton"):
                listFile.write(split[0] + "\n")
            line = readStatus.readline()

        listFile.close()


        #convert the sff file to a fasta file
        if(args.fileType == 'sff'):
            process = subprocess.Popen(["python", SCRIPTPATH + "sff_extract.py", "-o", directory + "/tempfasta", inputFile])
            process.wait()

        #convert the fastq file to a fasta file
        elif(args.fileType == 'fastq'):
            process = subprocess.Popen("perl " + SCRIPTPATH + "fastqToFasta.pl " + inputFile + " " + directory + "/tempfasta", shell=True)
            process.wait()

        elif(args.fileType == 'fasta'):
            process = subprocess.Popen("mv " + inputFile + " " + directory + "/tempfasta.fasta", shell=True)
            process.wait()
            process = subprocess.Popen("mv " + qualFile + " " + directory + "/tempfasta.fasta.qual", shell=True)
            process.wait()


        tempSeqFile = directory + "/tempfasta.fasta"
        tempQualFile = directory + "/tempfasta.fasta.qual"

        #get the unassembled reads and qualities from the converted FASTA file
        process = subprocess.Popen("perl " + SCRIPTPATH + "getSequencesByID.pl " + tempSeqFile + " " + listPath + " fasta > " + directory + "/unusedReads.fa", shell=True)
        process.wait()
        process = subprocess.Popen("perl " + SCRIPTPATH + "getQualityByID.pl " + tempQualFile + " " + listPath + " fasta > "+ directory + "/unusedReads.qual", shell=True)
        process.wait()

        os.remove(listPath)
        os.remove(tempSeqFile)
        os.remove(tempQualFile)
        if(os.path.isfile(directory + "/tempfasta.xml")):
            os.remove(directory + "/tempfasta.xml")

        #combine the assembly products with the unused reads
        miraAssemblyID = args.outputPrefix

        #merge the resulting unused reads and quality with the contig reads and quality
        seqFiles = [directory + "/unusedReads.fa", directory + "/NewblerAssembly/454AllContigs.fna"]
        qualFiles = [directory + "/unusedReads.qual", directory + "/NewblerAssembly/454AllContigs.qual"]

        newFile = open(directory + "/" + miraAssemblyID + "_in.454.fasta", 'w')
        cat(seqFiles, newFile)
        newFile.close()

        newFile = open(directory + "/" + miraAssemblyID + "_in.454.fasta.qual", 'w')
        cat(qualFiles, newFile)
        newFile.close()

        os.remove(directory + "/unusedReads.fa")
        os.remove(directory + "/unusedReads.qual")

    #submit newbler output to MIRA
    log_file = args.directory + "/" + args.outputPrefix + "_mira.log"
    jobString = " --job=denovo,est,normal,454"
    if(args.noQual):
        jobString += " --noqualities"
        clString = ""

    fastq_in_path = args.directory + "/" + args.outputPrefix + ".fastq"
    if(args.miraArgumentList == ''):
        manifestFile = args.directory + "/" + "mira_manifest"
        with open(manifestFile,'w+') as f:
            f.write("project = " + args.outputPrefix +"\n")
            f.write("job = est,denovo,accurate\n")
            #mrpc - minimum reads per contig, sssip - save simple singlets in project
            f.write("parameters = -GE:not=" + args.numCPU + " IONTOR_SETTINGS -AS:mrl=30:mrpc=1 -OUT:sssip=yes")
            f.write('\n')
            f.write("readgroup = " + args.outputPrefix + '\n')
            f.write("data = " + fastq_in_path + "\n")# ".khmer.out\n")
            f.write("technology = iontor")

       # miraArgs = ["mira", "--project=" + miraAssemblyID, jobString,  "--notraceinfo", "-GE:not=" + args.numCPU,  "454_SETTINGS", "-SK:pr=70", "-AL:mo=40:ms=10:mrs=70", "-AS:mrl=40", clString, "-OUT:sssip=yes"]
        os.chdir(args.directory)
        process = subprocess.Popen("mira " + manifestFile + " > " + log_file, shell=True)

    else:
        raise Exception("custom arguments not supported for mira")
    #        process = subprocess.Popen("mira " + args.miraArgumentList, shell=True)

    process.wait()

    ################# Post-proessing of MIRA assembly files
    miraDir = directory + '/' + miraAssemblyID + "_assembly/"
    unusedReadsFile = directory + '/' + miraAssemblyID + "_unused_reads.fasta"
    debrisList = miraDir + miraAssemblyID + "_d_info/" + miraAssemblyID +"_info_debrislist.txt"
    miraContigs = miraDir + miraAssemblyID + "_d_results/" + miraAssemblyID +"_out.unpadded.fasta"

    #get the unused reads from MIRA
    #unused_ids = set()
    #debris = open(debrisList, 'r')
    #for line in debris:
    #    line = line.rstrip()
    #    unused_ids.add(line)

    #debris.close()

    ### STOPPED HERE #####
    #seq_fh = open()

    #generate a file of MIRA singlets
    #process = subprocess.Popen("perl " + SCRIPTPATH + "getSequencesByID.pl " + directory + "/" + miraAssemblyID + "_in.454.fasta "+ debrisList + " fasta > " + unusedReadsFile, shell=True)
    #process.wait()

    #recover NEWBLER contigs from the MIRA singlets file
    if not(args.miraOnly):
        process = subprocess.Popen("python " + SCRIPTPATH + "ContigRecovery.py " + debrisList + " " + unusedReadsFile, shell=True)
        process.wait()

    #generate a file of contigs for the final assembly
    contigs = open(directory + "/" + miraAssemblyID + "_final_assem_contigs.fa", 'w')
    if not(args.miraOnly):
        contigFiles = [miraDir +"recoveredContigs.fa", miraContigs]
    else:
        contigFiles = [miraContigs]
    cat(contigFiles, contigs)
    contigs.close()

    #combine the unused MIRA input seqs and the unpadded fasta
    assem_file = directory + "/" + miraAssemblyID + "_project.fasta"
    final = open(assem_file, 'w')
    if not(args.miraOnly):
        finalFiles = [unusedReadsFile, miraContigs]
    else:
        finalFiles = [miraContigs]

    print(finalFiles)
    print(final)
    cat(finalFiles, final)
    final.close()

    #convert all non-standard bases to gaps for downstream processing (translation breaks if this step is skipped)
    tmp_file = directory + "/replace_nonstandard_bases.fasta"
    tmp = open(tmp_file, 'w')
    assem = open(assem_file, 'r')
    for record in SeqIO.parse(assem, "fasta"):
        tmp.write("".join([">", record.id, "\n"]))
        tmp.write("".join([replace_nonstandard_bases(str(record.seq)), "\n"]))

    tmp.close()
    assem.close()

    shutil.copy(tmp_file, assem_file)
    os.remove(tmp_file)

    #create the quantification file from the contigstats file
    contig_stats = miraDir + miraAssemblyID + "_d_info/" + miraAssemblyID + "_info_contigstats.txt"
    final_quant = directory + "/" + miraAssemblyID + "_quantification.txt"
    in_file = open(contig_stats, 'r')
    out_file = open(final_quant, 'w')
    for line in in_file:
        if("# name" in line):
            continue
        line = line.rstrip("\n")
        el = line.split("\t")
        out_file.write(el[0] + "\t" + str(el[3]) + "\n")

    in_file.close()
    out_file.close()

    #remove the log file
    #os.remove(log_file)
    #if(miraDir != "/"):
        #
        #shutil.rmtree(miraDir)


######################## END ASSEMBLY SUBROUTINE #########################################

#def handleBarcodes(barcodeFile, fileType, assemblyFile):
#
    #generate a list of the barcodes
#    bcFile = open(barcodeFile, 'r')
#    barcodes = []
#    for line in bcFile:
#        line = line.strip()
#        elements = line.split("\t")
#        if(len(elements) > 1):
#            barcodes.append(Barcode(elements[0], elements[1], elements[2]))
#    bcFile.close()

    #append an entry in barcodes for the unmatched sequences
#    barcodes.append(Barcode('unmatched', 'XXXXXX', 'unmatched'))

#    if(fileType == 'fastq'):
        #use the fastx toolkit to parse the barcoded file, with the original path as the prefix and the extension as the suffix
#        process = subprocess.Popen("cat " + assemblyFile + " | /usr/local/bin/fastx_barcode_splitter.pl --bcfile \"" + args.directory + "/barcodes.txt\" --bol --prefix \""
#                                + os.path.splitext(assemblyFile)[0] + "_\" --suffix \"" + os.path.splitext(assemblyFile)[1] + "\"", shell=True)
#        process.wait()

        #make directories for the individual assemblies, and start the assembly
#        for barcode in barcodes:
#            barcodeDir = args.directory + "/" + barcode.sample + "_assembly"
#            if(not (os.path.isdir(barcodeDir))):
#                os.mkdir(barcodeDir)
#            thisFile = os.path.splitext(assemblyFile)[0] + "_" + barcode.label + os.path.splitext(assemblyFile)[1]
            #trim the barcode off the front of the sequence
#            if(barcode.label != "unmatched"):
#                process = subprocess.Popen("fastx_clipper -Q33 -a " + barcode.sequence + " -i " + thisFile + " -o " + args.directory + "/clipped_tmp.fastq", shell=True)
#                process.wait()
#                process = subprocess.Popen("mv " + args.directory + "/clipped_tmp.fastq " + thisFile, shell=True)
#                process.wait()

#            assemble(barcodeDir, thisFile)


if __name__ == "__main__":
    main()
