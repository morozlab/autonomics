import argparse
import math
import re
import sys
from files.sequences import Filter, Utils
from files.io import SamFile
from alignment.io import Reader
from Bio import SeqIO


def nonZeroAvg(l):
    f = [r for r in l if r != 0]
    return sum(f)/len(f)


def trimLens(lens, minLen):
    '''
    Eliminates any reads below a certain threshold.  Function assumes that input
    list lens is sorted smallest to largest.
    '''

    index = 0
    for i in range(len(lens)):
        if lens[i] < minLen:
            index += 1
        else:
            break

    return lens[index:len(lens)]


def separateByAnnotated(alnFile, alnType, seqFile, seqType):
    p = Reader(alnFile, alnType)
    p.parse()
    annotated = {}
    unannotated = {}
    for key, record in SeqIO.index(seqFile, seqType):
        #add condition for filtering
        if(record.id in p.annotations):
            annotated[key] = record
        else:
            unannotated[key] = record

    return dict(annotated=annotated, unannotated=unannotated)


def getLens(seqs):
    '''
    Parses FASTA file using screed to create a sorted list of contig lengths.
    '''
    lens = []
    for key, sequence in seqs.iteritems():
        lens.append(len(sequence.seq))
    return sorted(lens)


def calcCoverage(sequenceDict, alnFile, alnType, do_quantification=False):

    coverage = {}
    quantification = {}
    maxLen = 0

    if(do_quantification):
        r = Reader(alnFile, alnType)
        r.read()
        reads_mapped = 0
        for record in sequenceDict.values():
            if(record.id in r.annotations):
                seen_hits = set()
                annotations = r.annotations[record.id]
                count = 0
                for annot in annotations:
                    if(annot.reference_id in seen_hits):
                        continue
                    seen_hits.add(annot.reference_id)
                    if(annot.great_alignment()):
                        count += 1
                        reads_mapped += 1

                print("".join([record.id, "\t", str(count)]))

            else:
                print("".join([record.id, "\t0"]))
                #print(str(reads_mapped) + " total mappings (query, reference).")

    else:
        #generate a list for each sequence in the input file with the coverage at each base set to 0
        for record in sequenceDict:
            thisCoverage = [0 for i in range(0, len(record.seq))]
            if(len(record.seq) > maxLen):
                maxLen = len(record.seq)

            #handle some alignment format-specific id quirks
            seqID = record.id
            if(alnType == 'sam'):
                seqID = seqID.split(" ")[0]

            #set the initialized coverage in the dictionary
            coverage[seqID] = thisCoverage

        #iterate over the entries in the alignment file, updating the coverage at each base
        if(alnType == 'sam'):
            fh = open(alnFile, 'r')
            for line in fh:
                #ignore header lines
                if(line.startswith("@")):
                    continue
                elements = line.split("\t")
                if(len(elements) > 2):
                    if(elements[2] != "*"):
                        start = int(elements[3]) - 1
                        end = start + len(elements[9])
                        if(end > len(coverage[elements[2]])):
                            end = len(coverage[elements[2]])
                        for i in range(start, end):
                            coverage[elements[2]][i] += 1


        #print the ids for the columns
        for key in coverage.iterkeys():
            sys.stdout.write(key + "\t")

        sys.stdout.write("\n")

        #print the coverage for each sequence
        for i in range(0, maxLen):
            for key in coverage.iterkeys():
                sys.stdout.write(str(coverage[key][i]) + "\t")
            sys.stdout.write("\n")


def calcNXX(lens, percent):
    '''
    Calculates any NXX (e.g. N50, N90) statistic.
    '''

    lenSum = sum(lens)
    threshold = (float(percent) / 100) * lenSum

    runningSum = 0
    nxx = 0
    nxxLen = 0

    for i in range(len(lens)-1, -1, -1):
        myLen = lens[i]
        nxx += 1
        runningSum += myLen

        if runningSum >= threshold:
            nxxLen = myLen
            break

    return nxx, nxxLen


def calcStats(filename, seqDict, m):

    try:
        minLen = int(m)
    except ValueError:
        print("Minimum contig length (-l) must be an integer to perform stats calculations.")
        return

    print("filename sum n trim_n min med mean max n50 n50_len n90 n90_len var stdev")

    lens = getLens(seqDict)
    trimmedLens = trimLens(lens, minLen)

    if(len(trimmedLens) == 0):
        print(filename + " - no sequences longer than threshold")
        return

    statN = len(lens)
    statTrimmedN = len(trimmedLens)
    statSum = sum(trimmedLens)
    statMin = min(trimmedLens)
    statMax = max(trimmedLens)
    statMed = trimmedLens[(len(trimmedLens)-1)/2]
    statMean = int(statSum / float(statTrimmedN))
    params = calcVarDev(lens, statMean)
    statN50, statN50Len = calcNXX(trimmedLens, 50)
    statN90, statN90Len = calcNXX(trimmedLens, 90)

    print(str (filename) + " " + str(statSum) + " " + str(statN) + " " + \
           str(statTrimmedN) + " " + str(statMin) + " " + str(statMed) + \
           " " + str(statMean) + " " + str(statMax) + " " + str(statN50) \
           + " " + str(statN50Len) + " " + str(statN90) + " " + \
           str(statN90Len) + " " + str(params[0]) + " " + str(params[1]) + "\n")


def calcVarDev(lens, avg):
    variation = 0;
    for length in lens:
        variation += math.pow((length - avg), 2)


    var = variation/len(lens) - 1

    params = [var, math.sqrt(var)]

    return params


def determine_accuracy_sam(genome_file, sam_aln):

    genome_dict = SeqIO.index(genome_file, 'fasta')

    sam_file = SamFile(sam_aln)

    num_bases = 0
    num_errors = 0
    matcher = re.compile("(\d+M)|(\d+I)|(\d+D)")
    num_inserts = 0
    num_dels = 0
    for r in sam_file.iterrecords():
        ref = r.reference
        #check if this read couldn't be mapped
        if(ref == "*" or ((r.flag & 4) != 0)):
            continue

        #get the appropriate version of the reference sequence
        ref_record = genome_dict[ref]
        ref_seq = str(ref_record.seq)

        #set up some initial numbers and counters
        ref_pos = r.position - 1
        query_seq = r.query_seq
        query_pos = 0
        num_bases += len(query_seq)

        #parse the CIGAR string
        cig_elements = matcher.findall(r.cigar)
        for el in cig_elements:
            matches = el[0].replace("M", "")
            inserts = el[1].replace("I", "")
            dels = el[2].replace("D", "")
            if(inserts != ""):
                inserts = int(inserts)
                num_errors += inserts
                query_pos += inserts
                num_inserts += inserts
            elif(dels != ""):
                dels = int(dels)
                num_errors += dels
                ref_pos += dels
                num_dels += dels
            elif(matches != ""):
                matches = int(matches)
                for i in range(0, matches):
                    if(query_seq[query_pos] != ref_seq[ref_pos]):
                        num_errors += 1

                    query_pos += 1
                    ref_pos += 1
                    #stop checking for mismatches if we're at the end of the reference sequence
                    if(ref_pos == len(ref_seq)):
                        break


    return (float(num_errors)/float(num_bases), num_inserts, num_dels, num_bases)



def mapping_accuracy(genome_file, aln_file, aln_type='sam'):

    if(aln_type=='sam'):
        return determine_accuracy_sam(genome_file, aln_file)


def main():

    parser = argparse.ArgumentParser(description = "This script has a bunch of utilities for dealing with sequence files!")

    parser.add_argument('-f', '--file', dest='seqFile', required=True, help="The input file containing the sequences")
    parser.add_argument('-fs', '--multiple-files', dest='multi_files', default=None, nargs="+")
    parser.add_argument('-ft', '--file-type', dest='seqFormat', default='fasta')
    parser.add_argument('-ot', '--output-type', dest='outFormat', default='fasta')
    parser.add_argument('-a', '--alignment-file', dest='alnFile', default=None)
    parser.add_argument('-at', '--alignment-type', dest='alnType', default=None)
    parser.add_argument('-q', '--min-quality-score', dest='minQual')
    parser.add_argument('-p', '--percent', dest='percent', help='Minimum percent for various calculations (percentage of bases having a certain quality score, etc) - enter as whole number, not fraction (e.g. 50, not .50')
    parser.add_argument('--filter-quality', dest='filterQuality', const=True, default=False, action='store_const', help="Filter the sequence file for quality.\nRequires the '-q' and '-p' arguments.")
    parser.add_argument('--annotated', dest='annotated', const=True, default=False, action='store_const', help='Perform actions only on annotated sequences as determined by alnFile')
    parser.add_argument('--combine-files', dest='combine', default=False, const=True, action='store_const')
    parser.add_argument('--get-by-id', dest='getbyid', const=True, default=False, action='store_const', help='Create a new sequence file in the specified output format from sequences matching IDs in the provided list')
    parser.add_argument('--coverage', dest='coverage', const=True, default=False, action='store_const', help="Generate the coverage of this sequence file in some other data set. Requires an alignment file specified with -a")
    parser.add_argument('--histogram', dest='histogram', const=True, default=False, action='store_const', help="Generate a histogram of sequence lengths")
    parser.add_argument('--mapping-accuracy', dest='mapping_accuracy', default=False, const=True, action='store_const', help="Calculate the mapping accuracy of short reads aligned to target sequences. Requires the speficiation of a sam-formatted alignment file with -a and a FASTA-formatted sequence file with -f.")
    parser.add_argument('--num-contigs', dest='num_contigs', const=True, default=False, action='store_const', help='Determine the number of singlets/contigs in the sequence file.')
    parser.add_argument('--quantification', dest='quantification', default=False, const=True, action='store_const')
    parser.add_argument('--remove-duplicates', dest='remove_dups', default=False, const=True, action='store_const')
    parser.add_argument('-l', '--length-cutoff', dest='minLen', default=0, help="")
    parser.add_argument('--sort' , dest='sort', default=None)
    parser.add_argument('--split-paired-end', dest='splitPairedEnd', default=False, const=True, action='store_const', help="Split an interleaved paired-end file into a file for each end.")
    parser.add_argument('--statistics', dest='stats', const=True, default=False, action='store_const', help="Generate sequence statistics for this file")
    parser.add_argument('--unannotated', dest='unannotated', const=True, default=False, action='store_const', help='Perform actions only on unannotated sequences as determined by alnFile')

    args = parser.parse_args()


    if(args.combine):
        for f in args.multi_files:
            fh = open(f, 'r')
            for record in SeqIO.parse(fh, args.seqFormat):
                sys.stdout.write(record.format(args.seqFormat))

            fh.close()

    if(args.coverage):
        if(args.annotated or args.unannotated):
            splitSeqs = separateByAnnotated(args.alnFile, args.alnType, args.seqFile, args.seqFormat)
            if(args.annotated):
                calcCoverage(splitSeqs['annotated'], args.alnFile, args.alnType)
            else:
                calcCoverage(splitSeqs['unannotated'], args.alnFile, args.alnType)

        else:
            #generate the sequence dictionary, call calcCoverage
            calcCoverage(SeqIO.index(args.seqFile, args.seqFormat), args.alnFile, args.alnType)


    if(args.stats):
        '''
        Outputs assembly statistics for provided FASTA files.
        '''
        if(args.annotated or args.unannotated):
            splitSeqs = separateByAnnotated(args.alnFile, args.alnType, args.seqFile, args.seqFormat)
            if(args.annotated):
                calcStats(args.seqFile, splitSeqs['annotated'], args.minLen)
            else:
                calcStats(args.seqFile, splitSeqs['unannotated'], args.minLen)
        elif(args.seqFormat == 'fasta'):
            calcStats(args.seqFile, SeqIO.index(args.seqFile, args.seqFormat), args.minLen)
        elif(args.seqFormat == 'fastq'):
            calcStats(args.seqFile, Utils.getDictFromFile(args.seqFile, args.seqFormat), args.minLen)


    if(args.histogram):
        '''
        This section of code generates histograms of sequence lengths! Well, right now it only generates a text file that can be used to generate histograms in another program. I will fix this, one day
        '''
        #get a size-sorted list of sequence lengths
        if(args.seqFormat == 'fasta'):
            lens = getLens(Utils.getDictFromFile(args.seqFile, args.seqFormat))

        #NEED TO ADD FASTQ support

        trimmedLens = trimLens(lens, int(args.minLen))

        hist = [0 for i in range(0, max(trimmedLens) + 1)]
        for length in lens:
            hist[length] += 1

        for i in range(0, len(hist)):
            sys.stdout.write(str(hist[i]) + "\n")

    if(args.filterQuality):
        '''
        This section of code is responsible for filtering sequences out of a set that fail to meet a quality threshold of minQual for percent bases
        '''
        f = Filter()
        if(args.seqFormat == 'fastq'):
            Utils.printFastqFromDict(f.filterFastqByQuality(args.seqFile, int(args.minQual), float(args.percent)))
        else:
            sys.stderr.write("We only support filtering of FASTQ files currently, check back later!")
            sys.exit()


    if(args.mapping_accuracy):

        print(mapping_accuracy(args.seqFile, args.alnFile))

    if(args.num_contigs):
        matcher = re.compile("contig")

        num_singlets = [0]
        num_contigs = [0]

        def increment_singlets():
            num_singlets[0] = num_singlets[0] + 1

        def increment_contigs():
            num_contigs[0] = num_contigs[0] + 1

        fh = open(args.seqFile, 'r')
        for record in SeqIO.parse(fh, args.seqFormat):
            if(matcher.search(record.id)):
                increment_contigs()
            else:
                increment_singlets()

        fh.close()
        f_singlets = num_singlets[0]
        f_contigs = num_contigs[0]
        print("n num_contigs num_singlets")
        print(" ".join([str(f_singlets + f_contigs), str(f_singlets), str(f_contigs)]))


    if(args.sort):

        sortedSeqs =  Utils.sort(SeqIO.parse(args.seqFile, args.seqFormat), args.sort)

        for record in sortedSeqs:
            if(args.seqFormat == 'fastq'):
                sys.stdout.write("@" + record.id + "\n")
                sys.stdout.write(str(record.seq) + "\n")
                sys.stdout.write("+" + record.id + "\n")
                sys.stdout.write(record.letter_annotations['phred_quality'])
            elif(args.seqFormat == 'fasta'):
                sys.stdout.write(">" + record.id + "\n")
                sys.stdout.write(str(record.seq) + "\n")

    if(args.remove_dups):
        Utils.remove_duplicates(args.seqFile, args.seqFormat)

    if(args.splitPairedEnd):

        firstSeq = None
        end1 = open(args.seqFile + ".1", 'w')
        end2 = open(args.seqFile + ".2", 'w')
        orphan = open(args.seqFile + ".orphans", 'w')
        for record in SeqIO.parse(args.seqFile, args.seqFormat):
            if(firstSeq):
                if(firstSeq.id.split(" ")[0] in record.id):
                    if(not firstSeq.id.endswith("/1")):
                        firstSeq.id += "/1"
                    if(not record.id.endswith("/2")):
                        record.id += "/2"
                    end1.write(firstSeq.format(args.seqFormat) + "\n")
                    end2.write(record.format(args.seqFormat) + "\n")
                    firstSeq = None

                else:
                    orphan.write(firstSeq.format(args.seqFormat) + "\n")
                    firstSeq = record

            else:
                firstSeq = record

        if(firstSeq):
            orphan.write(firstSeq.format(args.seqFormat) + "\n")

        end1.close()
        end2.close()
        orphan.close()

    if(args.quantification):
        if(args.alnFile is None or args.alnType is None):
            parser.print_help()
            sys.exit()

        calcCoverage(Utils.getDictFromFile(args.seqFile, args.seqFormat), args.alnFile, args.alnType, do_quantification=True)


if __name__ == "__main__":
    main()
