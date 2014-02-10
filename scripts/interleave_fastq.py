from Bio import SeqIO #Biopython 1.54 or later needed
import argparse
import sys
#######################################################
#
# Change the following settings to suit your needs
#

parser = argparse.ArgumentParser(description = "Script for interleaving fastq paired-end files into a single file")

parser.add_argument('-f1', '--first-file', dest='input_forward_filename', required=True)
parser.add_argument('-f2', '--second-file', dest='input_reverse_filename', required=True)
parser.add_argument('--format', dest='fileFormat', default='fastq')
parser.add_argument('-sf', '--forward-suffix', dest='f_suffix', default=None)
parser.add_argument('-sr', '--reverse-suffix', dest='r_suffix', default=None)
parser.add_argument('-o', '--output-file', dest='outputFile', required=True)

args = parser.parse_args()

input_forward_filename = str(args.input_forward_filename)
input_reverse_filename = str(args.input_reverse_filename)

output_pairs_filename = str(args.outputFile)
output_orphan_filename = str(args.outputFile) + "_unpaired_orphans.fastq"

f_suffix = args.f_suffix
r_suffix = args.r_suffix

#python /home/pwilliams/autonomics//scripts/interleave_fastq.py -f1 /srv/data2/pipeline//Aphrocallistes_vastus_4/Aphrocallistes_vastus_4.fastq -f2 /srv/data2/pipeline//Aphrocallistes_vastus_4/Aphrocallistes_vastus_4.fastq.end2 -o /srv/data2/pipeline//Aphrocallistes_vastus_4/Aphrocallistes_vastus_4.fastq_interleaved.fastq  


#######################################################
if f_suffix:
    f_suffix_crop = -len(f_suffix)
    def f_name(name):
        """Remove the suffix from a forward read name."""
        assert name.endswith(f_suffix), name
        return name[:f_suffix_crop]
else:
    f_name = None

if r_suffix:
    r_suffix_crop = -len(r_suffix)
    def r_name(name):
        """Remove the suffix from a reverse read name."""
        assert name.endswith(r_suffix), name
        return name[:r_suffix_crop]
else:
    r_name = None


print("Indexing forward file...")
forward_dict = SeqIO.index(input_forward_filename, args.fileFormat, key_function=f_name)

print("Indexing reverse file...")
reverse_dict = SeqIO.index(input_reverse_filename, args.fileFormat, key_function=r_name)

print("Ouputing pairs and forward orphans...")
pair_handle = open(output_pairs_filename, "w")
orphan_handle = open(output_orphan_filename, "w")
for key in forward_dict:
    if key in reverse_dict:
        forward_seq = forward_dict.get_raw(key)
        reverse_seq = reverse_dict.get_raw(key)
        pair_handle.write(forward_dict.get_raw(key))
        pair_handle.write(reverse_dict.get_raw(key))
    else:
        orphan_handle.write(forward_dict.get_raw(key))
pair_handle.close()
print("Ouputing reverse orphans...")
for key in reverse_dict:
    if key not in forward_dict:
        orphan_handle.write(reverse_dict.get_raw(key))
orphan_handle.close()
print("Done")
