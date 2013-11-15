#! /usr/freeware/bin/gawk -f

# print.awk
# 
# This script prints the contents of the final SecretomeP 1.0 output table.
# The input is the result of pasting of 'secretome.howlin' (just seqid),
# 'secretome.nn' (average NN score) and 'pretty.awk' (odds etc.). The input
# lines (one per sequence) are:
# 
# <seqid> <NN-score> "Secreted" ["=>"] <odds> <post.prob>
# 
# The script also catches SignalP output for the sequences.
# 
# VERSION:	Sep  9, 2004	launch, K. Rapacki
# 


BEGIN { # read 'signalp' output file and store the answers in array
	while (getline <"signalp.out") {
	      if ($2=="Y")
	         sp[$1]="signal peptide predicted by SignalP";
	      else if ($2=="0")
	         sp[$1]="SignalP not run";
	      else
	         sp[$1]="-";
	}
}

{
  printf("%-15s\t%5.3f\t%5.3f\t  %5.3f\t   %s\n", \
                $1,$2,$NF,$(NF-1),sp[$1]);
}
