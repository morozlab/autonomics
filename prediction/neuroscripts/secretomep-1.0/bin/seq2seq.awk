#! /usr/freeware/bin/gawk -f


# Command line parameters:
# 	OUTDIR		output directory
# 	OUTPUTFORMAT	output format
# 	MAXNAMELEN      max name length
# 	N		flag: replace non-nucleotide symbols with this char
# 	P		flag: replace non-aa symbols with this char
# 	V		verbosity level
# 	X		do not tolerate non-standard symbols?
# 
# Functions:
# 	consider()
# 	reset()
# 	report()
# 	printhow(OUTDIR,len,name,seq,ass)
# 	printfasta(OUTDIR,len,name,seq) 
# 

BEGIN {

  # allowed sequence symbols
  NON_NUCL_ALPH	= "[^ACGTU]";
  NON_PROT_ALPH	= "[^ACDEFGHIKLMNPQRSTVWY]";

  # prepare output directory
  if (OUTDIR!="/dev/stdout") {
     system("rm -rf " OUTDIR);
     system("mkdir " OUTDIR);
  }

  # prevent printing of non-existent previous entry
  corruption = 1;

}

# HOW entry header ============================================================
/^[ 1-9]?[ 0-9][ 0-9][ 0-9][ 0-9][0-9] [^ \t]/ {

  consider();			# print previous entry if applicable

  reset();			# reset flags and entry data

  eformat="how";
  ec[eformat]++;		# increment input entry count

  len = $1;			# sequence length
  name = $2;			# sequence name

  i = 0;			# running index

  ishowseq = 1;			# flag: picking up HOW sequence

  next;

}

# FASTA entry header ==========================================================
/^>/ {

  consider();			# print previous entry if applicable

  reset();			# reset flags and entry data

  eformat="fasta";
  ec[eformat]++;		# increment entry count

  len = 0;			# sequence length (not specified in FASTA)
  name = substr($1,2);		# sequence name

  if (name=="")			# special case: nameless FASTA entry
     name = ("seq." ec["how"]+ec["fasta"]);

  i = 0;			# running index

  isfastaseq = 1;		# flag: picking up FASTA sequence

  next;

}

###############################################################################
corruption { next; }
###############################################################################

# picking up HOW sequence =====================================================
ishowseq {

  # get and check line length .................................................
  llen = length($1);		# sequence line length (not more than 80)
  if (llen>80) {
     if (V)
        print "seq2seq: HOW sequence line too long in " name "\"" \
		>"/dev/stderr";
     corruption++;
     next;
  }

  # convert to upper-case and replace non-standard symbols ....................
  $1 = toupper($1);
  if (N!="")
     nonstd += gsub(NON_NUCL_ALPH,N,$1);
  else if (P!="")
     nonstd += gsub(NON_PROT_ALPH,P,$1);

  # load sequence .............................................................
  for (j=1; j<=llen; j++)
      seq[++i] = substr($1,j,1);

  # check total sequence length so far ........................................
  if (i>len) {
     if (V)
        print "seq2seq: HOW sequence too long in \"" name "\"" \
		>"/dev/stderr";
     corruption++;
     next;
  }
  else if (i==len) {
     ishowseq = 0;
     ishowass = 1;
     seqlen = i;
     i = 0;
     next;
  }
  
}

# picking up HOW assignment ===================================================
ishowass {

  # get and check line length .................................................
  llen = length($1);		# assignment line length (not more than 80)
  if (llen>80) {
     if (V)
        print "seq2seq: HOW assignment line too long in \"" name "\"" \
		>"/dev/stderr";
     corruption++;
     next;
  }
  
  # load assignment ...........................................................
  for (j=1; j<=llen; j++)
      ass[++i] = substr($1,j,1);

  # check total assignment length so far ......................................
  if (i>len) {
     if (V)
        print "seq2seq: HOW assignment too long in \"" name "\"" \
		>"/dev/stderr";
     corruption++;
     next;
  }

  seqass = i;
  
}

# picking up FASTA sequence ===================================================
isfastaseq {

  # remove whitespace and get line length .....................................
  gsub("[ \t]","",$0); 
  llen = length($0);

  # convert to upper-case and replace non-standard symbols ....................
  $0 = toupper($0);
  if (N!="")
     nonstd += gsub(NON_NUCL_ALPH,N,$1);
  else if (P!="")
     nonstd += gsub(NON_PROT_ALPH,P,$1);

  # load sequence .............................................................
  for (j=1; j<=llen; j++)
      seq[++i] = substr($0,j,1);

  len=i;				# keep it updated ...
  seqlen = i;

}

# =============================================================================
END {	consider(); 
        if (V>1)
	   report();
}


###############################################################################
function consider() {
###############################################################################

  if (corruption) {			# count corrupted entries
     cec[eformat]++;
     return 0;
  }
  else if (seqlen==0) {			# count zero length entries
     if (V)
        print "seq2seq: entry \"" name "\" has zero length" \
		>"/dev/stderr";
     zec[eformat]++;
     return 0;
  }
  else if (T && (seqlen>T)) {		# count too long entries
     if (V)
        print "seq2seq: entry \"" name "\" is too long" \
		>"/dev/stderr";
     Tec[eformat]++;
     return 0;
  }
  else if (seqlen<B) {			# count too short entries
     if (V)
        print "seq2seq: entry \"" name "\" is too short" \
		>"/dev/stderr";
     Bec[eformat]++;
     return 0;
  }

  # fix name ..................................................................
  origname = name;			# special characters
  namesub = gsub("[^A-Za-z0-9+,-._]","_",name);
  if (namesub) {
     if (V)
        print "seq2seq: entry name altered from \"" origname "\" to \"" \
		name "\"" >"/dev/stderr";
     mec[eformat]++;
  }

  if (MAXNAMELEN) {			# length
     if (length(name)>MAXNAMELEN) {
        name = substr(name,1,MAXNAMELEN);
        if (V)
           print "seq2seq: entry name truncated from \"" origname "\" to \"" \
		name "\"" >"/dev/stderr";
        tec[eformat]++;
     }
  }

  if (name in present) {		# non-unique name
     if (V)
        print "seq2seq: occurrence",present[name]+1,"of \"" name "\"," \
		" discarded"  >"/dev/stderr";
     dec[eformat]++;
     present[name]++;
     return 0;
  }
  present[name]=1;

  # end of name fixing ........................................................

  if (nonstd) {				# count nonstd symbol entries
     xec[eformat]++;
     if (V) {
        if (N!="")
           print "seq2seq: non-standard nucleotide symbols in \"" name "\"" \
	   	>"/dev/stderr";
        else if (P!="")
           print "seq2seq: non-standard protein symbols in \"" name "\"" \
	   	>"/dev/stderr";
     }
     if (X)
        return 0;
  }

  if (eformat=="how") {			# HOW: seq and ass length check
     if (seqlen!=len) {
        if (V)
           print "seq2seq: HOW sequence length error in \"" name "\"" \
	   	>"/dev/stderr";
        cec[eformat]++;
        return 0;
     }
     else if (seqass!=len) {
        if (V)
           print "seq2seq: HOW assignment length error in \"" name "\"" \
	   	>"/dev/stderr";
        cec[eformat]++;
        return 0;
     }
  }

  # print ...
  oec[OUTPUTFORMAT]++;
  if (OUTPUTFORMAT=="how") {
     if (eformat=="fasta") {
        for (k=1; k<=len; k++)
	    ass[k]=".";     
     }
     printhow(OUTDIR,len,name,seq,ass);
  }
  else
     printfasta(OUTDIR,len,name,seq);

  return 0;
}

###############################################################################
function reset() {
###############################################################################

  name		= "";
  seqlen	= 0;		# actual sequence length
  asslen	= 0;		# actual assignment length

  corruption	= 0;		# new entry
  nonstd	= 0;		# flag: non-standard symbols detected

  isfastaseq	= 0;		# flag: picking up FASTA sequence
  ishowass	= 0;		# flag: picking up HOW assignment
  ishowseq	= 0;		# flag: picking up HOW sequence


  return 0;

}

###############################################################################
function report() {
###############################################################################

  print "seq2seq 1.0 report" >"/dev/stderr";
  print "----------------------------------------------" >"/dev/stderr";
  print "# entries\t\t   how\t fasta\t total" >"/dev/stderr";
  print "----------------------------------------------" >"/dev/stderr";
  printf("INPUT\t\t\t%6d\t%6d",ec["how"],ec["fasta"]) >"/dev/stderr";
  printf("\t%6d\n",ec["how"]+ec["fasta"]) >"/dev/stderr";
  printf("Corrupted\t\t%6d\t%6d",cec["how"],cec["fasta"]) >"/dev/stderr";
  printf("\t%6d\n",cec["how"]+cec["fasta"]) >"/dev/stderr";
  printf("Zero-length\t\t%6d\t%6d",zec["how"],zec["fasta"]) >"/dev/stderr";
  printf("\t%6d\n",zec["how"]+zec["fasta"]) >"/dev/stderr";
  printf("Too long\t\t%6d\t%6d",Tec["how"],Tec["fasta"]) >"/dev/stderr";
  printf("\t%6d\n",Tec["how"]+Tec["fasta"]) >"/dev/stderr";
  printf("Too short\t\t%6d\t%6d",Bec["how"],Bec["fasta"]) >"/dev/stderr";
  printf("\t%6d\n",Bec["how"]+Bec["fasta"]) >"/dev/stderr";
  if ((N!="") || (P!="")) {
     printf("With non-std symbols\t%6d\t%6d",xec["how"],xec["fasta"]) \
     	>"/dev/stderr";
     printf("\t%6d\n",xec["how"]+xec["fasta"]) >"/dev/stderr";
  }
  printf("With modified name\t%6d\t%6d",mec["how"],mec["fasta"]) >"/dev/stderr";
  printf("\t%6d\n",mec["how"]+mec["fasta"]) >"/dev/stderr";
  if (MAXNAMELEN) {
     printf("With truncated name\t%6d\t%6d",tec["how"],tec["fasta"]) \
     	>"/dev/stderr";
     printf("\t%6d\n",tec["how"]+tec["fasta"]) >"/dev/stderr";
  }
  printf("With non-unique name\t%6d\t%6d",dec["how"],dec["fasta"]) >"/dev/stderr";
  printf("\t%6d\n",dec["how"]+dec["fasta"]) >"/dev/stderr";
  printf("OUTPUT\t\t\t%6d\t%6d",oec["how"],oec["fasta"]) >"/dev/stderr";
  printf("\t%6d\n",oec["how"]+oec["fasta"]) >"/dev/stderr";
  print "----------------------------------------------" >"/dev/stderr";

  return;
}

###############################################################################
function printhow(OUTDIR,len,name,seq,ass) {
###############################################################################

  # define output file name ...................................................
  if (OUTDIR=="/dev/stdout")
     outfile = OUTDIR;
  else
     outfile = (OUTDIR "/" name);

  # PRINT .....................................................................
  printf("%6d %s\n",len,name) >outfile;		# header

  for (j=1; j<=len; j++) {			# sequence
      printf(seq[j]) >outfile;
      if (!(j%80))
         printf(" %7d\n",j) >outfile;
  }
  if (--j%80)
     print "" >outfile;

  for (j=1; j<=len; j++) {			# assignment
      printf(ass[j]) >outfile;
      if (!(j%80))
         printf(" %7d\n",j) >outfile;
  }
  if (--j%80)
     print "" >outfile;

  # close output file .........................................................
  if (OUTDIR!="/dev/stdout")
     close(outfile);

  return 0;
}

###############################################################################
function printfasta(OUTDIR,len,name,seq) {
###############################################################################

  # define output file name ...................................................
  if (OUTDIR=="/dev/stdout")
     outfile = OUTDIR;
  else
     outfile = (OUTDIR "/" name);

  # PRINT .....................................................................
  print ">" name >outfile;			# header

  for (j=1; j<=len; j++) {			# sequence
      printf(seq[j]) >outfile;
      if (!(j%60))
         print "" >outfile;
  }
  if (j%60)
     print "" >outfile;

  # close output file .........................................................
  if (OUTDIR!="/dev/stdout")
     close(outfile);

  return 0;

}
