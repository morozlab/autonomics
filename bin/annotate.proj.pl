#!/usr/bin/perl -w
use strict;
use Sys::Hostname;
use File::Basename;
use Cwd;
use Getopt::Long;

my($proj, $mira, $data, $cpe, $paired_end, $noass, $bo, $qo, $assonly);

my $cmd_line = "$0 @ARGV";
my $args0 = "@ARGV";

GetOptions(
  "proj=s" => \$proj,
  "mira" => \$mira,
  "data=s" => \$data,
  "noass" => \$noass,
  "convert_paired_end" => \$cpe,
  "bo" => \$bo,
  "qo" => \$qo,
  "assonly" => \$assonly,
  "paired_end" => \$paired_end,
);

sub printUsage{
  print "$0\n";
  print "usage:\n";
  print "	-proj :project_name (REQUIRED)\n";
  print "       -mira   (OPTIONAL)\n";
  print "       -paired_end   (OPTIONAL)\n";
  print '       -convert_paired_end   (OPTIONAL. Use to convert fastq file that have \1 , \2 at end of "@" line to just 1 , 2)\n';
  print "       -assonly   (OPTIONAL)\n";
  print "       -noass   (OPTIONAL)\n";
  print "       -bo   (run blast_nr only, requires -data to be set also  OPTIONAL)\n";
  print "       -qo   (run quantification only  OPTIONAL)\n";
  print "       -data   (OPTIONAL [ NT | AA ] - REQUIRED if -noass used)\n";
};

# python /home/pwilliams/autonomics/scripts/systemtools.py --add-project --assign-workflow --add-jobs adapter_trim quality_trim read_normalization assemble quantification blast_nr blast_swissprot pfam kegg go --set-config adapter_trim:+ quality_trim:+ read_normalization:+ assemble:+ quantification:+ blast_nr:+ blast_swissprot:+ pfam:+ kegg:+ go:+ --set-args "quality_trim|pipeline_args;--paired-end" "adapter_trim|pipeline_args;--paired-end" "read_normalization|pipeline_args;--paired-end" "assemble|pipeline_args;--paired-end" --project-names Aplysia_californica_hatchling_veliger_12day_HiSeq_trans

# plus see end of this file

my $assemble;

if(defined($bo)) { $noass = 1; $bo = 1 ;  } else { $bo = 0; }
if(defined($qo)) { $qo = 1 ;  } else { $qo = 0; }
if(defined($cpe)) { $cpe = 1; } else { $cpe = 0; }

if(!defined($proj)) { printUsage(); exit 0; }
if(!defined($mira)) { $mira = 0; }  else { $mira = 1; }
if(!defined($noass)) { $noass = 0; $assemble = 1; }  else { $noass = 1; $assemble = 0;}

if(!defined($paired_end)) { $paired_end = 0; } else { $paired_end = 1; }

my $projdir = $ENV{ PROJECT_PATH };
$projdir .= "/";
$projdir .= $proj;
$projdir .= "/";

if ($noass && (!defined($data))) {
  print "-data NT | AA must be specified when -noass or -bo is set\n";
  exit 0;
}

if ($noass && (($data ne "NT") && ($data ne "AA"))){
  print "-data must be NT or AA\n";
  exit 0;
}

my $blastp = 0;

if ($noass) {
  if ($data eq "AA") { $blastp = 1; }
}

if ($paired_end && $mira) {
  print " can not set paired_end and mira flags at same time\n";
  exit 0;
}

if ($cpe && ((not $paired_end) || $bo || $mira || $noass)) { 
  print " can not set convert_paired_end with these args\n";
  exit 0;
}

if ($mira && $noass) {
  print " can not set mira flag with noass flag\n";
  exit 0;
}

if (!-d $projdir) {
  print "$projdir does not exist\n";
  exit 0;
}

my $fastq;
my $fq =  $projdir . $proj . ".fastq";
if (-e $fq) {
    print "this project has a fastq file\n";
    $fastq = 1;
} else {
    print "this project does not have a fastq file ... skipping quantification!\n";
    $fastq = 0;
}

my $fq2 = $projdir . $proj . ".fastq.end2";
if ((! $noass) && (-e $fq2) && (!$paired_end)) {
  print "this project has paired end data but the paired_end flag is not set!\n";
  exit 0;
}

my $fa;
if ($noass) {
  $fa = $projdir . $proj . "_project.fasta";
} else {
  $fa = $projdir . $proj . ".fastq";
}

if (not -e $fa) { 
  print "$fa is REQUIRED but does not exist!\n";
  exit 0;
}

if ($qo) {
  if (not $fastq) { 
    print "$fq is REQUIRED but does not exist!\n";
    exit 0;
  }
}

if ($paired_end) {
  if (not -e $fq2) {
    print "$fq2 is REQUIRED but does not exist!\n";
    exit 0;
  }
}

my $quant;
if ($paired_end) {
    $quant = "\"quantification|pipeline_args;--paired_end --aligner bowtie --db-type NT |resources;cpu:400\"";
} else {
    $quant = "";
}
# --aligner bowtie --db-type NT --paired-end | cpu:400 
#job_type|arg_name1;arg_val1|arg_name2;arg_val2

my $cmd;

my $proj_proj =  $projdir . $proj;
if ($cpe) {
  $cmd = "convert.paired.end.format $proj_proj";
  print "\n\n $cmd\n\n\n WAIT A MINUTE! \n\n";
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

my $arg = "";
if ($mira) {
  $arg = "--set-args \"assemble|pipeline_args;--assembler mira\"";
}

# mira, paired_end, noass, blastp

if ($paired_end && $assemble) {
  $arg = "--set-args \"quality_trim|pipeline_args;--paired-end\" \"adapter_trim|pipeline_args;--paired-end\" \"read_normalization|pipeline_args;--paired-end\" \"assemble|pipeline_args;--paired-end\" $quant";
}

if ($paired_end && $noass) {
  $arg = "--set-args $quant";
}

if ($blastp) {
  $arg = "--set-args \"blast_nr|pipeline_args;--aligner blastp\" \"blast_swissprot|pipeline_args;--aligner blastp\" \"pfam|pipeline_args;\"  $quant";
}



my $scripts_path = $ENV{AUTONOMICS_PATH};
$scripts_path .= '/scripts';

my $systools = $scripts_path . "/systemtools.py";

if ($noass) {
    if ($bo) { #blast_nr only
      $cmd = "python $systools --add-project --assign-workflow" .
             " --add-jobs blast_nr --set-config blast_nr:+ $arg --project-names $proj";
    } elsif ($qo) { #quantification only
      $cmd = "python $systools --add-project --assign-workflow" .
             " --add-jobs quantification --set-config quantification:+ $arg --project-names $proj";
    } else { # run everything except assembly
      if (($data eq "AA") || (not $fastq)) {  # omit quantification since no reads either for genemodels or no reads
        $cmd = "python $systools --add-project --assign-workflow" .
               " --add-jobs blast_nr blast_swissprot pfam kegg go" . 
               " --set-config blast_nr:+ blast_swissprot:+ pfam:+ kegg:+ go:+" .
               " $arg --project-names $proj";
      } else {  # no assy but run quantification too
        $cmd = "python $systools --add-project --assign-workflow" .
               " --add-jobs blast_nr blast_swissprot pfam kegg go quantification" .
               " --set-config blast_nr:+ blast_swissprot:+ pfam:+ kegg:+ go:+ quantification:+ $arg --project-names $proj";
      }
    }
} elsif ($assonly) {
  $cmd = "python $systools --add-project --assign-workflow" .
          " --add-jobs adapter_trim quality_trim read_normalization assemble" .
          " --set-config adapter_trim:+ quality_trim:+ read_normalization:+ assemble:+" .
          " $arg --project-names $proj";
} else { # assemble and everything else
  $cmd = "python $systools --add-project --assign-workflow" .
          " --add-jobs adapter_trim quality_trim read_normalization assemble quantification blast_nr blast_swissprot pfam kegg go" .
          " --set-config adapter_trim:+ quality_trim:+ read_normalization:+ assemble:+ quantification:+ blast_nr:+ blast_swissprot:+ pfam:+ kegg:+ go:+" .
          " $arg --project-names $proj";
}

print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

print "job submitted successfully\n";
