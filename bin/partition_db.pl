#!/usr/bin/perl -w
use strict;
use File::Basename;

#####################################################################
#
#                    partition_db.pl
#                    ==========
#
#           Peter Williams    LLNL    May 2007
#
# This code assumes that the database has been processed by sort_db.pl
#
# Call: partition_db.pl <number_of_partitions> <directory_with_sorted_database>
#
# So if your call to sort_db.pl looked like:
#     sort_db.pl <path>/NT <output_dir_with_path>
# Then <directory_with_sorted_database> = <output_dir_with_path>/NT
#
# For example for 4 partitions, the directory: 
# /usr/local/partitioned_databases/NT/NT_4_parts is created and the 
# formatted partitions NT_frag0, NT_frag0.nhr, .. , NT_frag3, NT_frag3.nsq
# are placed therein.
#
# this code assumes that the path to the blast command 'formatdb' is
# in your $PATH environment variable.  If not, edit this file and add 
# the path to formatdb where it appears towards the end of this file.
#
#####################################################################


my ($nfrags,$sorted_db, $part_file, $ext) = @ARGV;
if ((not defined $sorted_db) ||
    (not defined $part_file) ||  # <proj>_project_blast_nr need to add _NN.fasta
    (not defined $ext) ||
    (not defined $nfrags)) {
    print "\nUsage: $0 <nfrags> <path_sorted_database> <part_file_name> <ext>\n";
    exit 0;
}

my $proj_dir = dirname($sorted_db);
my $proj_name = basename($proj_dir);

my @out;
for (my $i=1 ; $i <= $nfrags ; $i++) {
    my $fh;
    my $fn = $proj_dir . '/' . $part_file . "_" . $i . $ext;
    open( $fh, ">$fn") || die ("could not open $fn\n");
    push(@out, $fh);
}

my $fragment = -1;
my $dir = 1;

open(FILE,  $sorted_db) || die ("could not open input  $sorted_db\n");
while( my $line = <FILE> ) {
    if ($line =~ /^>/) {
	$fragment += $dir;
	if ($fragment == $nfrags) {
	    $fragment = $nfrags-1;
	    $dir = -1;
	} elsif ($fragment == -1) {
	    $fragment = 0;
	    $dir = 1;
	}
    }

    my $fh = $out[$fragment];
    print $fh ($line);
}
close(FILE) || die ("could not close input fa\n");

for (my $i=0 ; $i < $nfrags ; $i++) {
    close($out[$i]) || die ("could not close fragment $i\n");
}

=stop
for (my $i=0; $i < $nfrags; $i++) {
  my $fn = $part_dir . '/' . $db_name ."_frag" . $i;
  my $cmd = "formatdb -l $part_dir/formatdb.log -i $fn -p F";
  print "$cmd\n";
  `$cmd`;
}
=cut
