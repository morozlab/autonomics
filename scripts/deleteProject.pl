#!/apps/perl/perls/perl-5.16.0/bin/perl -w

#  #!/usr/bin/perl -w

use warnings;
use strict;
use DBI;
use File::Basename;
use File::Path;
# use Cwd;
use POSIX;

my %nb_abrevs = (
'misc' => 'misc',
'apt' => 'aplysia',
'apn' => 'aplysia2',
'pb' => 'pleurobrachia',
'sandbox' => 'sandbox',
'genomes' => 'genomes',
'sponges' => 'porifera',
'nb1' => 'nb1',
'molluscs' => 'molluscs',
#'snails' => 'gastropoda',
'secr' => 'secretoryMolecules',
'verts' => 'vertebrates',
'squid' => 'cephalopods',
'combs' => 'ctenophora',
    );


my ($pname, $nb_abrev) = @ARGV;
if ((not defined $pname) || (not defined $nb_abrev)) {
   print "\nUsage: $0 <project_name> <nb_abrev>\n";
   print "NB_abrev  NB_name\n";
   while( my( $key, $value ) = each %nb_abrevs ){
	print "$key\t $value\n";
   }
    exit 0;
}

my $nb_name;
if( exists($nb_abrevs{$nb_abrev} ) ){
    $nb_name = "nb_".$nb_abrevs{$nb_abrev};
} else {
    print "\nUsage: $0 <nb_abrev>  <project_id>\n";
    print "\nINVALID NB_ABREV: $nb_abrev    here are the valid ones\n\n";
    print "NB_abrev  NB_name\n";
    while( my( $key, $value ) = each %nb_abrevs ){
	print "$key\t $value\n";
    }
    print "\n";
    exit 0;
}

my $log = "log.del." . $nb_abrev . "." . $pname . "." . $$;
my $cmd = "delete_project.pl $nb_abrev $pname |& tee $log";
# print "$cmd\n";
system ($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }
