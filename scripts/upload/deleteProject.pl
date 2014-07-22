#!/apps/perl/perls/perl-5.16.0/bin/perl -w
use warnings;
use strict;
use DBI;
use File::Basename;
use File::Path;
use POSIX;
use nb_names;

my ($pname, $nb_abrev) = @ARGV;
if ((not defined $pname) || (not defined $nb_abrev)) {
   print "\nUsage: $0 <project_name> <nb_abrev>\n";
   print "NB_abrev  NB_name\n";
   nb_names::print_valid_nb_abrevs();
   exit 0;
}

my $nb_name = nb_names::nb_abrev_lookup($nb_abrev);
if ($nb_name) {
    print "nb_name: $nb_name\n";
} else {
    print "\nUsage: $0 <nb_abrev>  <project_name>\n";
    print "\nINVALID NB_ABREV: $nb_abrev    here are the valid ones\n\n";
    print "NB_abrev  NB_name\n";
    nb_names::print_valid_nb_abrevs();
    print "\n";
    exit 0;
}

my $log = "log.del." . $nb_abrev . "." . $pname . "." . $$;
my $cmd = "delete_project.pl $nb_name $pname |& tee $log";
# print "$cmd\n";
system ($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }
