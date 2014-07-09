#!/usr/bin/perl -w
use strict;

my ($old_name, $new_name) = @ARGV;
if ((not defined $old_name) ||
    (not defined $new_name)) {
    print "\nUsage: $0 <existing_proj_name> <new_proj_name>\n";
    print " You need to be in dir where old proj name exists\n";
    exit 0;
}

$old_name =~ s{/\z}{};
$new_name =~ s{/\z}{};

print "copying $old_name to $new_name\n";

my $cmd = "mkdir $new_name";
print "$cmd\n";
 system($cmd);
 if ( $? ) { die "Command failed: $cmd: $!"; }

$cmd = "cp -r $old_name/* $new_name";
print "$cmd\n";
  system($cmd);
 if ( $? ) { die "Command failed: $cmd: $!"; }

opendir(DIR,  $new_name) or die "Unable to opendir $new_name:";
while (my $file = readdir(DIR)) {
    next if ($file =~ m/^\./);
    if ($file =~ m/$old_name.*/) {
      my $old_file = $file;
      print "file before: $old_file\n";
      $file =~ s/$old_name/$new_name/;
      print "file after: $file\n";
      $cmd = "mv $new_name/$old_file $new_name/$file";
      print "$cmd\n";
      system($cmd);
      if ( $? ) { die "Command failed: $cmd: $!"; }
    }
}
close DIR;

chdir $new_name or die "cant change dirs";

my $gocats_name = $new_name . "_gocats.txt";
my $bu_name = $gocats_name .".orig";
$cmd = "cp $gocats_name $bu_name";
print "$cmd\n";
 system($cmd);
 if ( $? ) { die "Command failed: $cmd: $!"; }

$cmd = "sed -i \'s/$old_name/$new_name/g\' $gocats_name";
print "$cmd\n";
 system($cmd);
 if ( $? ) { die "Command failed: $cmd: $!"; }




