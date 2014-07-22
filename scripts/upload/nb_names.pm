package nb_names;
use strict;
use warnings;

#require Exporter;
#@ISA = qw(Exporter);
#@EXPORT = qw( nb_abrev_exists print_valid_nb_abrevs);
 
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

sub nb_abrev_lookup {
   my ($nb_abrev) = @_;
   if( exists($nb_abrevs{$nb_abrev})) { return "nb_" . $nb_abrevs{$nb_abrev}; }
   else { return 0; }
}

sub print_valid_nb_abrevs {
    while( my( $key, $value ) = each %nb_abrevs ){
	print "$key\t $value\n";
    }
}

1;

