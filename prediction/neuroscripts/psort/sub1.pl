
############################################################
##  PSORT II:  subroutine package 1 (980607)
##                      copyright: Kenta Nakai 1996-1998
############################################################
## 
## 971205: bug fix in 'mtop'
## 980508: small change in 'gvh'
## 980508: neural network for cyt/nuc discrimination by A.Reinhardt

## Sequence Numbering convention:
##   argument: 1, 2, 3,...   internal: 0, 1, 2,...


package sub1;

## definitions ###################################################

%hydro = (
     'A'=> 1.8, 'R'=>-4.5, 'N'=>-3.5, 'D'=>-3.5, 'C'=> 2.5,
     'Q'=>-3.5, 'E'=>-3.5, 'G'=>-0.4, 'H'=>-3.2, 'I'=> 4.5,
     'L'=> 3.8, 'K'=>-3.9, 'M'=> 1.9, 'F'=> 2.8, 'P'=>-1.6,
     'S'=>-0.8, 'T'=>-0.7, 'W'=>-0.9, 'Y'=>-1.3, 'V'=> 4.2,
     'B'=>-3.5, 'Z'=>-3.5, 'X'=>-0.5
);

### The GES scale (Engelman et al, 1986) X -1 ###
%hydro2 = (
     'A'=>  6.7, 'R'=>-51.5, 'N'=>-20.1, 'D'=>-38.5, 'C'=>  8.4,
     'Q'=>-17.2, 'E'=>-34.3, 'G'=>  4.2, 'H'=>-12.6, 'I'=> 13.0,
     'L'=> 11.7, 'K'=>-36.8, 'M'=> 14.2, 'F'=> 15.5, 'P'=> -0.8,
     'S'=>  2.5, 'T'=>  5.0, 'W'=>  7.9, 'Y'=> -2.9, 'V'=> 10.9,
     'B'=>  0.0, 'Z'=>  0.0, 'X'=>  0.0
);

%aanum = (
     'A'=> 0, 'R'=> 1, 'N'=> 2, 'D'=> 3, 'C'=> 4,
     'Q'=> 5, 'E'=> 6, 'G'=> 7, 'H'=> 8, 'I'=> 9,
     'L'=>10, 'K'=>11, 'M'=>12, 'F'=>13, 'P'=>14,
     'S'=>15, 'T'=>16, 'W'=>17, 'Y'=>18, 'V'=>19,
     'B'=>20, 'Z'=>21, 'X'=>22
);

## new McGeoch's method ###################################################
sub psg {
    $CRMAX = 11;  ## maximum position of CR_end (external numbering )
    my($SEGLEN) = 8;
    my($PSGTHR) = 4.4;

    my($cr_leng, $poschg, $negchg, $ur_leng) = &regions;

    $ur_leng2 = ($cr_leng+$ur_leng > 30) ? 30 - $cr_leng : $ur_leng;
    return -10000 if $cr_leng+$SEGLEN > $::leng;

    my($ur_peak) = &hydseg2($SEGLEN, $cr_leng+1, $ur_leng2);

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#esig\">PSG:  a new signal peptide prediction method</A>\n";
    }elsif($::OPT_V) {
	print "PSG:  a new signal peptide prediction method\n";
    }
## ad hoc rules ##
    if($ur_leng >= 60 || ($poschg-$negchg) < 0) {
	$ur_peak = 0;
	print ">>> too long uncharged segment: $ur_leng\n"
			if($ur_leng >= 60 && $OPT_V);
	printf ">>> negative net charge in the N-region: %3.0f\n", 
		$poschg-$negchg if($poschg-$negchg < 0 && $OPT_V);
    }
    if($::OPT_V) {
	printf("      N-region:  length %d;  pos.chg %d;  neg.chg %d\n",
						$cr_leng, $poschg, $negchg);
	printf("      H-region:  length %d;  peak value %6.2f\n",
              					$ur_leng, $ur_peak);
	printf("      PSG score: %6.2f\n\n", $ur_peak - $PSGTHR);
    }
	
    $ur_peak - $PSGTHR;
}

# specifically for gram-positive bacteria
sub psggp {
    $CRMAX = 11;  ## maximum position of CR_end (external numbering )
    my($SEGLEN) = 8;
    my($PSGTHR) = 4.4;

    my($cr_leng, $poschg, $negchg, $ur_leng) = &regions;

    $ur_leng2 = ($cr_leng+$ur_leng > 30) ? 30 - $cr_leng : $ur_leng;
    return -10000 if $cr_leng+$SEGLEN > $::leng;

    my($ur_peak) = &hydseg2($SEGLEN, $cr_leng+1, $ur_leng2);

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#esig\">PSGGP:  a new signal peptide prediction method for Gram-positive bacteria</A>\n";
    } elsif($::OPT_V) {
	print "PSGGP:  a new signal peptide prediction method for Gram-positive bacteria\n";
    }
## ad hoc rules ##
    if($ur_leng >= 60 || ($poschg-$negchg) < 0) {
	$ur_peak = 0;
	print ">>> too long uncharged segment: $ur_leng\n"
		if($ur_leng >= 60 && $OPT_V);
	printf ">>> negative net charge in the N-region: %3.0f\n", 
			$poschg-$negchg if($poschg-$negchg < 0 && $OPT_V);
    }
    if($::OPT_V) {
	printf("      N-region:  length %d;  pos.chg %d;  neg.chg %d\n",
						$cr_leng, $poschg, $negchg);
	printf("      H-region:  length %d;  peak value %6.2f\n",
              					$ur_leng, $ur_peak);
	printf("      PSG score: %6.2f\n\n", $ur_peak - $PSGTHR);
    }
	
    $ur_peak - $PSGTHR;
}

# return the postion of most C-term. charged residue
# (output 0 means 'no charge'),
# the charges within a region, and the ur_leng;

sub regions {
# calculate the CR length
    my($i, $aa);
    for ($i=$CRMAX-1; $i>=0; $i--) {
                $aa = $::seq[$i];
                last if($aa eq 'R' || $aa eq 'D' || $aa eq 'E' || $aa eq 'K');        }
    my($crend) = $i;

# count charged residues
    my($poschg, $negchg) = (0, 0);
    for ($i=0; $i<=$crend; $i++) {
                $aa = $::seq[$i];
                if($aa eq 'R' || $aa eq 'K') {$poschg++}
                elsif($aa eq 'D' || $aa eq 'E') {$negchg++}
    }

# calculate the UR length
    my($urleng) = 0;
    for($i=$crend+1; $i<$::leng; $i++, $urleng++) {
                $aa = $::seq[$i];
                last if($aa eq 'R' || $aa eq 'D' || $aa eq 'E' || $aa eq 'K');        }
        #print substr($::seq, $start, $urleng), "  ";

        ($crend+1, $poschg, $negchg, $urleng);
}

### most hydrophobic segment
### Args: SEGLENG, START, RANGE0
sub hydseg2 {
        my($segleng, $start, $urleng) = @_;
        die "*** too short sequence\n" if $start + $segleng - 1 > $::leng;
        $start--;  ### internal notation
        my($range) = ($urleng - $segleng < 0) ? 0 : ($urleng - $segleng);

        my(@aas) = split(//, substr($::seq, $start, $segleng+$range));
        my($i);
        my($sum) = 0;
        foreach $i (1..$segleng) {
#			die ">>>$id: $aas[$i-1]} undefined" 
#				unless defined($hydro2{$aas[$i-1]});
			$sum += $hydro2{$aas[$i-1]};
        }
        my($hmax, $imax) = ($sum, 0);
        foreach $i (1..$range) {
#                die ">>>$id: $aas[$i+$segleng-1]} undefined" 
#                        unless defined($hydro2{$aas[$i+$segleng-1]});
                $sum = $sum - $hydro2{$aas[$i-1]} + $hydro2{$aas[$i+$segleng-1]};
                if($sum > $hmax) {
                        $hmax = $sum;
                        $imax = $i;
                }
        }
        $hmax /= $segleng;
}


## von Heijne's method for signal seq prediction ##################

sub gvh {

    package GvH;

    $flg = $_[0];
    my ($i, $j, $end);

# parameter for prokaryotes
    @prm_prk = (
    [ 1.14, 0.92, 0.92, 1.03, 0.63, 0.78, 0.45, 0.63, 0.78, 0.78,
      2.01,-0.47, 2.27, 1.73, 0.22 ],
    [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     -3.58, 0.00,-3.58, 0.00, 0.00 ],
    [-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,
     -3.58,-0.69,-3.58, 0.00, 1.39 ],
    [-0.79,-0.79,-0.79,-0.79,-0.79,-0.79,-0.79,-0.79,-0.79,-0.79,
     -3.58,-0.79,-3.58, 0.60, 1.29 ],
    [ 0.43, 1.12, 0.84, 1.12,-0.26,-0.26, 1.82,-0.26, 1.12,-0.26,
     -3.58, 1.68,-3.58,-0.26,-0.26 ],
    [ 0.39,-0.30,-0.30,-0.30, 0.11, 0.62,-0.30, 0.39,-0.30,-0.30,
     -3.58,-0.30,-0.30,-0.99,-0.99 ],
    [ 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22,
     -3.58, 2.17,-3.58, 0.22, 0.22 ],
    [ 0.57,-0.53, 1.08,-0.53, 1.08,-0.53,-0.53, 0.57,-0.53,-0.53,
     -3.58,-0.53,-3.58,-0.53, 0.16 ],
    [-0.92,-0.92,-0.92,-0.92,-0.92,-0.92,-0.92,-0.92,-0.92,-0.92,
     -3.58,-0.22,-3.58, 0.18,-0.92 ],
    [ 1.09, 1.40, 1.20, 1.09, 1.20, 1.57,-0.99,-0.99,-0.30,-0.30,
     -0.99,-0.30,-3.58,-0.99,-0.99 ],
    [ 0.51, 1.20, 0.51, 0.51, 1.61, 1.20, 1.61, 0.51, 0.51, 1.20,
     -3.58, 1.90,-3.58, 0.51, 0.51 ],
    [-0.47,-0.47,-0.47,-0.47,-0.47,-0.47,-0.47,-0.47,-0.47,-0.47,
     -3.58, 0.63,-3.58,-0.47, 0.92 ],
    [-0.53,-0.53,-0.53,-0.53,-0.53,-0.53, 0.16, 0.57, 1.08, 0.16,
     -3.58,-0.53,-3.58,-0.53, 1.08 ],
    [-0.34,-0.34,-0.34,-0.34,-0.34,-0.34,-0.34,-0.34, 0.36, 0.36,
     -3.58, 0.76,-3.58,-0.34,-0.34 ],
    [-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,
     -3.58,-0.53,-3.58,-0.53,-0.53 ],
    [-0.96,-0.96,-0.96, 0.43, 0.43,-0.96, 0.65, 1.75, 0.65, 1.12,
      0.65,-0.26,-0.26,-0.96,-0.96 ],
    [-0.10,-0.79, 0.60,-0.10,-0.10,-0.10,-0.10,-0.10, 0.82,-0.79,
      0.31,-0.79,-0.79,-0.79,-0.10 ],
    [ 0.69, 1.03,-0.92, 0.18,-0.92, 0.47, 1.03,-0.92,-0.92, 0.47,
      0.18,-0.92,-3.58,-0.22,-0.92 ],
    [ 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92,
     -3.58, 0.92,-3.58, 0.92, 0.92 ],
    [-0.26,-0.26,-0.26,-0.26,-0.26,-0.26,-0.26,-0.26,-0.26, 0.84,
     -3.58,-0.26,-3.58,-0.26,-0.26 ]
   );

# parameter for eukaryotes
    @prm_euk = (
    [ 0.10,-0.11,-0.04, 0.03, 0.32, 0.22, 0.22, 0.16, 0.54, 0.03,
      1.18,-0.88, 1.71, 0.22,-0.88 ],
    [-0.41, 0.29, 0.69, 0.44, 0.69, 1.13, 0.29, 0.58, 0.11, 0.29,
      1.44,-0.41, 0.69, 0.58,-0.41 ],
    [-2.19,-2.19,-2.19,-2.19,-2.19,-2.19,-2.19,-2.19,-0.58,-1.09,
     -5.08,-0.58,-5.08, 0.12, 0.21 ],
    [-2.30,-2.30,-2.30,-2.30,-2.30,-2.30,-2.30,-2.30,-1.20,-0.36,
     -5.08,-0.36,-5.08, 0.26, 0.34 ],
    [ 0.84, 0.47, 0.68, 0.68, 0.07, 0.22, 1.17, 0.84,-0.34,-0.11,
     -5.08, 0.84,-5.08, 0.07,-0.34 ],
    [-1.11,-1.11,-1.39,-0.70,-1.39, 0.07,-1.39,-1.80, 0.45, 1.03,
     -0.88,-0.55, 1.17,-0.19,-0.55 ],
    [-1.22,-1.22,-1.22,-1.22,-1.22,-1.22,-1.22,-1.22, 0.39,-1.22,
     -5.08, 0.57,-5.08, 0.16,-0.53 ],
    [ 0.71, 0.71, 0.08,-0.21, 0.40,-0.39,-0.62, 0.08,-0.39,-2.00,
      0.30,-0.39,-5.08, 0.08,-0.06 ],
    [-2.42,-2.42,-2.42,-2.42,-2.42,-2.42,-2.42,-2.42,-2.42,-1.04,
     -5.08,-1.73,-5.08,-0.03,-0.23 ],
    [ 1.77, 1.73, 1.78, 1.88, 1.86, 1.31, 1.67, 1.40,-0.19, 0.64,
     -0.41, 0.50,-2.49,-0.41,-1.11 ],
    [-0.99, 0.11, 0.95, 0.39,-0.99, 0.80,-0.30,-0.30,-0.99,-0.99,
     -5.08,-0.99,-5.08,-0.99,-0.30 ],
    [-1.96,-1.96,-1.96,-1.96,-1.96,-1.96,-1.96,-1.96,-0.86,-0.86,
     -5.08, 0.34,-5.08,-0.57,-0.01 ],
    [-1.31,-2.00,-1.31,-2.00,-2.00,-0.62,-2.00, 0.08, 0.99, 0.64,
     -5.08,-2.00,-0.90,-2.00, 1.09 ],
    [-1.84,-1.84,-1.84,-1.84,-1.84,-0.05,-1.84,-1.84, 0.46, 0.24,
     -5.08, 1.05,-0.74, 1.10, 0.46 ],
    [-1.34,-2.03,-2.03,-2.03,-2.03,-2.03,-2.03,-2.03,-0.08,-0.64,
     -5.08, 0.68,-5.08, 0.46, 0.17 ],
    [-0.24,-1.34,-0.35,-0.64, 0.13,-0.13, 0.27, 0.34, 0.82,-0.04,
      0.70, 0.40, 0.56, 0.27,-0.13 ],
    [-1.58, 0.03,-0.66,-0.89,-0.66, 0.29,-0.33,-0.33, 0.21,-0.48,
      0.56,-0.19,-0.48,-1.17, 0.03 ],
    [ 0.59, 0.81, 0.30, 0.48, 0.16, 0.30,-0.01, 0.89,-2.41, 0.08,
      1.06,-1.31,-5.08,-0.33, 0.43 ],
    [ 0.80, 0.51, 0.51,-0.59,-0.59, 0.11, 1.20, 0.51,-0.59, 0.51,
     -5.08, 1.61,-5.08, 0.11,-0.59 ],
    [-1.72,-1.72,-0.34,-1.72,-1.72,-1.72,-0.62,-1.72,-1.72,-1.03,
     -5.08,-0.11,-5.08,-1.72, 0.22 ]
   );

## GvH numbering method
    %aagvh = (
     'A'=> 0, 'C'=> 1, 'D'=> 2, 'E'=> 3, 'F'=> 4,
     'G'=> 5, 'H'=> 6, 'I'=> 7, 'K'=> 8, 'L'=> 9,
     'M'=>10, 'N'=>11, 'P'=>12, 'Q'=>13, 'R'=>14,
     'S'=>15, 'T'=>16, 'V'=>17, 'W'=>18, 'Y'=>19
    );

    my($imax, $smax) = (-1, -1000);
    	## position of max score  ## max score
    $end = ($::leng - 15 > 48) ? (48) : ($::leng - 15);
#    my($end) = ($::leng - 15 > 35) ? (35) : ($::leng - 15);
    my($score);

    foreach $i (0..$end) {
	$score = 0;
	foreach $j (0..14) {
	    next unless defined($aagvh{$::seq[$i+$j]});
	    if ($flg eq "prokaryote") {
		$score += $prm_prk[$aagvh{$::seq[$i+$j]}][$j];
	    } elsif($flg eq "eukaryote") {
		$score += $prm_euk[$aagvh{$::seq[$i+$j]}][$j];
	    } else {die "gvh: invalid argument"}
	}
#	printf "%2d: %5.2f\n", $i+13, $score-7.5;
	if ($score > $smax) {
	    $smax = $score;
	    $imax = $i;
	}
    }
    $gvhscore = $smax - 7.5;
    $cleavsite = $imax + 13;

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#esig\">GvH:  von Heijne's method for signal seq. recognition</A>\n";
    } elsif($::OPT_V) {
	print "GvH:  von Heijne's method for signal seq. recognition\n";
    }
    if($::OPT_V) {
	printf "      GvH score (threshold: -2.1): %6.2f\n", $gvhscore;
	printf "      possible cleavage site: between %2d and %2d\n\n", 
						$cleavsite, $cleavsite+1;
    }

    ($gvhscore, $cleavsite);
}


## amino acid composition ###########################################
sub mkaac {
    local($start) = $_[0];
        # numbering of $start, $::seq: 0, 1, 2, ..
    local($i, @aac);
    local($len) = $leng - $start;

    foreach $i (0..22) {$aac[$i]=0}

    for($i=$start; $i<=$#seq; $i++) {
	if(defined($aanum{$::seq[$i]})) {
	    $aac[$aanum{$::seq[$i]}]++;
	} else {
	    warn "*** mkaac: illegal char $::seq[$i]: skipped\n";
	}
    }

#    foreach $i (0..9) {printf " %5.1f", $aac[$i]}; print "\n";
#    foreach $i (10..20) {printf " %5.1f", $aac[$i]}; print "\n";
    foreach $i (0..$#aac) {$aac[$i] = $aac[$i] / $len * 100}
    @aac;
}


## Klein et al's ALOM ###############################################

# new version (2 threshold parameters)
sub alom2 {
    my($start, $threshold_loose, $threshold_strict) = @_;
	## threshold_loose (larger value) > threshold_strict (smaller)
    $start -= 1; ## internal residue numbering

    ($A1, $A0, $LSEG) = (-9.02, 14.27, 17);
    return (-1, 999, -1) if ($::leng < 20 || $start > $::leng-$LSEG);

    my($i, $j, $count, @hyd, $xmax, $sum, $imax, $hmax, $x0);
    undef %tms;
	###  %tms is used by other routenes

    foreach $i ($start..$::leng-1) {
	unless(defined $hydro{$::seq[$i]}) {
	    warn "** alom2: illegal char $::seq[$i] at $i+1\n";
	    $hyd[$i] = -100;
	    next;
	}
	$hyd[$i] = $hydro{$::seq[$i]};
    }

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#ealom\">ALOM: Klein et al's method for TM region allocation</A>\n";
    } elsif($::OPT_V) {
	print "ALOM: Klein et al's method for TM region allocation\n";
    }
    printf "      Init position for calculation: %d\n", $start+1 if $::OPT_V;

    for($count=0; ; $count++) {

#### hydseg
	$sum = 0;
	$imax = $start;
	foreach $i ($start..$start+$LSEG-1) {
	    $sum += $hyd[$i];
	}
	$hmax = $sum;

	foreach $i ($start+1..$::leng-$LSEG) {
	    $sum = $sum - $hyd[$i-1] + $hyd[$i+$LSEG-1];
	    if($sum > $hmax) {
		$hmax = $sum;
		$imax = $i;
	    }
	}
	$hmax /= $LSEG;

#### allocation
	$x0 = $A1 * $hmax + $A0;
	$xmax = $x0 if $count==0;  ### maximaum hydrophobicity value

	last if $x0 > $threshold_loose;
	$tms{$imax+1} = $x0;

#### erase the segment
	foreach $i (0..$LSEG-1) {
	    $hyd[$imax+$i] = -10;
	}
    }

#### check the strict threshold value
    printf "      Tentative number of TMS(s) for the threshold %4.1f:  %2d\n", 
			$threshold_loose, $count if $::OPT_V;
    if($count > 0) {
	my($count2) = 0;
	foreach (keys %tms) {
	    $count2++ if $tms{$_} <= $threshold_strict;
	}
### caution!  this may ignore uncleaved signals
	if($count2 < 2) {
	    foreach (keys %tms) {
		delete $tms{$_} if $tms{$_} > $threshold_strict;
	    }
	    $count = $count2;
	    printf "      Number of TMS(s) for threshold %4.1f:  %2d\n", 
		$threshold_loose, $count if $::OPT_V;
	}
    } else {
	printf "      number of TMS(s) .. fixed\n" if $::OPT_V;
    }

#### output
    $mostn = -1;
    unless($count == 0) {
	foreach (sort {$a <=> $b} keys %tms) {
	    if($::OPT_V) {
		printf "      INTEGRAL    Likelihood =%6.2f", $tms{$_};
		printf "   Transmembrane %4d -%4d\n", $_, $_+$LSEG-1;
	    }
	    $mostn = $_ if $mostn < 0;
	}
    }
    if($::OPT_V) {
	printf("      PERIPHERAL  Likelihood =%6.2f (at %d)\n", 
		$x0, $imax+1);
	printf "      ALOM score: %6.2f  (number of TMSs: %d)\n\n", 
			$xmax, $count;
    }

    ($count, $xmax, $mostn);
}

## Membrane Topology ####################################################

sub mtop {
    my($i_middle) = $_[0];
    my($i, $k, $c);
    my($start, $end);
    return -100
		if ($::leng < 20 || $i_middle < 1 || $i_middle > $::leng-1);

  # extention of maximum hydrophobicity segment
    for ($i = $i_middle; $i >= 0; $i--) {
	$c = $::seq[$i];
	last if ( $c eq 'R' || $c eq 'D' || $c eq 'E'
			|| $c eq 'H' || $c eq 'K');
    }
    $start = $i;
    for ($i = $i_middle; $i < $::leng; $i++) {
	$c = $::seq[$i];
	last if ( $c eq 'R' || $c eq 'D' || $c eq 'E'
			|| $c eq 'H' || $c eq 'K');
    }
    $end = $i;

  # charge difference
    my($nchg, $cchg) = (0, 0);
    for($i = $start, $k=0; $i>=0 && $k<15; $i--, $k++) {
	$c = $::seq[$i];
	if($c eq 'D' || $c eq 'E') {
	    $nchg -= 1.0;
	} elsif($c eq 'R' || $c eq 'K') {
	    $nchg += 1.0;
	} elsif($c eq 'H') {
	    $nchg += 0.5;
	}
    }
    $nchg += 1.0 if $i < 0;  ### for unmodified N terminus

    for($i = $end, $k=0; $i<$::leng && $k<15; $i++, $k++) {
	$c = $::seq[$i];
	if($c eq 'D' || $c eq 'E') {
	    $cchg -= 1.0;
	} elsif($c eq 'R' || $c eq 'K') {
	    $cchg += 1.0;
	} elsif($c eq 'H') {
	    $cchg += 0.5;
	}
    }

    $mtopscr = $cchg - $nchg;
    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#mtop\">MTOP: Prediction of membrane topology (Hartmann et al.)</A>\n";
    } elsif($::OPT_V) {
	print "MTOP: Prediction of membrane topology (Hartmann et al.)\n";
    }
    if($::OPT_V) {
	print ' 'x6, "Center position for calculation: $i_middle\n";
	printf "      Charge difference: %4.1f   C(%4.1f) - N(%4.1f)\n", 
						$mtopscr, $cchg, $nchg;
	print ' ' x 6;
	($mtopscr>0) ? print "C > N: C-terminal side will be inside\n\n" :
		   print "N >= C: N-terminal side will be inside\n\n";
    }
    $mtopscr;
}

## assignment of membrane topology

##   output: mtype, start, end
##
###  caution: uses $mostn
sub mtype_assign {
    my($tmsnum, $sig, $mtopscr) = @_;

    my($pos, $c);

    if($tmsnum == 0) {
	('__', 0, 0);
    } elsif($tmsnum == 1) {
	if($sig eq 'cleavable') {
	    $pos = (keys %tms)[0] + $LSEG;
	    $pos=$::leng if $pos > $::leng;  ##### ad hoc
	    ('1a', $pos, $::leng);
	} elsif(($mostn / $::leng) > 0.8) {  ### exceptional topology
	    print ">>> Single TMS is located near the C-terminus\n\n"
				if $::OPT_V;
	    $pos = (keys %tms)[0] - 1;
	    ('Nt', 1, $pos);
	} else {
	    if($mtopscr > 0) {
		('1b', $mostn, $::leng);
	    } else {
		('2 ', 1, $mostn);
	    }
	}
    } else {
	if($mtopscr > 0) {
	    ('3b', 0, 0);
	} else {
	    ('3a', 0, 0);
	}
    }
}

## Contribution from Dr. Astrid Reinhardt ############################
#
#####  CytoplasmicNuclearPredictionNetwork ######
#
#  Edited by Astrid Reinhardt (Translated from C to perl by K. Nakai)
#
#  Input:   Takes a sequence file in FASTA format
#
#  Output:  Prints the predicted subcellular location ("cytoplasmic" or 
#           "nuclear") and the expected prediction accuracy as determined 
#	    by the reliability index in percent (XX.X) to STDOUT.
#

sub mkaac2 {
    my($start) = $_[0];
        # numbering of $start, $::seq: 0, 1, 2, ..
    my($i, @aac);
    my($len) = $::leng - $start;

    my %aanum2 = (  ## different from %aanum !!

		 'A'=> 0, 'C'=> 1, 'D'=> 2, 'E'=> 3, 'F'=> 4,
		 'G'=> 5, 'H'=> 6, 'I'=> 7, 'K'=> 8, 'L'=> 9,
		 'M'=>10, 'N'=>11, 'P'=>12, 'Q'=>13, 'R'=>14,
		 'S'=>15, 'T'=>16, 'V'=>17, 'W'=>18, 'Y'=>19,

		 );

    foreach $i (0..20) {$aac[$i]=0}

    for($i=$start; $i<$::leng; $i++) {
	my($c) = $::seq[$i];
	
        if(defined($aanum2{$c})) {
            $aac[$aanum2{$c}]++;
        } else {
	    my($j) = $i -1;
            warn "*** mkaac: illegal char $c at position $j: skipped\n";
        }
    }    

#    foreach $i (0..9) {printf " %5.1f", $aac[$i]}; print "\n";
#    foreach $i (10..20) {printf " %5.1f", $aac[$i]}; print "\n";
    foreach $i (0..$#aac) {$aac[$i] = $aac[$i] / $len}
    @aac;
}

sub areinha {
#   $NumberOfInputNodes = 20;
#   $NumberOfOutputNodes = 2;
    my($output) = -1000000;

    my @Sources =  (
		   1,  2,  3,  4,  5, 
		   6,  7,  8,  9, 10, 
		   11, 12, 13, 14, 15, 
		   16, 17, 18, 19, 20, 

		   1,  2,  3,  4,  5, 
		   6,  7,  8,  9, 10, 
		  11, 12, 13, 14, 15, 
		  16, 17, 18, 19, 20, 
		    );

#  /* Weigths definition section */
    my @Weights =  (
		    3.941450, -16.521629, -4.772880, 1.039410, 45.400211, 
		    13.965260, -4.350150, 15.738930, -23.394470, 13.119190, 
		    4.275200, -30.573271, -28.866810, -20.460581, -25.836269,
		    -39.424400, -7.444780, 44.204250, 53.092739, 9.683040, 

		    -3.941440, 16.521629, 4.772880, -1.039410, -45.400219, 
		    -13.965260, 4.350160, -15.738940, 23.394480, -13.119180, 
		    -4.275210, 30.573280, 28.866800, 20.460590, 25.836281, 
		    39.424381, 7.444790, -44.204239, -53.092709, -9.683030, 
		    );

    my @Sources_latter = @Sources[20..39];
    my @Weights_latter = @Weights[20..39];

# /* unit definition section */
#          float act;         /* Activation       */
#          float Bias;        /* Bias of the Unit */
#          int   NoOfSources; /* Number of predecessor units */
#   struct UT   **sources; /* predecessor units */
#          float *weights; /* weights from predecessor units */

    @Units = (
	      { ## unit dummy
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => [], weights => [],
	      },
	      { ## unit 1
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 2
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 3
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 4
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 5
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 6
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 7
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 8
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 9
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 10
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 11
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 12
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 13
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 14
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 15
  		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 16
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 17
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 18
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 19
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 20
		  act    => 0, Bias   => 0, NoOfSources => 0,
		  sources => \@Sources, weights => \@Weights,
	      },
	      { ## unit 21
		  act    => 0, Bias   => 2.100050,
		  NoOfSources => 20,
		  sources => \@Sources,
		  weights => \@Weights,
	      },
	      { ## unit 22
		  act    => 0, Bias   => -2.100050,
		  NoOfSources => 20,
		  sources => \@Sources_latter,
		  weights => \@Weights_latter,
	      }
      );

##
    my(@NetInput) = &mkaac2(0);
    my @NetOutput = &RunNeuralNetworkCytoplasmicNuclear(@NetInput); 

#    print "@NetOutput\n";

#  /*-----------------------------------------------------------
#    determine and print the predicted subcellular location 
#    and calculate the difference between both output nodes 
#    (needed to determine the reliability index).
#  -----------------------------------------------------------*/

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#nncn\">NNCN: Reinhardt's method for Cytplasmic/Nuclear discrimination</A>\n";
    } elsif($::OPT_V) {
	print "NNCN: Reinhardt's method for Cytoplasmic/Nuclear discrimination\n";
    }

    if ($NetOutput[0] > $NetOutput[1]){
	print ' 'x6, "Prediction: cytoplasmic\n" if $::OPT_V;
	$OutPutNodeDifference = $NetOutput[0] - $NetOutput[1];
    } else {
	print ' 'x6, "Prediction: nuclear\n" if $::OPT_V;
	$OutPutNodeDifference = $NetOutput[1] - $NetOutput[0];
    }

#  /*-----------------------------------------------------------
#    dermine the reliability index and print the expected 
#    prediction accuracy associated with it
#  -----------------------------------------------------------*/

    if ($OutPutNodeDifference<0.2){
	$output = 55.5;
    }elsif ($OutPutNodeDifference<0.4) {
	$output = 70.6;
    }elsif ($OutPutNodeDifference<0.6) {
	$output = 76.7;
    }elsif ($OutPutNodeDifference<0.8)  {
	$output = 89.0;
    } else  {
	$output = 94.1;
    }
    print ' 'x6, "Reliability: $output\n\n" if $::OPT_V;
    
    ($NetOutput[1] > $NetOutput[0]) ? $output : (-$output);
}
#####################################

sub Act_Logistic {
    my($sum, $bias) = @_;
    (($sum+$bias)<10000.0)?(1.0/(1.0+exp(-$sum-$bias))):0.0;
}

#/*=====================================================================
#  function "RunNeuralNetworkCytoplasmicNuclear" actually runs the 
#  prediction by the neural network.
#====================================================================*/

sub RunNeuralNetworkCytoplasmicNuclear {

    my(@in) = @_;

#  /* layer definition section (names & member units) */

    my @Input = (
	     1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 
	    11, 12, 13, 14, 15, 16, 17, 18, 19, 20
    );

    my @Output1 = ($Units[21], $Units[22]); #/* members */
    my @Output = (21, 22);

    for($member = 0; $member < 20; $member++) {
	$Units[$Input[$member]]{act} = $in[$member];
    }

    for ($member = 0; $member < 2; $member++) {
	$unitref = $Output1[$member];
	$sum = 0.0;
	for ($source = 0; $source < $$unitref{NoOfSources}; $source++) {
	    $sum +=  $Units[ $$unitref{sources}[$source] ]{act}
	         * $$unitref{weights}[$source];
        }
        $$unitref{act} = &Act_Logistic($sum, $$unitref{Bias});
    };

    for($member = 0; $member < 2; $member++) {
        $out[$member] = $Units[$Output[$member]]{act};
    }

#    print "@out\n";
    @out;

}


1;
