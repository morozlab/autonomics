
############################################################
##  PSORT II:  subroutine package 2
##                        copyright: Kenta Nakai 1996-1998
##                               last update: June 8, 1998
############################################################

## Sequence Numbering convention:
##   argument: 1, 2, 3,...   internal: 0, 1, 2,...

package sub2;

$LSEG2 = 8;

%aanum = (
     'A'=> 0, 'R'=> 1, 'N'=> 2, 'D'=> 3, 'C'=> 4,
     'Q'=> 5, 'E'=> 6, 'G'=> 7, 'H'=> 8, 'I'=> 9,
     'L'=>10, 'K'=>11, 'M'=>12, 'F'=>13, 'P'=>14,
     'S'=>15, 'T'=>16, 'W'=>17, 'Y'=>18, 'V'=>19,
     'B'=>20, 'Z'=>21, 'X'=>22
);


## mitochondrial #############################################

sub mitdisc {
    my($i, $j, $range, $m75, $m95, $limit, $score, $n);
    my(@freq, @iseq);

    foreach $i (0..22) {$freq[$i] = 0}
	$range = &getrange($::seq);
	foreach $j (0..$range-1) {
	    die unless defined $aanum{$::seq[$j]};
	    $n = $aanum{$::seq[$j]};
	    $freq[$n]++;
	    $iseq[$j] = $n;
	}
    $limit = ($range>$LSEG2) ? $LSEG2 : $range;
    splice(@iseq, $limit);
#print join(" ", @iseq), "\n";
    $m95 = &hmom(95, @iseq);
    $m75 = &hmom(75, @iseq);

    my(@a) = 
		( 2.0102 - 1.0973,  ## 0: R
		  0.3922 - 0.2538,  ## 1: M75
		  0.4737 - 0.3566,  ## 2: M95
		 -0.7417 + 0.2655,  ## 3: G 
		  9.2710 -11.3204,  ## 4: AC 
		  0.7522 - 0.4503,  ## 5: ST 
		-17.5993 +13.0901 ); ##6: const

    $score = $a[0] * $freq[1]
		   + $a[1] * $m75
		   + $a[2] * $m95
		   + $a[3] * $freq[7]
		   + $a[4] * ($freq[3]+$freq[6])
		   + $a[5] * ($freq[15]+$freq[16])
		   + $a[6];

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#mit\">MITDISC: discrimination of mitochondrial targeting seq</A>\n";
    } elsif($::OPT_V) {
	print "MITDISC: discrimination of mitochondrial targeting seq\n";
    }
    if($::OPT_V) {
	printf "      R content:   %5d       Hyd Moment(75): %5.2f\n", 
						$freq[1], $m75;
	printf "      Hyd Moment(95): %5.2f    G content:   %5d\n", 
						$m95, $freq[7];
	printf "      D/E content: %5d       S/T content: %5d\n",
					$freq[3]+$freq[6], $freq[15]+$freq[16];
	printf "      Score: %5.2f\n\n", $score;
    }
    $score;
}

sub hmom {
	my($deg, @iseg) = @_;
	my($scr, $i);
	my($ARG) = 3.14159 * $deg / 180;
	my($sums, $sumc) = (0, 0);
	my($len) = scalar(@iseg);

### The GES scale (Engelman et al, 1986) ###
	my(@hyd) = 
#   A/L    R/K    N/M    D/F    C/P    Q/S    E/T    G/W    H/Y    I/V
 ( -6.70, 51.50, 20.10, 38.50, -8.40, 17.20, 34.30, -4.20, 12.60,-13.00,
  -11.70, 36.80,-14.20,-15.50,  0.80, -2.50, -5.00, -7.90,  2.90,-10.90, 0.);

	foreach $i (0..$len-1) {
		$sums += $hyd[$iseg[$i]] * sin($ARG * $i);
		$sumc += $hyd[$iseg[$i]] * cos($ARG * $i);
	}
	$scr = sqrt($sums ** 2 + $sumc ** 2) / $LSEG2;
}

sub getrange {
        my($sq) = $_[0];
#		print $sq, "\n";

        my($POS1) = 15;
        my($pos, $i, $c); ## probable transit peptide 1-$pos
        my($flg) = "not found"; ## acidic residue before $POS1
        my($len) = length($sq);
        if($len < $POS1) {return 0}
        for($i=1; $i<$POS1; $i++) {
                $c = substr($sq, $i, 1);
                if($c eq 'D' || $c eq 'E') {
                        $flg = "found";
                        $pos = $i+1; # from the next residue
                        last;
                }
        }
        if($flg eq "not found") {$pos = $POS1}
        for($i=$pos; $i<$len; $i++) {
                $c = substr($sq, $i, 1);
                if($c eq 'D' || $c eq 'E') {
                        return $i+1; ## pos is counted from 1
                }
        }
        $len;  # if no acidic resudue at all
}

## gavel et al's method
sub gavel {
    my($ipos, $isite, $motif);

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#mit\">Gavel: prediction of cleavage sites for mitochondrial preseq</A>\n";
    } elsif($::OPT_V) {
	print "Gavel: prediction of cleavage sites for mitochondrial preseq\n";
    }

	if(($ipos = &r3) > 0) {
		$isite = $ipos + 4;
		$motif = 'R-3';
		if($::OPT_V) {
			print "      R-3 motif at $isite  ";
			print substr($::seq, $ipos-1, 4),"|",
				substr($::seq, $ipos+3, 1),"\n\n";
		}
	} elsif (($ipos = &r10) > 0) {
		$isite = $ipos + 11;
		$motif = 'R-10';
		if($::OPT_V) {
			print "      R-10 motif at $isite  ";
			print substr($::seq, $ipos-1, 3)," ",
				substr($::seq, $ipos+2, 2),"\n\n";
		}
	} elsif (($ipos = &r2) > 0) {
		$motif = 'R-2';
		$isite = $ipos + 11;
		if($::OPT_V) {
			print "      R-2 motif at $isite  ";
			print substr($::seq, $ipos-1, 3),"|",
					substr($::seq, $ipos+2, 2),"\n\n";
		}
	} else {
		$isite = 0;
		$motif = '___';
		if($::OPT_V) {
			print "      cleavage site motif not found\n\n";
		}
	}
	($isite, $motif);
}

sub abs {
	($_[0]>=0) ? $_[0] : (-1 * $_[0]);
}

## RXY(S/A) motif
sub r3 {
	my($max) = ($::leng < 67) ? ($::leng - 10) : 67;
	my($min) = 11;  ## 12-1
	my($ipos) = -1000;
	foreach $i ($min..$max-3) {  ## $max-1 -> $max-3
		if($::seq[$i] eq 'R' && $::seq[$i+2] eq 'Y'
			&& ($::seq[$i+3] eq 'S' || $::seq[$i+3] eq 'A')) {
			$ipos = $i if &abs($i-23) < &abs($ipos-23);
		}
	}
	$ipos;
}

## RXFS motif
sub r10 {
	my($max) = ($::leng < 60) ? ($::leng-10) : 60;
	my($min) = 4;  ## 5-1
	my($ipos) = -1000;
	foreach $i ($min..$max-3) {  ## $max-1 -> $max-3
		if($::seq[$i] eq 'R' && $::seq[$i+2] eq 'F' && $::seq[$i+3] eq 'S') {
			$ipos = $i if &abs($i-23) < &abs($ipos-23);
		}
	}
	$ipos;
}

## transition point in terms of negative charges
sub r2 {
	my($i, $j, $ipos);
	my($k) = -1000;
	LOOP: foreach $i (0..$::leng-13) {
		if($::seq[$i] eq 'D' || $::seq[$i] eq 'E') {
			foreach $j ($i+1..$i+12) {
				if($::seq[$j] eq 'D' || $::seq[$j] eq 'E') {
					$k = 0;  ### found
					$ipos = $i;
					last LOOP;
				}
			}
		}
	}
	if($k == 0) {
		for ($k = $ipos-2; $k > 0; $k--) {
			last if $::seq[$k] eq 'R';
		}
	}
	return $k;
}

# nuclear ##############################################################

sub nucdisc {
    my($i, $j, $score);

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#nuc\">NUCDISC: discrimination of nuclear localization signals</A>\n";
    } elsif($::OPT_V) {
	print "NUCDISC: discrimination of nuclear localization signals\n";
    }

### new parameter 
#$start = (times)[0];
    $score = 0;
    $score += ( 0.0901 - 0.0274) * &nls1;
    $score += ( 0.1648 - 0.0786) * &nls2;
    $score += ( 0.1063 - 0.0241) * &bipartite;
    $score += ( 1.1162 - 0.1665) * &nucaa;
    $score += (-1.2642 + 0.7904);
#$end = (times)[0];
#printf "(old)%.2f\t", $end - $start;

    printf "      NLS Score: %5.2f\n\n", $score if $::OPT_V;
    $score;
}

# nls1: 4 residue pattern
sub nls1 {
	my($i, $j, $scr, $c, $iscr);
	my($nbas, $np, $nh, $flg) = (0, 0, 0, 'none');
# init value
	foreach $j (0..3) {
		$c = $::seq[$j];
		if($c eq 'R' || $c eq 'K') {$nbas++;}
		elsif($c eq 'P') {$np++;}
		elsif($c eq 'H') {$nh++;}
	}
	if($nbas == 4) {$scr = 5}
	elsif($nbas == 3 && $np == 1) {$scr = 4}
	elsif($nbas == 3 && $nh == 1) {$scr = 3}
	else {$scr =  0}
	if($scr > 0 && $::OPT_V) {
		$flg = 'found';
		printf "      pat1: %s at    1\n", substr($::seq, 0, 4);
	}

	foreach $i (0..$::leng-5) {
		$c = $::seq[$i+4];
		if($c eq 'R' || $c eq 'K') {$nbas++;}
		elsif($c eq 'P') {$np++;}
		elsif($c eq 'H') {$nh++;}
		$c = $::seq[$i];
		if($c eq 'R' || $c eq 'K') {$nbas--;}
		elsif($c eq 'P') {$np--;}
		elsif($c eq 'H') {$nh--;}

		if($nbas == 4) {$iscr = 5}
		elsif($nbas == 3 && $np == 1) {$iscr = 4}
		elsif($nbas == 3 && $nh == 1) {$iscr = 3}
		else {$iscr = 0}
		if($iscr > 0) {
			$scr += $iscr;
			$flg = 'found' if $flg eq 'none';
			printf "      pat4: %s (%1d) at %4d\n", 
				substr($::seq, $i+1, 4), $iscr, $i+2 if $::OPT_V;
		}
	}
	printf "      pat4: none\n" if $::OPT_V && ($flg eq 'none');
	
	$scr;
}

# nls2: (max) 7 residue pattern
sub nls2 {
	my($i, $iscr);
	my($scr, $flg) = (0, 'none');
	foreach $i (0..$::leng-8) {
		if($::seq[$i] eq 'P') {
			$iscr = &scr7($i);
			if($iscr != -1) {
				$scr += $iscr;
				$flg = 'found' if $flg eq 'none';
				printf "      pat7: %s (%d) at %4d\n", 
					substr($::seq, $i, 7), $iscr, $i+1 if $::OPT_V;
			}
		}
	}
	printf "      pat7: none\n" if($::OPT_V && ($flg eq 'none'));
	$scr;
}

sub scr7 {
	my($i) = $_[0];
	return -1 unless $::seq[$i] eq 'P';
	my($max) = -1;
	my($j, $k, $nbas);
	foreach $k (0..2) {
		$nbas = 0;
		foreach $j (1..4) {
			$nbas++ if $::seq[$i+$j+$k]eq'R'||$::seq[$i+$j+$k]eq'K'; 
		}
		if($nbas == 4) {$max = 5}
		elsif($nbas == 3) {$max = ($max > 5-$k) ? $max : (5-$k)}
	}
	$max;
}

# robbins: bipartite NLS (Robbins ,, Dingwall)

sub bipartite {
	my($i, $cnt);
	my($scr, $flg) = (0, 'none');
	foreach $i (0..$::leng-17) {
		next unless (substr($::seq, $i, 2) =~ /[RK][RK]/);
		$cnt = 0;
		foreach $j ($i+12..$i+16) {
			$cnt++ if $::seq[$j] eq 'R' || $::seq[$j] eq 'K';
		}
		$scr += $cnt if $cnt >= 3;
		if($cnt >= 3) {
			$scr += $cnt;
			$flg = 'found' if $flg eq 'none';
			printf "      bipartite: %s at %4d\n", 
				substr($::seq, $i, 17), $i+1 if $::OPT_V;
		}
	}
	printf "      bipartite: none\n" if $::OPT_V && $flg eq 'none';
	$scr;
}

# nucaac
sub nucaa {
	my($cnt) = 0;
	my($i);
	foreach $i (0..$::leng-1) {
		if($::seq[$i]eq'R' || $::seq[$i]eq'K') {
			$cnt++;
		}
	}
	$cnt /= $::leng;
	printf "      content of basic residues: %5.1f%%\n", $cnt * 100
			if $::OPT_V;
	if($cnt > 0.2) {
		int($cnt * 10 - 1);
	} else {
		0;
	}
}

# RNA-binding motif #######################################################

sub rnp1 {
    my(@match);
    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#nuc\">RNA-binding motif:</A>" 
    } elsif($::OPT_V) {
	print "RNA-binding motif:";
    }
    @match = ($::seq =~ /[RK]G[^EDRKHPCG][AGSCI][FY][LIVA].[FYM]/g);
    if($::OPT_V) {
	&report_motif(0, @match);
	print "\n";
    }
    ($#match < 0) ? 0 : ($#match+1);
}

# actinin-type actin-binding ##############################################

sub actin {
    my(@match);
    my($hit) = 0;
    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#motif\">Actinin-type actin-binding motif:</A>\n";
    }elsif($::OPT_V) {
	print "Actinin-type actin-binding motif:\n";
    }
    print "      type 1:" if $::OPT_V;
    @match = ($::seq =~ /[EQ]..[ATV]F..W.N/g);
    if($::OPT_V) {
	&report_motif(0, @match);
    }
    $hit += scalar(@match);
    print "      type 2:" if $::OPT_V;
    @match = ($::seq =~ /[LIVM].[SGN][LIVM][DAGHE][SAG].[DEAG][LIVM].[DEAG]....[LIVM].L[SAG][LIVM][LIVM]W.[LIVM][LIVM]/g);
    if($::OPT_V) {
	&report_motif(0, @match);
	print "\n";
    }
    $hit += scalar(@match);
}

# lumen pf ER ##############################################################

sub hdel {
    my($seg) = substr($::seq, $::leng-4);
    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#er\">KDEL</A>: ER retention motif in the C-terminus: ";
    } elsif($::OPT_V) {
	print "KDEL: ER retention motif in the C-terminus: ";
    }

    if ($seg =~ /HDEL$/ || $seg =~ /KDEL$/) {
	print "$seg\n\n" if $::OPT_V;
	return 1;
    } else {
	print "none\n\n" if $::OPT_V;
	return 0;
    }
}

# peroxisome ##############################################################

sub pts1 {
    my($scr, $pat);
    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#pox\">SKL</A>: peroxisomal targeting signal in the C-terminus: ";
    } elsif($::OPT_V) {
	print "SKL: peroxisomal targeting signal in the C-terminus: ";
    }
		
    if($::seq    =~ /SKL$/) {$pat = 'SKL'; $scr = 10/12}
    elsif($::seq =~ /SKF$/) {$pat = 'SKF'; $scr = 1/2}
    elsif($::seq =~ /AKL$/) {$pat = 'AKL'; $scr = 1/2}
    elsif($::seq =~ /([SAGCN][RKH][LIVMAF])$/) {$pat = $1; $scr = 1/4}
    else {$scr = 0}
    if($::OPT_V) {
	if($scr > 0) {
	    print "$pat\n\n";
	} else {
	    print "none\n\n";
	}
    }
    $scr;
}

sub pts2 {
    my($scr, @match);
    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#pox\">SKL2</A>: 2nd peroxisomal targeting signal: ";
    } elsif($::OPT_V) {
	print "PTS2: 2nd peroxisomal targeting signal: ";
    }
    @match = ($::seq =~ /[RK][LI].....[HQ]L/g);
    if($::OPT_V) {
	&report_motif(0, @match);
	print "\n";
    }
    ($#match < 0) ? 0 : ($#match+1);
}

# vacuolar ##############################################################

sub vaccalc {
    my(@match);
    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#lys\">VAC</A>: possible vacuolar targeting motif:";
    } elsif($::OPT_V) {
	print "VAC: possible vacuolar targeting motif:";
    }
#    print "" if $::OPT_V;
    @match = ($::seq =~ /[TIK]LP[NKI]/g);
    if($::OPT_V) {
	&report_motif(0, @match);
	print "\n";
    }
    ($#match < 0) ? 0 : ($#match+1);
}

# N-myristoyl ############################################################

sub nmyr {
    my($pat, $i);
    my($scr, $num_k) = (0, 0);

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#anc\">NMYR</A>: N-myristoylation pattern : ";
    }elsif($::OPT_V) {
	print "NMYR: N-myristoylation pattern : ";
    }

    if($::seq =~ /^(M?G[^EDRKHPFYW]..[STAGCN][^P])/) {
	$pat = $1;
	$scr++;
	print $pat, "\n" if $::OPT_V;
	if($::seq =~ /^M?GC/) {
	    $scr++;
	    print "      3rd aa is cysteine (may be palmitylated)\n"
				if $::OPT_V;
	}
	foreach $i (3..9) {
	    $num_k++ if $::seq[$i] eq 'K'
	}
	if($num_k >= 2) {
	    $scr++;
	    print "      additional alternating lysine motif\n" if $::OPT_V;
	}
	print "\n" if $::OPT_V;
    } else {
	print "none\n\n" if $::OPT_V;
    }
    $scr;
}

# isoprenyl ##############################################################

sub isoprenyl {
    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#anc\">Prenylation motif:</A> ";
    }elsif($::OPT_V) {
	print "Prenylation motif: ";
    }
    if($::seq =~ /(C[^DENQ][LIVM].)$/) {
	if($::OPT_V) {
	    print "      CaaX motif in the C-terminus: $1\n";
	    print "         if X is S, A, or M, it will be farnesylated\n";
	    print "         otherwise, it will be geranylgeranylated\n\n";
	}
	2;
    } elsif($::seq =~ /(C.C)$/) {
	print "      CXC motif in the C-terminus: $1\n\n" if $::OPT_V;
	1;
    } elsif($::seq =~ /(CC..)$/) {
	print "      CC motif near the C-terminus: $1\n\n" if $::OPT_V;
	1;
    } else {
	print "none\n\n" if $::OPT_V;
	0;
    }
}

# ER membrane ##############################################################

sub erm {
    my($mtype) = $_[0];
    my($count, $scr) = (0, 0);
    my(@pat);
    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#er\">ER Membrane Retention Signals:</A> ";
    } elsif($::OPT_V) {
	print "ER Membrane Retention Signals: ";
    }
    my($nseg, $cseg) = (substr($::seq, 1, 4), substr($::seq, $::leng-5, 4));

    @pat = ($nseg =~ /R/g);
    $count = scalar(@pat);
    if($count != 0) {
	print "\n      XXRR-like motif in the N-terminus: $nseg\n\n"
	    if $::OPT_V;
	$scr = $count;
	$scr += 2 if $mtype eq '2 ' || $mtype eq 'Nt';
    }

    @pat = ($cseg =~ /K/g);
    $count = scalar(@pat);
    if($count != 0) {
	print "\n      KKXX-like motif in the C-terminus: $cseg\n\n"
	    if $::OPT_V;
	$scr += $count;
	$scr += 2 if $mtype eq '1a' || $mtype eq '1b';
    } else {
	print "none\n\n" if $::OPT_V;
    }
    $scr;
}

# YQRL motif ##############################################################

sub yqrl {
    my($tmsnum, $start) = @_;

    my(@match);
    my($scr) = 0;

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#ves\">memYQRL</A>: transport motif from cell surface to Golgi:";
    } elsif($::OPT_V) {
	print "memYQRL: transport motif from cell surface to Golgi:";
    }
    if($tmsnum == 0) {
	print " none\n\n" if $::OPT_V;
    } elsif($tmsnum == 1) {
#	@match = (substr($::seq, $start-1, $end-$start+1) =~ /YQRL/g);
	@match = ($::seq =~ /YQRL/g);
	$scr += scalar(@match) + 3;
	if($::OPT_V) {
	    &report_motif($start-1, @match);
	    print "\n";
	}
    } else {
	@match = ($::seq =~ /YQRL/g);
	$scr += scalar(@match);
	if($::OPT_V) {
	    &report_motif($start-1, @match);
	    print "\n";
	}
    }
    $scr;
}

# tyrosine(s) in the tail #################################################

sub tyros {
    my($tmsnum, $start, $end) = @_;
    my(@ylist) = ();
    my($i);
    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#motif\">Tyrosines in the tail:</A>";
    } elsif($::OPT_V) {
	print "Tyrosines in the tail:";
    }
    if($tmsnum != 1) {
	print " none\n\n" if $::OPT_V;
	return 0;
    } elsif($end - $start + 1 > 50) {
	print " too long tail\n\n" if $::OPT_V;
	return 0;
    } else {
	for($i=$start-1; $i<=$end-1; $i++) {
	    push(@ylist, $i) if $::seq[$i] eq 'Y';
	}
	if(scalar(@ylist) == 0) {
	    print " none\n\n" if $::OPT_V;
	    return 0;
	}
	print join(",", @ylist), "\n\n" if $::OPT_V;
	return 10 * scalar(@ylist) / ($end - $start + 1);
    }
}

# dileucine motif ########################################################

sub dileu {
    my($tmsnum, $start, $end) = @_;
    my(@match);
    my($scr) = 0;
    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#motif\">Dileucine motif in the tail:</A>";
    } elsif($::OPT_V) {
	print "Dileucine motif in the tail:";
    }
    if($tmsnum != 1) {
	print " none\n\n" if $::OPT_V;
    } else {
	@match = (substr($::seq, $start-1, $end-$start+1) =~ /LL/g);
	$scr = 10 * scalar(@match) / ($end - $start + 1);
	if($::OPT_V) {
	    &report_motif($start-1, @match);
	    print "\n";
	}
    }
    $scr;
}

##########################################################################
## must be called when $::OPT_V
sub report_motif {
	my($init, @pat) = @_;
##	print "### ", join(" ", @pat), "\n";
	if($#pat < 0) {
		print " none\n";
		return;
	}
	my($pos);
	print " found\n";
	foreach (@pat) {
		$pos = index($::seq, $_, $init);
		$init = $pos + 1;
		printf "      %s at %d\n", $_, $pos+1;
	}
}

# coilded coil #########################################################
# 
#### A. Lupas's algorithm to detect coiled-coil regions
##           it returns the number of residues predicted to be in 
##           the coiled-coil structure

sub lupas {

    my $WINLEN = 28;

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#lupas\">COIL: Lupas's algorithm to detect coiled-coil regions\n</A>";
    } elsif($::OPT_V) {
	print "COIL: Lupas's algorithm to detect coiled-coil regions\n";
    }

    if($::leng < $WINLEN) {
	print "      length: $::leng too short; skipped\n\n" if $::OPT_V;
	return 0;
    }

    my %param = (
	L => [3.167, 0.297, 0.398, 3.902, 0.585, 0.501, 0.483],
	I => [2.597, 0.098, 0.345, 0.894, 0.514, 0.471, 0.431],
        V => [1.665, 0.403, 0.386, 0.949, 0.211, 0.342, 0.360],
        M => [2.240, 0.370, 0.480, 1.409, 0.541, 0.772, 0.663],
        F => [0.531, 0.076, 0.403, 0.662, 0.189, 0.106, 0.013],
        Y => [1.417, 0.090, 0.122, 1.659, 0.190, 0.130, 0.155],
        G => [0.045, 0.275, 0.578, 0.216, 0.211, 0.426, 0.156],
        A => [1.297, 1.551, 1.084, 2.612, 0.377, 1.248, 0.877],
        K => [1.375, 2.639, 1.763, 0.191, 1.815, 1.961, 2.795],
        R => [0.659, 1.163, 1.210, 0.031, 1.358, 1.937, 1.798],
        H => [0.347, 0.275, 0.679, 0.395, 0.294, 0.579, 0.213],
        E => [0.262, 3.496, 3.108, 0.998, 5.685, 2.494, 3.048],
        D => [0.030, 2.352, 2.268, 0.237, 0.663, 1.620, 1.448],
        Q => [0.179, 2.114, 1.778, 0.631, 2.550, 1.578, 2.526],
        N => [0.835, 1.475, 1.534, 0.039, 1.722, 2.456, 2.280],
        S => [0.382, 0.583, 1.052, 0.419, 0.525, 0.916, 0.628],
        T => [0.169, 0.702, 0.955, 0.654, 0.791, 0.843, 0.647],
        C => [0.824, 0.022, 0.308, 0.152, 0.180, 0.156, 0.044],
        W => [0.240, 0.,    0.,    0.456, 0.019, 0.,    0.   ],
        P => [0.,    0.008, 0.,    0.013, 0.,    0.,    0.   ],
      );

    my($product, $bias, @max);

    foreach $i (0..$::leng-1) {
	$max[$i] = -1;
#	$maxbias[$i] = -1;
    }

    foreach $bias (0..6) {
	$product = $param{$::seq[0]}[$bias];
	foreach $j (1..$WINLEN-1) {
	    $product *= $param{$::seq[$j]}[($bias+$j)%7];
	}
	foreach $j (0..$WINLEN-1) {
	    if($product > $max[$j]) {
		$max[$j] = $product;
#		$maxbias[$j] = $bias;
	    }
	}
	foreach $i (1..$::leng-$WINLEN) {
	    if(($p = $param{$::seq[$i-1]}[($bias+$i-1)%7]) != 0) {
		$product /= $p;
	    } else {
		$product = 1;
		foreach $k (1..$WINLEN-1) {
		    $product *= $param{$::seq[$i-1+$k]}[($bias+$i-1+$k)%7];
		}
	    }
	    $product *= $param{$::seq[$i-1+$WINLEN]}[($bias+$i-1+$WINLEN)%7];
	    foreach $j (0..$WINLEN-1) {
		if($product > $max[$i+$j]) {
		    $max[$i+$j] = $product;
#		    $maxbias[$i+$j] = $bias;
		}
	    }
	}
    }

    my $count = 0;    
    foreach $i (0..$::leng-1) {
	if($max[$i] > 4750) {
	    $val = $max[$i] ** (1/28);
#	    printf "%4d  %s  %5.4f  %1d\n", 
#	    $i, $seq[$i], &gauss($val), $maxbias[$i];
	    printf "%4d %s  %4.2f\n", $i, $::seq[$i], &gauss($val) if $::OPT_V;
	    $count++;
	}
    }
    print "      total: $count residues \n\n" if $::OPT_V;
    $count;
}

sub gauss {
    my $val = $_[0];
    return 0 if $val == 0;

    my $gcc = &gaussian($val, 1.63, 0.24);
    my $gg  = &gaussian($val, 0.77, 0.20);

    $gcc / (30 * $gg + $gcc);
}

sub gaussian {
    my ($score, $mean, $sigma) = @_;

    my $a = $sigma * 2 * 3.141592653589;
    my $b = -0.5 * (($score - $mean)/$sigma) ** 2;
    (exp $b) / $a;
}

	
1;
