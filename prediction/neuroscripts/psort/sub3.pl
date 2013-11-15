
############################################################
##  PSORT II:  subroutine package 3
##                       copyright: Kenta Nakai 1996, 1997
##                       last revised: Nov. 12, 1997
############################################################

## Sequence Numbering convention:
##   argument: 1, 2, 3,...   internal: 0, 1, 2,...

package sub3;

##########################################################################
## must be called when $::OPT_V
sub report_motif2 {
	my($msg, @pat) = @_;
	my($pos);
	my($init) = 0;
	if($#pat >= 0) {
		print "\n   $msg *** found ***\n" if $::OPT_V;
		foreach (@pat) {
			$pos = index($::seq, $_, $init);
			$init = $pos + 1;
			printf "         %s at %d\n", $_, $pos+1 if $::OPT_V;
		}
	}
	scalar(@pat);
}

sub dnabind {
    my($pat, @pat);
    my($hit) = 0;

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#motif\">checking 63 PROSITE DNA binding motifs:</A> ";
    } elsif($::OPT_V) {
	print "checking 63 PROSITE DNA binding motifs: ";
    }

## 'Homeobox' domain signature (PS00027)
	$msg = "   'Homeobox' domain signature (PS00027): " if $::OPT_V;
	$pat = "[LIVMFYG][ASLVR].{2}[LIVSTACN].[LIVM].{4}[LIV][RKNQESTAIY][LIVFSTNKH]W[FYVC].[NDQTAH].{5}[RKNAIMW]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## 'Homeobox' antennapedia-type protein signature (PS00032)
	$msg = "   'Homeobox' antennapedia-type protein signature (PS00032): " if $::OPT_V;
	$pat = "[LIVMFE][FY]PWM[KRQTA]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## 'Homeobox' engrailed-type protein signature (PS00033)
	$msg = "   'Homeobox' engrailed-type protein signature (PS00033): " if $::OPT_V;
	$pat = "LMAQGLYN";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## 'Paired box' domain signature (PS00034)
	$msg = "   'Paired box' domain signature (PS00034): " if $::OPT_V;
	$pat = "RPC.{11}CVS";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## 'POU' domain signature 1 (PS00035)
	$msg = "   'POU' domain signature 1 (PS00035): " if $::OPT_V;
	$pat = "[RKQ]R[LIM].[LF]G[LIVMFY].Q.[DNQ]VG";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## 'POU' domain signature 2 (PS00465)
	$msg = "   'POU' domain signature 2 (PS00465): " if $::OPT_V;
	$pat = "SQ[ST][TA]I[SC]RFE.[LS].[LI][ST]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Zinc finger, C2H2 type, domain (PS00028)
	$msg = "   Zinc finger, C2H2 type, domain (PS00028): " if $::OPT_V;
	$pat = "C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Zinc finger, C3HC4 type, signature (PS00518)
	$msg = "   Zinc finger, C3HC4 type, signature (PS00518): " if $::OPT_V;
	$pat = "C.H.[LIVMFY]C.{2}C[LIVMYA]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Nuclear hormones receptors DNA-binding region signature (PS00031)
	$msg = "   Nuclear hormones receptors DNA-binding region signature (PS00031): " if $::OPT_V;
	$pat = "C.{2}C.[DE].{5}H[FY].{4}C.{2}C.{2}FF.R";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## GATA-type zinc finger domain (PS00344)
	$msg = "   GATA-type zinc finger domain (PS00344): " if $::OPT_V;
	$pat = "C.NC.{4}T.LWR[RK].{3}G.{3}CNAC";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Poly(ADP-ribose) polymerase zinc finger domain (PS00347)
	$msg = "   Poly(ADP-ribose) polymerase zinc finger domain (PS00347): " if $::OPT_V;
	$pat = "C[KR].C.{3}I.K.{3}[RG].{16,18}W[FYH]H.{2}C";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Fungal Zn(2)-Cys(6) binuclear cluster domain (PS00463)
	$msg = "   Fungal Zn(2)-Cys(6) binuclear cluster domain (PS00463): " if $::OPT_V;
	$pat = "[GASTPV]C.{2}C[RKHSTACW].{2}[RKH].{2}C.{5,9}C.{2}C.{6,8}C";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Copper-fist domain (PS01119)
	$msg = "   Copper-fist domain (PS01119): " if $::OPT_V;
	$pat = "M[LIVMF]{3}.{3}K[MY]AC.{2}CI[KR].H[KR].{3}C.H.{8}[KR].[KR]GRP";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

########### SKIP-FLAG: True
## Leucine zipper pattern (PS00029)
	$msg = "   Leucine zipper pattern (PS00029): " if $::OPT_V;
	$pat = "L.{6}L.{6}L.{6}L";
	@pat = $::seq =~ /$pat/go;
	if(&report_motif2($msg, @pat) > 0 && $OPT_V) {
		print "*** this ambiguous motif was not counted as hit\n";
	}

## bZIP transcription factors basic domain signature (PS00036)
	$msg = "   bZIP transcription factors basic domain signature (PS00036): " if $::OPT_V;
	$pat = "[KR].{1,3}[RKSAQ]N.{2}[SAQ]{2}.[RKTAENQ].R.[RK]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Myb DNA-binding domain repeat signature 1 (PS00037)
	$msg = "   Myb DNA-binding domain repeat signature 1 (PS00037): " if $::OPT_V;
	$pat = "W[ST].{2}E[DE].{2}[LIV]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Myb DNA-binding domain repeat signature 2 (PS00334)
	$msg = "   Myb DNA-binding domain repeat signature 2 (PS00334): " if $::OPT_V;
	$pat = "W.{2}[LI][SAG].{4,5}R.{8}[YW].{3}[LIVM]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Myc-type, 'helix-loop-helix' dimerization domain signature (PS00038)
	$msg = "   Myc-type, 'helix-loop-helix' dimerization domain signature (PS00038): " if $::OPT_V;
	$pat = "[DENSTAP]K[LIVMWAGSN][^FYWCPHKR][LIVT][LIV].{2}[STAV][LIVMSTAC].[VMFYH][LIVMTA][^P][^P][LIVMSR]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## p53 tumor antigen signature (PS00348)
	$msg = "   p53 tumor antigen signature (PS00348): " if $::OPT_V;
	$pat = "MCNSSCMGGMNRRP";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## CBF-A/NF-YB subunit signature (PS00685)
	$msg = "   CBF-A/NF-YB subunit signature (PS00685): " if $::OPT_V;
	$pat = "CVSE.ISF[LIVM]T[SG]EA[SC][DE][KRQ]C";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## CBF-B/NF-YA subunit signature (PS00686)
	$msg = "   CBF-B/NF-YA subunit signature (PS00686): " if $::OPT_V;
	$pat = "YVNAKQY.RILKRR.ARAKLE";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## 'Cold-shock' DNA-binding domain signature (PS00352)
	$msg = "   'Cold-shock' DNA-binding domain signature (PS00352): " if $::OPT_V;
	$pat = "[FY]GFI.{6,7}[DER][LIVM]F.H.[STK].[LIVMFY]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## CTF/NF-I signature (PS00349)
	$msg = "   CTF/NF-I signature (PS00349): " if $::OPT_V;
	$pat = "RKRKYFKKHEKR";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ets-domain signature 1 (PS00345)
	$msg = "   Ets-domain signature 1 (PS00345): " if $::OPT_V;
	$pat = "L[FYW][QEDH]F[LI][LVQK].[LI]L";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ets-domain signature 2 (PS00346)
	$msg = "   Ets-domain signature 2 (PS00346): " if $::OPT_V;
	$pat = "[RKH].{2}M.Y[DENQ].[LIVM][STAG]R[STAG][LI]R.Y";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Fork head domain signature 1 (PS00657)
	$msg = "   Fork head domain signature 1 (PS00657): " if $::OPT_V;
	$pat = "KP[PT][FYVQH]S[FY].{2}[LIVM].{3}[AC]I";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Fork head domain signature 2 (PS00658)
	$msg = "   Fork head domain signature 2 (PS00658): " if $::OPT_V;
	$pat = "W[QK][NS]S[LIV]RH";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## HSF-type DNA-binding domain signature (PS00434)
	$msg = "   HSF-type DNA-binding domain signature (PS00434): " if $::OPT_V;
	$pat = "L.{3}[FY]KH.N.[STA]SF[LIVM]RQLN.Y.[FYW][RKH]K[LIVM]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## IRF family signature (PS00601)
	$msg = "   IRF family signature (PS00601): " if $::OPT_V;
	$pat = "P.{2}WK[TA].{2}RCA[LIVM]N.{6}EV.[DE]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## LIM domain signature (PS00478)
	$msg = "   LIM domain signature (PS00478): " if $::OPT_V;
	$pat = "C.{2}C.{15,21}[FYWH]H.{2}C.{2}C.{2}C.{3}[LIVMF]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## SRF-type transcription factors DNA-binding and dimerization domain (PS00350)
	$msg = "   SRF-type transcription factors DNA-binding and dimerization domain (PS00350): " if $::OPT_V;
	$pat = "R.[RK].{5}I.[DN].{3}R.{2}T[FY].[RK]{3}.{2}[LIVM].K{2}A.E[LIVM][ST].L.{4}[LIVM].[LIVM]{3}.{6}[LIVMF].{2}[FY]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## TEA domain signature (PS00554)
	$msg = "   TEA domain signature (PS00554): " if $::OPT_V;
	$pat = "GRNELI.{2}YI.{3}T.{3}RT[RK]{2}Q[LIVM]SSH[LIVM]QV";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Transcription factor TFIIB repeat signature (PS00782)
	$msg = "   Transcription factor TFIIB repeat signature (PS00782): " if $::OPT_V;
	$pat = "G[KR].{3}[STAGN].[LIVMA][GSTA]{2}[CSAV][LIVM][LIVMFY][LIVMA][GSA][STAC]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Transcription factor TFIID repeat signature (PS00351)
	$msg = "   Transcription factor TFIID repeat signature (PS00351): " if $::OPT_V;
	$pat = "Y.P.{2}[IF].{2}[LIVM]{2}.[KR].{3}P[RKQ].{3}L[LIVM]F.[STN]G[KR][LIVM]{2}.{2}G[TAG][KR].{7}[AGC].{6}[LIVM]{2}";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## TFIIS zinc ribbon domain signature (PS00466)
	$msg = "   TFIIS zinc ribbon domain signature (PS00466): " if $::OPT_V;
	$pat = "C.{2}C.{9}[LIVM]QTR[STA].DEP.{6}C.{2}C";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## DEAD-box subfamily ATP-dependent helicases signature (PS00039)
	$msg = "   DEAD-box subfamily ATP-dependent helicases signature (PS00039): " if $::OPT_V;
	$pat = "[LIVMF]{2}DEAD[RKEN].[LIVMFYGSTN]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## DEAH-box subfamily ATP-dependent helicases signature (PS00690)
	$msg = "   DEAH-box subfamily ATP-dependent helicases signature (PS00690): " if $::OPT_V;
	$pat = "[GSAH].[LIVMF]{3}DE[ALIV]H[NECR]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Eukaryotic putative RNA-binding region RNP-1 signature (PS00030)
#	$msg = "   Eukaryotic putative RNA-binding region RNP-1 signature (PS00030): " if $::OPT_V;
#	$pat = "[RK]G[^EDRKHPCG][AGSCI][FY][LIVA].[FYM]";
#	@pat = $::seq =~ /$pat/go;
#	$hit += &report_motif2($msg, @pat);

## Fibrillarin signature (PS00566)
	$msg = "   Fibrillarin signature (PS00566): " if $::OPT_V;
	$pat = "[GST][LIVMP]VYAVEF";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## MCM2/3/5 family signature (PS00847)
	$msg = "   MCM2/3/5 family signature (PS00847): " if $::OPT_V;
	$pat = "IDEFDKM";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## XPA protein signature 1 (PS00752)
	$msg = "   XPA protein signature 1 (PS00752): " if $::OPT_V;
	$pat = "C.[DE]C.{3}[LIVMF].{1,2}D.{2}L.{3}F.{4}C.{2}C";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## XPA protein signature 2 (PS00753)
	$msg = "   XPA protein signature 2 (PS00753): " if $::OPT_V;
	$pat = "[LIVM]{2}T[KR]TE.K.[DE]Y[LIVMF]{2}.D.[DE]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## XPG protein signature 1 (PS00841)
	$msg = "   XPG protein signature 1 (PS00841): " if $::OPT_V;
	$pat = "I[KR]P.[FY]VFDG.{2}P.LK";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## XPG protein signature 2 (PS00842)
	$msg = "   XPG protein signature 2 (PS00842): " if $::OPT_V;
	$pat = "G[LIVM]P[FY][LIVM].AP.EAEA[QS]C";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Single-strand binding protein family signature 1 (PS00735)
	$msg = "   Single-strand binding protein family signature 1 (PS00735): " if $::OPT_V;
	$pat = "[LIVMF][NS][KR][LIVM].[LIVM]{2}G[NHR][LIVM][GT].[DE]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Single-strand binding protein family signature 2 (PS00736)
	$msg = "   Single-strand binding protein family signature 2 (PS00736): " if $::OPT_V;
	$pat = "T.W[HY][RNS][LIVM].[LIVMF][FY][NGR]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Histone H2A signature (PS00046)
	$msg = "   Histone H2A signature (PS00046): " if $::OPT_V;
	$pat = "[AC]GL.FPV";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Histone H2B signature (PS00357)
	$msg = "   Histone H2B signature (PS00357): " if $::OPT_V;
	$pat = "[KR]E[LIVM][EQ]T.{2}[KR].[LIVM]{2}.[PAG][DE]L.[KR]HA[LIVM][STA]EG";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Histone H3 signature 1 (PS00322)
	$msg = "   Histone H3 signature 1 (PS00322): " if $::OPT_V;
	$pat = "KAPRKQL";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Histone H3 signature 2 (PS00959)
	$msg = "   Histone H3 signature 2 (PS00959): " if $::OPT_V;
	$pat = "PF.[RA]LV[KRQ][DEG][IV]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Histone H4 signature (PS00047)
	$msg = "   Histone H4 signature (PS00047): " if $::OPT_V;
	$pat = "GAKRH";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## HMG1/2 signature (PS00353)
	$msg = "   HMG1/2 signature (PS00353): " if $::OPT_V;
	$pat = "[FI]S[KR]KCS[EK]RWKTM";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## HMG-I and HMG-Y DNA-binding domain (A+T-hook) (PS00354)
	$msg = "   HMG-I and HMG-Y DNA-binding domain (A+T-hook) (PS00354): " if $::OPT_V;
	$pat = "[AT].{1,2}[RK]{2}[GP]RGRP[RK].";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## HMG14 and HMG17 signature (PS00355)
	$msg = "   HMG14 and HMG17 signature (PS00355): " if $::OPT_V;
	$pat = "RRSARLSA[RK]P";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bromodomain signature (PS00633)
	$msg = "   Bromodomain signature (PS00633): " if $::OPT_V;
	$pat = "[STANF].{2}F.{4}[DNS].{5,7}[DENQTF]Y[HFY].{2}[LIVMFY].{3}[LIVM].{4}[LIVM].{6,8}Y.{12,13}[LIVM].{2}N[SAC].{2}[FY]N";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bromodomain profile (PS50014)
	$msg = "   Bromodomain profile (PS50014): " if $::OPT_V;
	$pat = "[FYL].[LIVMC][KR]W.[GDNR][FYWLE].{5,6}[ST]W[ES][PSTDN].{3}[LIVMC]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Chromo and chromo shadow domain profile (PS50013)
	$msg = "   Chromo and chromo shadow domain profile (PS50013): " if $::OPT_V;
	$pat = "G.ND.{2}ALGR.T";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Regulator of chromosome condensation (RCC1) signature 2 (PS00626)
	$msg = "   Regulator of chromosome condensation (RCC1) signature 2 (PS00626): " if $::OPT_V;
	$pat = "[LIVMFA][STAGC]{2}G.{2}H[STAGLI][LIVMFA].[LIVM]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Protamine P1 signature (PS00048)
	$msg = "   Protamine P1 signature (PS00048): " if $::OPT_V;
	$pat = "[AV]R[NFY]R.{2,3}[ST].S.S";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Nuclear transition protein 1 signature (PS00541)
	$msg = "   Nuclear transition protein 1 signature (PS00541): " if $::OPT_V;
	$pat = "SKRKYRK";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Nuclear transition protein 2 signature 1 (PS00970)
	$msg = "   Nuclear transition protein 2 signature 1 (PS00970): " if $::OPT_V;
	$pat = "H.{3}HS[NS]S.PQS";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Nuclear transition protein 2 signature 2 (PS00971)
	$msg = "   Nuclear transition protein 2 signature 2 (PS00971): " if $::OPT_V;
	$pat = "K.RK.{2}EGK.{2}K[KR]K";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## DNA mismatch repair proteins mutL / hexB / PMS1 signature (PS00058)
	$msg = "   DNA mismatch repair proteins mutL / hexB / PMS1 signature (PS00058): " if $::OPT_V;
	$pat = "GFRGEAL";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## DNA mismatch repair proteins mutS family signature (PS00486)
	$msg = "   DNA mismatch repair proteins mutS family signature (PS00486): " if $::OPT_V;
	$pat = "[ST][LIVM].[LIVM].DE[LIVM][GC][RH]G[ST].{4}G";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## mutT domain signature (PS00893)
	$msg = "   mutT domain signature (PS00893): " if $::OPT_V;
	$pat = "G.{5}E.{4}[STAGC][LIVMA].RE[LIVMF].EE";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

	if($::OPT_V) {
		print $hit == 0 ? " none\n\n" : "\n";
	}
	$hit;
}


sub ribosomal {
    my($pat, @pat);
    my($hit) = 0;

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#motif\">checking 71 PROSITE ribosomal protein motifs:</A> ";
    } elsif($::OPT_V) {
	print "checking 71 PROSITE ribosomal protein motifs: ";
    }

## Ribosomal protein L2 signature (PS00467)
	$msg = "   Ribosomal protein L2 signature (PS00467): " if $::OPT_V;
	$pat = "P.{2}RG[STAIV]{2}.N[APK].[DE]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L3 signature (PS00474)
	$msg = "   Ribosomal protein L3 signature (PS00474): " if $::OPT_V;
	$pat = "F.{6}[DN].{2}[AGS].[ST].G[KRH]G.{2}G.{3}R";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L5 signature (PS00358)
	$msg = "   Ribosomal protein L5 signature (PS00358): " if $::OPT_V;
	$pat = "[LIVM].{2}[LIVM][STAC][GE][QV].{2}[LIVMA].[STC].[STAG][KR].[STA]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L6 signature 1 (PS00525)
	$msg = "   Ribosomal protein L6 signature 1 (PS00525): " if $::OPT_V;
	$pat = "[PS][DENS].YK[GA]KG[LIVM]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L6 signature 2 (PS00700)
	$msg = "   Ribosomal protein L6 signature 2 (PS00700): " if $::OPT_V;
	$pat = "Q.{3}[LIVM].{2}[KR].{2}R.F.DG[LIVM]Y[LIVM].{2}[KR]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L9 signature (PS00651)
	$msg = "   Ribosomal protein L9 signature (PS00651): " if $::OPT_V;
	$pat = "G.{2}G.{4}V.{2}G[FY].{2}N[FY]L.{3}G.[GA].{3}[STN]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L10 signature (PS01109)
	$msg = "   Ribosomal protein L10 signature (PS01109): " if $::OPT_V;
	$pat = "[DEH].{2}G[LIVMF][STN][VA].[DEQK][LIVMA].{2}[LM]R";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L11 signature (PS00359)
	$msg = "   Ribosomal protein L11 signature (PS00359): " if $::OPT_V;
	$pat = "[RKN].[LIVM].G[ST].{2}[SNQ][LIVM]G.{2}[LIVM].{0,1}[DENG]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L13 signature (PS00783)
	$msg = "   Ribosomal protein L13 signature (PS00783): " if $::OPT_V;
	$pat = "[LIVM][KRV]GM[LV][PS].{4,5}[GS][QEKR].{5}[LIVM].[AV][LFY].G";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L14 signature (PS00049)
	$msg = "   Ribosomal protein L14 signature (PS00049): " if $::OPT_V;
	$pat = "A[LIV]{3}.{9,10}[DNS]G.{4}[FY].{2}[NT].{2}V[LIV]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L15 signature (PS00475)
	$msg = "   Ribosomal protein L15 signature (PS00475): " if $::OPT_V;
	$pat = "K[LIVM]{2}.{2}G.[LIVM].{3,4}[LIVM].[LIVMF].{4}[LIVMF][ST].{2}A.{3}[LIVM].{3}G";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L16 signature 1 (PS00586)
	$msg = "   Ribosomal protein L16 signature 1 (PS00586): " if $::OPT_V;
	$pat = "[KR]R.[GSC][KQV][LIVM]W[LIVM][KR][LIVM][LFY][AP]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L16 signature 2 (PS00701)
	$msg = "   Ribosomal protein L16 signature 2 (PS00701): " if $::OPT_V;
	$pat = "RMG.[GR]KG.{4}[FWK]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L19 signature (PS01015)
	$msg = "   Ribosomal protein L19 signature (PS01015): " if $::OPT_V;
	$pat = "R[KR]G.VR[KR]AKLYYLR[ED]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L20 signature (PS00937)
	$msg = "   Ribosomal protein L20 signature (PS00937): " if $::OPT_V;
	$pat = "K.{3}[KR].[LIVM]WI[STNALV]R[LIVM]N.{3}[RK]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L22 signature (PS00464)
	$msg = "   Ribosomal protein L22 signature (PS00464): " if $::OPT_V;
	$pat = "[RKQN].{4}[RH][GAS].G[KRQS].{9}[HDN][LIVM].[LIVMS].[LIVM]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L23 signature (PS00050)
	$msg = "   Ribosomal protein L23 signature (PS00050): " if $::OPT_V;
	$pat = "[RK]{2}[AM][IVFYT][IV][RKT]L[STANQK].{7}[LIVMFT]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L24 signature (PS01108)
	$msg = "   Ribosomal protein L24 signature (PS01108): " if $::OPT_V;
	$pat = "V.[LIVM][LIVMA][STA][GS].{2}[KR].{3}G.[LIVMA]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L27 signature (PS00831)
	$msg = "   Ribosomal protein L27 signature (PS00831): " if $::OPT_V;
	$pat = "G.[LIVM]{2}.RQRG.{5}G";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L29 signature (PS00579)
	$msg = "   Ribosomal protein L29 signature (PS00579): " if $::OPT_V;
	$pat = "[KRNQS].{3}[LIVMFYA][KRS].[LIVMFYT][KR]{2}[DESTAN][LIVM]A[RC]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L30 signature (PS00634)
	$msg = "   Ribosomal protein L30 signature (PS00634): " if $::OPT_V;
	$pat = "[LIVM][LIVMT].{2}[RKAP][GS].{6,8}[KDN].{3}[LIVMT][LIVM].{2}[LIVMF].[LIVM].{2}[LIVMT].{3}[LIVMT]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L31 signature (PS01143)
	$msg = "   Ribosomal protein L31 signature (PS01143): " if $::OPT_V;
	$pat = "C[SG].CHPF[FY]TG[KR]Q[KR]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L33 signature (PS00582)
	$msg = "   Ribosomal protein L33 signature (PS00582): " if $::OPT_V;
	$pat = "Y.[ST].[KR][NS][KRS].{3}[PA].{1,2}[LIVM]E.{2}K[FY]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L34 signature (PS00784)
	$msg = "   Ribosomal protein L34 signature (PS00784): " if $::OPT_V;
	$pat = "KRT[FYW]QP.{3}[KR][RH].[KR].{2}GF.{2}RM";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L35 signature (PS00936)
	$msg = "   Ribosomal protein L35 signature (PS00936): " if $::OPT_V;
	$pat = "K[LIVM]KT.{2}[GA].{2}KR[LIVMFY][KR]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L36 signature (PS00828)
	$msg = "   Ribosomal protein L36 signature (PS00828): " if $::OPT_V;
	$pat = "C.{2}C.{2}[LIVM].R.{3}[LIVMN].[LIVM].C.{3,4}[KR]H.Q.Q";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L1e signature (PS00939)
	$msg = "   Ribosomal protein L1e signature (PS00939): " if $::OPT_V;
	$pat = "AG.{2}T.AES[FYWQ]G[ST]GR";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L7Ae signature (PS01082)
	$msg = "   Ribosomal protein L7Ae signature (PS01082): " if $::OPT_V;
	$pat = "[CA].{4}VP[FY].{2}[LIVM].[GSQ][KRQ].{2}LG";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L13e signature (PS01104)
	$msg = "   Ribosomal protein L13e signature (PS01104): " if $::OPT_V;
	$pat = "[KR]Y.{2}K[LIVM]R[TA]G[KR]GF[ST]L.E";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L18e signature (PS01106)
	$msg = "   Ribosomal protein L18e signature (PS01106): " if $::OPT_V;
	$pat = "S[KR].{2}R.P[LIVM]S[LIVM]SR[LIVM]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L19e signature (PS00526)
	$msg = "   Ribosomal protein L19e signature (PS00526): " if $::OPT_V;
	$pat = "R.[KR].{5}[KR].{3}[KRH].{2}G.G.R.G.{3}AR.{3}[KQ].{2}W.{7}R.{2}L.{3}R";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L24e signature (PS01073)
	$msg = "   Ribosomal protein L24e signature (PS01073): " if $::OPT_V;
	$pat = "[FY].G.{2}I.PG.G.{2}[FY].[KRH].D";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L27e signature (PS01107)
	$msg = "   Ribosomal protein L27e signature (PS01107): " if $::OPT_V;
	$pat = "GKN.WFF.KLRF\$";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L30e signature 1 (PS00709)
	$msg = "   Ribosomal protein L30e signature 1 (PS00709): " if $::OPT_V;
	$pat = "[STA].{5}G.[QKR].{2}[LIVM][KQT].{2}[KR].G.{2}K.[LIVM]{3}";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L30e signature 2 (PS00993)
	$msg = "   Ribosomal protein L30e signature 2 (PS00993): " if $::OPT_V;
	$pat = "[DE]LG[STA].{2}G[KR].{6}[LIVM].[LIVM].[DEN].G";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L31E signature (PS01144)
	$msg = "   Ribosomal protein L31E signature (PS01144): " if $::OPT_V;
	$pat = "VR[LIVM].{3}[LIVM]N.A.W.[KR]G";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L32e signature (PS00580)
	$msg = "   Ribosomal protein L32e signature (PS00580): " if $::OPT_V;
	$pat = "F.R.{4}[KR].{2}[KR][LIVM].{3}WR[KR]P.G";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L34E signature 1 (PS01145)
	$msg = "   Ribosomal protein L34E signature 1 (PS01145): " if $::OPT_V;
	$pat = "QR[LIVM]T.[KR]{3}.{2}Y.T.SN";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L34E signature 2 (PS01146)
	$msg = "   Ribosomal protein L34E signature 2 (PS01146): " if $::OPT_V;
	$pat = "RI.RAF[LIVM]{2}[DE]{2}QK[LIVM]{2}";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L35Ae signature (PS01105)
	$msg = "   Ribosomal protein L35Ae signature (PS01105): " if $::OPT_V;
	$pat = "GK[LIVM].R.HG.{2}G.V.A.F.{3}LP";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L37e signature (PS01077)
	$msg = "   Ribosomal protein L37e signature (PS01077): " if $::OPT_V;
	$pat = "GT.S.G.[KR].{3}[ST]H.{2}C.RCG";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein L46e signature (PS00051)
	$msg = "   Ribosomal protein L46e signature (PS00051): " if $::OPT_V;
	$pat = "N.{2}RR[NH]WRR";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S2 signature 1 (PS00962)
	$msg = "   Ribosomal protein S2 signature 1 (PS00962): " if $::OPT_V;
	$pat = "[LIVMFA].{2}[LIVMFYC]{2}.[STAC][GSTANER][STALV][HY][LIVMF]G";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S2 signature 2 (PS00963)
	$msg = "   Ribosomal protein S2 signature 2 (PS00963): " if $::OPT_V;
	$pat = "P.{2}[LIVMF]{2}[LIVMS].[GDN].{3}[DENL].{3}[LIVM].E.{4}[GNQR][LIVM][AP]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S3 signature 1 (PS00548)
	$msg = "   Ribosomal protein S3 signature 1 (PS00548): " if $::OPT_V;
	$pat = "[LIVMF][RE].G.{2}[KRQA].{3}[DNSA].{2}[FYW][SANV][NQDER]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S3 signature 2 (PS00734)
	$msg = "   Ribosomal protein S3 signature 2 (PS00734): " if $::OPT_V;
	$pat = "[GQR]R[LIVM].G.E[LIVM]A[KR][LIVMSTAR]E";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S4 signature (PS00632)
	$msg = "   Ribosomal protein S4 signature (PS00632): " if $::OPT_V;
	$pat = "[LIVM][DE].RL.{3}[LIVM][VMFYH][KR].{3}[AGTCF].[ST].{3}[SAI][KR].[LIVMF]{2}";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S5 signature (PS00585)
	$msg = "   Ribosomal protein S5 signature (PS00585): " if $::OPT_V;
	$pat = "G[KRQ].{3}[FY].[ACV].{2}[LIVMA][LIVM][AG][DN].{2}G.[LIVM]G.[SAG].{5,6}[DEQ][LIVM].{2}A[LIVM]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S6 signature (PS01048)
	$msg = "   Ribosomal protein S6 signature (PS01048): " if $::OPT_V;
	$pat = "[WL]G.[KR].LAY.I";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S7 signature (PS00052)
	$msg = "   Ribosomal protein S7 signature (PS00052): " if $::OPT_V;
	$pat = "[DENSK].[LIVME].{3}[LIVMFT]{2}.{6}GK[KR].{5}[LIVMF]{2}.{2}[STA]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S8 signature (PS00053)
	$msg = "   Ribosomal protein S8 signature (PS00053): " if $::OPT_V;
	$pat = "[GE].{2}[LIV]{2}[STY]T.{2}G[LIVM]{2}.{4}A[KRHA]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S9 signature (PS00360)
	$msg = "   Ribosomal protein S9 signature (PS00360): " if $::OPT_V;
	$pat = "GGG.{2}[GSA]Q.{2}[SA].{3}[GSA].[STA][KR][GSA]L";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S10 signature (PS00361)
	$msg = "   Ribosomal protein S10 signature (PS00361): " if $::OPT_V;
	$pat = "[AV].{3}[GDNS][LIVMSTA].{3}GP[LIVM].[LIVM]PT";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S11 signature (PS00054)
	$msg = "   Ribosomal protein S11 signature (PS00054): " if $::OPT_V;
	$pat = "[DNE]VTP.[PA].[DN]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S12 signature (PS00055)
	$msg = "   Ribosomal protein S12 signature (PS00055): " if $::OPT_V;
	$pat = "[RK].PNS[AR].R";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S13 signature (PS00646)
	$msg = "   Ribosomal protein S13 signature (PS00646): " if $::OPT_V;
	$pat = "[KRQ]G[LIVMFY]RH.{2}[GSN].{2}[LIVMC]RGQ";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S14 signature (PS00527)
	$msg = "   Ribosomal protein S14 signature (PS00527): " if $::OPT_V;
	$pat = "[RP].{0,1}C.{11,12}[LIVMF].[LIVMF][SC][RG].{3}[RN]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S15 signature (PS00362)
	$msg = "   Ribosomal protein S15 signature (PS00362): " if $::OPT_V;
	$pat = "[LIVM].{2}H[LIVMFY].{5}D.{2}[SAGN].{3}[LF].{9}[LIVM].{2}Y";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S16 signature (PS00732)
	$msg = "   Ribosomal protein S16 signature (PS00732): " if $::OPT_V;
	$pat = "[LIVM].[LIVM]RL[ASK]R.G[AR]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S17 signature (PS00056)
	$msg = "   Ribosomal protein S17 signature (PS00056): " if $::OPT_V;
	$pat = "GD.[LIV].[LIV].[QEK].[RK]P[LIV]S";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S18 signature (PS00057)
	$msg = "   Ribosomal protein S18 signature (PS00057): " if $::OPT_V;
	$pat = "I[DY]Y.{2}[LIVM].{2}[LIVM].{2}[FY][LIVM][ST][DE].GK[LIVM]{2}.{2}R[LIVM]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S19 signature (PS00323)
	$msg = "   Ribosomal protein S19 signature (PS00323): " if $::OPT_V;
	$pat = "[STDN]G[KRQM].{6}[LIVM].{4}[LIVM][GS].{2}[LF][GA]EF.{2}[ST]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S4e signature (PS00528)
	$msg = "   Ribosomal protein S4e signature (PS00528): " if $::OPT_V;
	$pat = "H.KR[LIVM][SAN].P.{2}W.[LIVM].[KR]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S6e signature (PS00578)
	$msg = "   Ribosomal protein S6e signature (PS00578): " if $::OPT_V;
	$pat = "[LIVM][STAM]GG.D.{2}G.PM";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S7e signature (PS00948)
	$msg = "   Ribosomal protein S7e signature (PS00948): " if $::OPT_V;
	$pat = "RLVRELEKKFSGKH";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S17e signature (PS00712)
	$msg = "   Ribosomal protein S17e signature (PS00712): " if $::OPT_V;
	$pat = "KLQE{3}RE[KR]{2}D";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S19e signature (PS00628)
	$msg = "   Ribosomal protein S19e signature (PS00628): " if $::OPT_V;
	$pat = "[SA].{2}[LIVMA].R.[AIV]LQ.L[EQ]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S21e signature (PS00996)
	$msg = "   Ribosomal protein S21e signature (PS00996): " if $::OPT_V;
	$pat = "LYVPRKCSA";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S24e signature (PS00529)
	$msg = "   Ribosomal protein S24e signature (PS00529): " if $::OPT_V;
	$pat = "T[RKHQ]FG.{2}K[ST].G[FY][GA].[LIVM]Y[DN]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S26e signature (PS00733)
	$msg = "   Ribosomal protein S26e signature (PS00733): " if $::OPT_V;
	$pat = "YCVSCAIH";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Ribosomal protein S28e signature (PS00961)
	$msg = "   Ribosomal protein S28e signature (PS00961): " if $::OPT_V;
	$pat = "ESEREARRLR\$";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

	if($::OPT_V) {
		print $hit == 0 ? " none\n\n" : "\n";
	}
	$hit;
}


sub bactdna {
    my($pat, @pat);
    my($hit) = 0;

    if($::OPT_W) {
	print "<A HREF=\"/psort/helpwww2.html#motif\">checking 33 PROSITE prokaryotic DNA binding motifs:</A> ";
    } elsif($::OPT_V) {
	print "checking 33 PROSITE prokaryotic DNA binding motifs: ";
    }

## Prokaryotic dksA/traR C4-type zinc finger (PS01102)
	$msg = "   Prokaryotic dksA/traR C4-type zinc finger (PS01102): " if $::OPT_V;
	$pat = "C[DES].C.{3}I.{3}R.{4}P.{4}C.{2}C";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Prokaryotic transcription elongation factors signature 1 (PS00829)
	$msg = "   Prokaryotic transcription elongation factors signature 1 (PS00829): " if $::OPT_V;
	$pat = "T.{2}[GS].{3}L.{2}EL.{2}L.{3,4}R";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Prokaryotic transcription elongation factors signature 2 (PS00830)
	$msg = "   Prokaryotic transcription elongation factors signature 2 (PS00830): " if $::OPT_V;
	$pat = "S.{2}SP[LIVM][AG].[AG][LIVM][LIVMY].{4}[DG][DE]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial regulatory proteins, araC family signature (PS00041)
	$msg = "   Bacterial regulatory proteins, araC family signature (PS00041): " if $::OPT_V;
	$pat = "[KRQ][LIVMA].{2}[GSTALIV][^FYWPGDN].{2}[LIVMSA].{4,9}[LIVMF].{2}[LIVMSTA][GSTACIL].{3}[GANQRF][LIVMFY].{4,5}[LFY].{3}[FYIVA][^FYWHCM].{3}[GSADENQKR].[NSTAPKL][PARL]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial regulatory proteins, araC family DNA-binding domain profile (PS01124)
	$msg = "   Bacterial regulatory proteins, araC family DNA-binding domain profile (PS01124): " if $::OPT_V;
	$pat = "C.{2}D[LIVM].{6}[ST].{4}S[HYR][HQ]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial regulatory proteins, asnC family signature (PS00519)
	$msg = "   Bacterial regulatory proteins, asnC family signature (PS00519): " if $::OPT_V;
	$pat = "[GSTAP].{2}[DE][LIVM][SA].{2}[LIVMFY][GN][LIVMS]S.{6}RV.{2}[LIVM].{3}G";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial regulatory proteins, crp family signature (PS00042)
	$msg = "   Bacterial regulatory proteins, crp family signature (PS00042): " if $::OPT_V;
	$pat = "[LIVM][STAG][RHNW].{2}[LI][GA].[LIVMFYA][LIVSC][GA].[STAC].{2}[MST].[GSTN]R.[LIVMF].{2}[LIVMF]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial regulatory proteins, deoR family signature (PS00894)
	$msg = "   Bacterial regulatory proteins, deoR family signature (PS00894): " if $::OPT_V;
	$pat = "R.{3}[LIVM].{3}[LIVM].{17}[STA].{2}T[LIVMA]R[KRNA]D[LIVMF]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial regulatory proteins, gntR family signature (PS00043)
	$msg = "   Bacterial regulatory proteins, gntR family signature (PS00043): " if $::OPT_V;
	$pat = "[EV].{2}[LIVM].{3}[LIVMFYK].[LIVMFT][NGSTK]R.{2}[LIVM].{3}[LIVMFY].{2}[LS]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial regulatory proteins, iclR family signature (PS01051)
	$msg = "   Bacterial regulatory proteins, iclR family signature (PS01051): " if $::OPT_V;
	$pat = "[GA].{3}[DS].{2}E.{6}[CSA][LIVM][GSA].{2}[LIVM][FYH][DN]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial regulatory proteins, lacI family signature (PS00356)
	$msg = "   Bacterial regulatory proteins, lacI family signature (PS00356): " if $::OPT_V;
	$pat = "[LIVM].[DE][LIVM]A.{2}[STAGV].V[GSTP].{2}[STAG][LIVMA].{2}[LIVMFYAN][LIVMC]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial regulatory proteins, luxR family signature (PS00622)
	$msg = "   Bacterial regulatory proteins, luxR family signature (PS00622): " if $::OPT_V;
	$pat = "[GDC].{2}[NSTAV].{2}[IV][GSTA].{2}[LIVMFYWC].[LIVMFYWCR].{3}[NST][LIVM].{5}[NRHSA][LIVMTA].{2}[KR]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial regulatory proteins, lysR family signature (PS00044)
	$msg = "   Bacterial regulatory proteins, lysR family signature (PS00044): " if $::OPT_V;
	$pat = "[LIVMFYT].{2}[STGALV][STAG].{5}[PSTAV][PNQHKRIV].{2}[LIVMFA][STA].{2}[LIVMFW].{2}[LIVMFW][RKEQA].{2}[LIVMFYNTE]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial regulatory proteins, marR family signature (PS01117)
	$msg = "   Bacterial regulatory proteins, marR family signature (PS01117): " if $::OPT_V;
	$pat = "[STNA][LIA].[RNGS].{4}[LM][EIV].{2}[GES][LFYW][LIVC].{7}[DN][RKQG][RK].{6}T.{2}[GA]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial regulatory proteins, merR family signature (PS00552)
	$msg = "   Bacterial regulatory proteins, merR family signature (PS00552): " if $::OPT_V;
	$pat = "[GSA].[LIVMFA][ASM].{2}[STACLIV][GSDENQR][LIVC][STANHK].{3}[LIVM][RHF].[YW][DEQ].{2,3}[GHDNQ][LIVMF]{2}";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial regulatory proteins, tetR family signature (PS01081)
	$msg = "   Bacterial regulatory proteins, tetR family signature (PS01081): " if $::OPT_V;
	$pat = "G[LIVMFYS].{2,3}[TS][LIVMT].{2}[LIVM].{5}[LIVQS][STAGENQH].[GPAR].[LIVMF][FYST].[HFY][FV].[DNST]K.{2}[LIVM]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Transcriptional antiterminators bglG family signature (PS00654)
	$msg = "   Transcriptional antiterminators bglG family signature (PS00654): " if $::OPT_V;
	$pat = "[ST].H.{2}[FA]{2}[LIVM][EQK]R.{2}[QNK]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Sigma-54 factors family signature 1 (PS00717)
	$msg = "   Sigma-54 factors family signature 1 (PS00717): " if $::OPT_V;
	$pat = "P[LIVM].[LIVM].{2}[LIVM]A.{2}[LIVM].[LIVM]H.ST[LIVM]SR";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Sigma-54 factors family signature 2 (PS00718)
	$msg = "   Sigma-54 factors family signature 2 (PS00718): " if $::OPT_V;
	$pat = "RRTV[AT]KYR";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Sigma-70 factors family signature 1 (PS00715)
	$msg = "   Sigma-70 factors family signature 1 (PS00715): " if $::OPT_V;
	$pat = "D[LIVMF]{2}[HEQS].G.[LIVMFA]GL[LIVMFYE].[GSAM][LIVMAP]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Sigma-70 factors family signature 2 (PS00716)
	$msg = "   Sigma-70 factors family signature 2 (PS00716): " if $::OPT_V;
	$pat = "[STN].{2}[DEQ][LIVM][GAS].{4}[LIVMF][STG].{3}[LIVMA].[NQR][LIVMA][EQH].{3}[LIVM].{2}[LIVM]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Sigma-70 factors ECF subfamily signature (PS01063)
	$msg = "   Sigma-70 factors ECF subfamily signature (PS01063): " if $::OPT_V;
	$pat = "[STAIV].D[LIVM][LIVMA]Q.[STAV][LIVMFY][LIVMA].[STAV][LIVMFYW]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Sigma-54 interaction domain ATP-binding region A signature (PS00675)
	$msg = "   Sigma-54 interaction domain ATP-binding region A signature (PS00675): " if $::OPT_V;
	$pat = "[LIVMFY]{3}.G[DEQ][STE]G[STAV]GK.{2}[LIVMFY]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Sigma-54 interaction domain ATP-binding region B signature (PS00676)
	$msg = "   Sigma-54 interaction domain ATP-binding region B signature (PS00676): " if $::OPT_V;
	$pat = "[GS].[LIVMF].{2}A[DNEQASH][GNEK]G[STIM][LIVMFY]{3}[DE][EK][LIVM]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Sigma-54 interaction domain C-terminal part signature (PS00688)
	$msg = "   Sigma-54 interaction domain C-terminal part signature (PS00688): " if $::OPT_V;
	$pat = "[FYW]P[GS]N[LIVM]R[EQ]L.[NHAT]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Bacterial histone-like DNA-binding proteins signature (PS00045)
	$msg = "   Bacterial histone-like DNA-binding proteins signature (PS00045): " if $::OPT_V;
	$pat = "[GSK]F.{2}[LIVMF].{4}[RKEQA].{2}[RST].[GA].[KN]P.T";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Dps protein family signature 1 (PS00818)
	$msg = "   Dps protein family signature 1 (PS00818): " if $::OPT_V;
	$pat = "HW.[LIVM].G.{5}[LIVM]H.{3}[DE]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Dps protein family signature 2 (PS00819)
	$msg = "   Dps protein family signature 2 (PS00819): " if $::OPT_V;
	$pat = "[DE].[LIVM]AER.{3}[LI]G.{2}[PA]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## DnaA protein signature (PS01008)
	$msg = "   DnaA protein signature (PS01008): " if $::OPT_V;
	$pat = "I[GA].{2}[LIVMF][SGDNK].{0,1}[KR].H[STP][STV][LIVM]{2}.[SA].{2}[KRE][LIVM]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## RecF protein signature 1 (PS00617)
	$msg = "   RecF protein signature 1 (PS00617): " if $::OPT_V;
	$pat = "P[ED].{3}[LIVM]{2}.G[GSAD]P.{2}RR.[FY][LIVM]D";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## RecF protein signature 2 (PS00618)
	$msg = "   RecF protein signature 2 (PS00618): " if $::OPT_V;
	$pat = "[LIVMFY]{2}.D.{2,3}[SA]ELD.{2}[KRH].{3}L";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Small, acid-soluble spore proteins, alpha/beta type, signature 1 (PS00304)
	$msg = "   Small, acid-soluble spore proteins, alpha/beta type, signature 1 (PS00304): " if $::OPT_V;
	$pat = "K.E[LIV]A.[DE][LIVMF]G[LIVMF]";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

## Small, acid-soluble spore proteins, alpha/beta type, signature 2 (PS00684)
	$msg = "   Small, acid-soluble spore proteins, alpha/beta type, signature 2 (PS00684): " if $::OPT_V;
	$pat = "[KR][SAQ].G.VGG.[LIVM].[KR]{2}[LIVM]{2}";
	@pat = $::seq =~ /$pat/go;
	$hit += &report_motif2($msg, @pat);

	if($::OPT_V) {
		print $hit == 0 ? " none\n\n" : "\n";
	}
	$hit;
}

1;
