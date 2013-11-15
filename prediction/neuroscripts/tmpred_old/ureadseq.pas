(*	Ureadseq.p
 *
 *	Copyright 1990 by d.g.gilbert
 *	dogStar Software && Indiana University Biology Dept.
 *	email: gilbertd@iubio.bio.indiana.edu
 *
 *	read nucleic/protein sequence in these formats:
 *   Stanford/IG, Genbank, NBRF, EMBL, UWGCG,
 *   DNAStrider, Fitch, Pearson, NCR/Zuker
 *	data files may have multiple sequences
 *
 * This program may be freely copied and used by anyone.
 * Developers are encourged to incorporate parts in their
 * programs, rather than devise their own private sequence
 * format.
 *
 * Unix and VAX-VMS pascal compilers barf on this for various
 * reasons, but any UCSD-flavor pascal should do okay.
 *)

 (*-----------------------------------------------------------------------*)
 (* modifications for use with TURBO-PASCAL 5.0 /5.5 introduced by
    K.Hofmann  (KHOFMANN@cipvax.biolan.uni-koeln.de) are:

    o Turbo-Pascal insists on ELSE instead of OTHERWISE as an excluding
      possibility of a CASE statement.

    o Turbo Pascal does not like procedure-variables as parameters of
      other functions if these procedures are defined locally as they
      were in the original MPW-pascal-version. So i moved the endXXX-
      functions out of the readSeq block and forced FAR-coding by
      including the appropriate compliler directives ($F+/$F-).

    o Turbo Pascal does the opening of files in two steps whereas most
      standard compilers do this in one step. An example:
      TurboPascal: assign(file,filename);
                   reset(file);

      misc.Pascal: reset(file,filename);

      VAX-Pascal : open(file,filename,readonly);
                   reset(file);
                                                                         *)
(*-----------------------------------------------------------------------*)

{$ifdef osf1 }
  [ENVIRONMENT ('ureadseq.pen'),
   INHERIT ('dep_dec.pen')]
  module ureadseq(output);
{$endif }

{$ifdef vms }
  [ENVIRONMENT ('ureadseq.pen'),
   INHERIT ('dep_vms.pen')]
  module ureadseq;
{$endif }

{$ifdef Dos}
  unit ureadseq;
  INTERFACE
  uses dep_tp;
{$endif }

{$ifdef Windows}
  unit ureadseq;
  INTERFACE
  uses dep_win;
{$endif }


CONST
	NEWLINE	= 13; {end of line char}

	{errors}
	eFileNotFound  	=	-1;
	eNoData 				=	-2;
	eMemFull 				=	-3;
	eItemNotFound  	=	-4;
	eWrongFormat		=	-5;
	
	{formats}
	kIG				= 1;
	kGenBank 	= 2;
	kNBRF			= 3;
	kEMBL			= 4;
	kGCG			= 5;
	kStrider 	= 6;
	kFitch		= 7;
	kPearson 	= 8;
	kZuker 		= 9;
	kUnknown 	= 10;
	kMinFormat = 1;
	kMaxFormat = 10;

        endIG      = 1;
	endStrider = 2;
	endGB      = 3;
	endNBRF    = 4;
	endPearson = 5;
	endEMBL    = 6;
	endZuker   = 7;
	endFitch   = 8;

TYPE
	Sequence     = PACKED ARRAY [0..32000] OF CHAR;
	SeqPtr       = ^Sequence;

{$ifdef Borland}

	FUNCTION readSeq(
			    choice: integer;     {0=list seqs in file, >0 = read item}
			    fname : String255;   {file name}
			    seq   : SeqPtr;      {sequence or list storage}
			VAR seqlen: longint; 	{in: max seq size, out: seq size}
			VAR nseq  : integer;    {out (list): # seqs in file}
			VAR format: integer;    {out (list), in (>0): seq format}
			VAR seqId : String255   {out (seq): seq info}
                                   ): integer;  {return error: 0=okay}

	PROCEDURE writeSeq(
		       VAR out    : text;
                           seq    : seqPtr;
                           seqlen : longint;
			   inform,              {format of input seq}
			   outform: integer;    {format for output seq}
                         VAR seqid: string255);
{$endif}

{$ifdef Borland}
  IMPLEMENTATION
{$endif}

CONST 
	NNB				= 18;
	aminos				= 'ABCDEFGHIKLMNPQRSTVWXYZ*';
	iubbase				= 'ACGTUMRWSYKVHDBXN.'; 
	igbase				= 'ACGTUJRLMYKNNNNNN?';
	primenuc			= 'ACGTU';
	protonly			= 'EFIPQZ';

TYPE
	basestr	= string20;

{  endfunc = function( VAR s: String; VAR addend: boolean): boolean;}


FUNCTION islower(c:char): boolean;
BEGIN 	islower:= (c >= 'a') AND (c <= 'z');
END;

FUNCTION min(a,b:longint): longint;
BEGIN  IF (a<b) THEN min:= a ELSE min:= b;
END;

{----
FUNCTION _toupper( c: char):char;
BEGIN  _toupper:= chr(ord(c)-32);
END;
FUNCTION _tolower( c: char): char;
BEGIN  _tolower:= chr(ord(c)+32);
END;
------}

FUNCTION toupper(c: char):char;
BEGIN
	IF (c >= 'a') AND (c <= 'z') THEN toupper:= chr(ord(c)-32)
	ELSE toupper:= c;
END;

FUNCTION f_endIG(s: string255; var addend:boolean): boolean;
BEGIN
	addend:= true;
	f_endIG := (indx(s,'1') > 0) OR (indx(s,'2') > 0);
END;

FUNCTION f_endStrider(s: string255; var addend:boolean): boolean;
BEGIN
	addend:= false;
	f_endStrider := (indx(s,'//') > 0);
END;

FUNCTION f_endGB(s: string255; var addend:boolean):boolean;
BEGIN
	addend:= false;
	f_endGB := (indx(s,'//') > 0) OR (indx(s,'LOCUS') = 1);
END;

FUNCTION f_endNBRF(s:string255; VAR addend:boolean):boolean;
var  k: integer;
BEGIN
	k:= indx(s,'*');
	if k > 0 then begin
		s[0]:= chr(k-1);
		addend:= true;
		f_endNBRF:= true;
		end
	else if (indx(s,'>') = 1) then begin
		addend:= false;
		f_endNBRF:= true;
		end
	else
		f_endNBRF:= false;
END;

FUNCTION f_endPearson(s:string255; VAR addend:boolean):boolean;
BEGIN
	addend:= false;
	f_endPearson := (indx(s,'>') = 1);
END;

FUNCTION f_endEMBL(s:string255; VAR addend:boolean):boolean;
BEGIN
	addend:= false;
	f_endEMBL := (indx(s,'ID   ') = 1);
END;

FUNCTION f_endZuker(s:string255; VAR addend:boolean):boolean;
BEGIN
	addend:= false;
	f_endZuker := (indx(s,'(') = 1);
END;

FUNCTION f_endFitch(s:string255; VAR addend:boolean):boolean;
BEGIN
	addend:= false;
	f_endFitch := (s[1] <> ' ');
END;


FUNCTION readSeq(
		   choice: integer;     {0=list seqs in file, >0 = read item}
                   fname : string255;   {file name}
                   seq   : SeqPtr;      {sequence or list storage}
               VAR seqlen: longint;     {in: max seq size, out: seq size}
               VAR nseq  : integer;     {out (list): # seqs in file}
               VAR format: integer;     {out (list), in (>0): seq format}
               VAR seqId : String255    {out (seq): seq info}
                        ): integer;     {return error: 0=okay}

VAR
	allDone, done, gotuw, isfitch, addit, doadd	: boolean;
	err, i, l  : integer;
	maxseq, ninfo	: longint;
	f	: text;
	s	: string255;


	FUNCTION isSeqChar(c: CHAR): boolean;
		{---
		Various programs allow all upper and lower case
		letters, periods (.), and various symbols. Change to
		suit your needs. Digits ['0'..'9'] are NOT allowed here.
		----}
	BEGIN
		IF (c <= ' ') THEN isSeqChar := false
		ELSE IF ((c >= 'A') AND (c <= 'Z'))
				 OR ((c >= 'a') AND (c <= 'z'))
				 OR (c = '_') OR (c = '?') OR (c = '*')
				 OR (c = '.') OR (c = '-')
		THEN isSeqChar := true
		ELSE isSeqChar := false;
	END;     {isSeqChar}

	PROCEDURE addseq( s: string255);
	VAR i : integer;
	BEGIN
		FOR i := 1 TO length(s) DO BEGIN
			IF isSeqChar(s[i]) AND addit THEN BEGIN
				if seqlen >= maxseq then err:= eMemFull
				else begin
					seq^[seqlen] := s[i];
					seqlen := seqlen + 1;
					END;
				END;
			END;
	END; {addseq}

	PROCEDURE itoc( n: integer; var s: string255; var i: integer);
	BEGIN
		IF (n < 0) THEN BEGIN
			i:= i+1; s[i]:= '-';
			itoc( -n, s, i);
			END
		ELSE BEGIN
			{! C hates recursion}
			IF (n >= 10) THEN itoc( n DIV 10, s, i);
			i:= i+1;
			s[i]:= chr( ord('0') + n MOD 10);
			END;
	END; { itoc }

	PROCEDURE addinfo(VAR s: string255);
	VAR  i, l: integer; snum: string255;
	BEGIN
                i:= 0; 
	     {  itoc( nseq, snum, i); 
		snum[0]:= chr(i);          this looks criminal}
		snum:=int2str(nseq,5);
		s:= snum+')  '+s;
	        l:= length(s) + 1;
             {  s[l] := chr(NEWLINE);      and this looks suspect}
		s:=s+chr(NEWLINE);
		FOR i:=1 TO l DO IF (ninfo < maxseq) THEN BEGIN
                  seq^[ninfo] := s[i];
                  ninfo := ninfo + 1;
                END;
	END;  {addinfo}

	PROCEDURE readLoop(margin   : integer; 
	                   addfirst : boolean; 
			   endtest  : integer);
	VAR
			addend : boolean;
	BEGIN
		nseq:= nseq+1;
		seqlen:= 0;
		if choice = 0 then doadd:= false
		else doadd:= nseq = choice;
		IF addfirst THEN addseq(s);  {! fitch 1st string}
		REPEAT
                  {check eof Before read !}
                  done := eof(f);
                  if not done then readln(f, s);
                  {done := done OR endTest(s, addend);}
                  if endTest=endIG then done:=done or f_endIG(s,addend)
                  else
                  if endTest=endStrider then done:=done or f_endStrider(s,addend)
                  else
                  if endTest=endGB then done:=done or f_endGB(s,addend)
                  else
                  if endTest=endNBRF then done:=done or f_endNBRF(s,addend)
                  else
                  if endTest=endPearson then done:=done or f_endPearson(s,addend)
                  else
                  if endTest=endEMBL then done:=done or f_endEMBL(s,addend)
                  else
                  if endTest=endZuker then done:=done or f_endZuker(s,addend)
                  else
                  if endTest=endFitch then done:=done or f_endFitch(s,addend);
		  IF doadd AND (addend OR NOT done) THEN BEGIN
                    IF margin > 0 THEN s:=substr(s,margin+1,length(s)-margin);
                    addseq(s);
                  END;
		UNTIL done;
		if choice = 0 then addinfo(seqID)
		else allDone:= nseq >= choice;
	END;     {readLoop}

	PROCEDURE readIG;         {IG -- many seqs/file }
	BEGIN
		format:= kIG;
		WHILE NOT allDone DO BEGIN
			REPEAT
				readln(f, s);
			UNTIL eof(f) OR ((s <> '') AND (s[1]<>';'));
			IF eof(f) THEN allDone:= TRUE
			ELSE BEGIN
				seqId := s;
				readLoop(0, false, endIG);
				END;
			END;
	END;     {readIG}




	PROCEDURE readStrider;
	BEGIN
		format:= kStrider;
		WHILE NOT allDone DO BEGIN
			{Strider -- 1 seq only? }
			readln(f, s);
			seqId := s;
			seqID:=substr(seqId,2,length(seqID)-1); 
                        REPEAT
				readln(f, s);
			UNTIL eof(f) OR (indx(s,';') <> 1);
			readLoop(0, true, endStrider);
			IF NOT allDone THEN
			 WHILE NOT (eof(f) OR ((s <> '') AND (indx(s,';') <> 1))) DO
				readln(f, s);
			IF eof(f) THEN allDone:= TRUE;
			END;
	END;     {readStrider}



	PROCEDURE readGenBank;
	BEGIN
		format:= kGenBank;
		WHILE NOT allDone DO BEGIN
			{GenBank -- many seqs/file }
			seqId := s;
			seqID:=substr(seqId,13,length(seqID)-12); 
			REPEAT
				readln(f, s);
			UNTIL eof(f) OR ((s <> '') AND (indx(s,'ORIGIN') = 1));
			readLoop(0, false, endGB);
			IF NOT allDone THEN
			 WHILE NOT (eof(f) OR ((s <> '') AND (indx(s,'LOCUS') = 1))) DO
				readln(f, s);
			IF eof(f) THEN allDone:= TRUE;
			END;
	END;     {readGenBank}
	


	PROCEDURE readNBRF;
	BEGIN
		format:= kNBRF;
		WHILE NOT allDone DO BEGIN
			{NBRF -- many seqs/file }
			seqId := s;
			seqID:=substr(seqId,5,length(seqID)-4); 
			readln(f, s);       {junk line}
			readLoop(0, false, endNBRF);
			IF NOT allDone THEN
			 WHILE NOT (eof(f) OR ((s <> '') AND (indx(s,'>DL;') = 1))) DO
				readln(f, s);
			IF eof(f) THEN allDone:= TRUE;
			END;
	END;     {readNBRF}




	PROCEDURE readPearson;
	BEGIN
		format:= kPearson;
                WHILE NOT allDone DO BEGIN
			{Pearson -- many seqs/file }
			seqId := s;
			seqID:=substr(seqId,2,length(seqID)-1); 
			readLoop(0, false, endPearson);
			IF NOT allDone THEN
			 WHILE NOT (eof(f) OR ((s <> '') AND (indx(s,'>') = 1))) DO
				readln(f, s);
			IF eof(f) THEN allDone:= TRUE;
			END;
	END;     {readPearson}




	PROCEDURE readEMBL;
	BEGIN
		format:= kEMBL;
		WHILE NOT allDone DO BEGIN
			seqId := s;
			seqID:=substr(seqId,6,length(seqID)-5); 
			REPEAT
				readln(f, s);
			UNTIL eof(f) OR (indx(s,'SQ   ') = 1);
			readLoop(0, false, endEMBL);
			IF NOT allDone THEN
			 WHILE NOT (eof(f) OR ((s <> '') AND (indx(s,'ID   ') = 1))) DO
				readln(f, s);
			IF eof(f) THEN allDone:= TRUE;
		END;
	END;     {readEMBL}




	PROCEDURE readZuker;
	BEGIN
		format:= kZuker;
		WHILE NOT allDone DO BEGIN
			{Zuker -- many seq/file ?}
			{! 1st string is Zuker's Fortran format }
			readln(f, s);       {s == "seqLen seqId string..."}
			seqId := s;
			seqID:=substr(seqId,7,length(seqID)-6); 
			readLoop(0, false, endZuker);
			IF NOT allDone THEN
			 WHILE NOT (eof(f) OR ((s <> '') AND (indx(s,'(') = 1))) DO
				readln(f, s);
			IF eof(f) THEN allDone:= TRUE;
		END;
	END;     {readZuker }



	PROCEDURE readFitch;
	VAR first : boolean;
	BEGIN
		format:= kFitch;
		first := true;
		WHILE NOT allDone DO BEGIN
			{Fitch -- many seqs/file }
			IF NOT first THEN seqId := s;
			readLoop(0, first, endFitch);
			IF eof(f) THEN allDone:= TRUE;
			first := false;
			END;
	END;     {readFitch }
	
	
	PROCEDURE readUnknown;
	BEGIN
		format:= kUnknown; {Unknown -- 1seq/file}
		addit := choice > 0;
		nseq:= nseq+1;
		seqlen:= 0;
		addseq( seqId);        {from above..}
		seqId := fname+'  [Unknown form]';
		REPEAT
			addseq(s);
			done := eof(f);
			readln(f, s);
		UNTIL done;
		if choice = 0 then addinfo(seqId);
		allDone:= true;
	END;     {readUnknown }
	

	PROCEDURE readUWGCG;
	BEGIN
		format:= kGCG;
		{ UWGCG -- 1seq/file}
		addit := choice > 0;
		nseq:= nseq+1;
		seqlen:= 0;
		seqId := s;
		REPEAT
			done := eof(f);
			if not eof(f) then readln(f, s);
			IF (NOT done) and (length(s)>0) THEN BEGIN
				{delete(s, 1, 9);}  {skip margin}
                                {s:=substr(s,10,length(s)-9);}
				addseq(s);
				END;
		UNTIL done;
		if choice = 0 then addinfo(seqID);
		allDone:= true;
	END;     {readUWGCG }
	
BEGIN
	maxseq:= seqlen;
	seqlen:= 0;
	ninfo := 0;
	nseq  := 0;
	seqId := '';
	readSeq:= 0;
	if choice < 0 then choice:= 0;
	if (choice = 0)
	 or (format < kMinFormat) or (format > kMaxFormat)
		then format:= 0;
	addit := choice > 0;
	allDone:= false;

	depend_err:='';
	fopen(f,fname,'r');
	if depend_err='' then err:=0 else err:=1;
	
	IF err <> 0 THEN begin
		err:= eFileNotFound;
		allDone:= TRUE;
		end
	
	ELSE IF (format > 0) THEN BEGIN
		{don't need to re-check format}
		REPEAT 
			readln(f, s);
			l:= length(s);
			while (l>0) and (s[l]=' ') do l:= l-1;
		UNTIL (l<>0) OR eof(f);
		IF (l = 0) THEN err:= eNoData
		ELSE CASE format OF
			kUnknown: readUnknown;
			kIG     : readIG;
			kStrider: readStrider;
			kGenBank: readGenBank;
			kNBRF   : readNBRF;
			kPearson: readPearson;
			kEMBL   : readEMBL;
			kZuker	: readZuker;
			kFitch	: BEGIN
						seqId := s;
						readln(f, s);
						readFitch;
						END;
			kGCG		: REPEAT
						gotuw := indx(s,'..') > 0;
						IF gotuw THEN readUWGCG;
						readln(f, s);
					 UNTIL eof(f) or allDone;
			END;
		END
		
	ELSE BEGIN {find format}
		{ check for ".." of uwgcg, since it can masquerade as any
			other format }
		i := 0;
		REPEAT
			readln(f, s);
			i := i + 1;
			gotuw := indx(s,'..') > 0;
			IF gotuw THEN gotuw := indx(s,'Check:') > 0;
			done := gotuw OR eof(f) OR (i > 800);
			{! ECOLAC UW/GenBank document header is 300 lines !}
			IF (i < 5) AND (indx(s,';') = 1) THEN BEGIN
				gotuw := false; 
				done := true; {fix for ToIG/ToNBRF/ToEMBL of UWGCG ?}
				END;
		UNTIL done;

		IF gotuw THEN readUWGCG
		ELSE BEGIN
			reset(f);
			REPEAT 
				readln(f, s);
				l:= length(s);
				while (l>0) and (s[l]=' ') do l:= l-1;
			UNTIL (l<>0) OR eof(f);
			IF (l=0) THEN err:= eNoData;
			END;

		IF (err<>0) OR allDone THEN {skip to end}
				
		ELSE IF indx(s,';') = 1 THEN BEGIN
			IF indx(s,'Strider') > 0 THEN readStrider
			ELSE readIG;
			END

		ELSE IF indx(s,'LOCUS') = 1 THEN readGenBank

		ELSE IF indx(s,'>DL;') = 1 THEN readNBRF

		ELSE IF indx(s,'>') = 1 THEN readPearson

		ELSE IF indx(s,'ID   ') = 1 THEN readEMBL

		ELSE IF indx(s,'(') = 1 THEN readZuker

		ELSE BEGIN
			seqId := s;
			readln(f, s);         {test for fitch format}
			i := 1;
			REPEAT
				isfitch := (((i - 1) MOD 4  = 0) AND (s[i] = ' ')) 
						    OR (((i - 1) MOD 4 <> 0) AND (s[i] <> ' '));
				i := i + 1;
			UNTIL (i >= length(s)) OR NOT isfitch;
			IF isfitch THEN 
				readFitch
			ELSE 
				readUnknown;
			END;
		END;

	IF choice = 0 THEN seqlen := ninfo;
	close(f);
	readSeq := err;
END; {readSeq}



function GCGchecksum( seq: seqPtr; seqlen: longint): integer;
VAR
	check, count, i: longint;
begin
	check:= 0; count:= 0;
	for i:= 0 to seqlen-1 do begin
		count:= count+1;
		check:= check + count * ord(toupper(seq^[i]));
		if (count = 57) then count:= 0;
		end;
	GCGchecksum:= check MOD 10000;
END; {GCGchecksum}	


PROCEDURE mapbase( seq: seqPtr; seqlen: longint;
									 frombase, tobase: basestr);
{ translate ambiguity codes b/n ig and iub/gcg }
LABEL 1,2;
VAR
  i         : longint;
  c         : char;
  isnuc, 
  islow     : boolean;
  pn,pp,
  j,maxtest : integer;
  s1        : char;
BEGIN
	s1:= '.';
	maxtest := min(200, seqlen);
	pn:= 0; pp:= 0; isnuc:= true;
	FOR i:= 0 TO maxtest-1 DO BEGIN
		c := seq^[i];
		IF islower(c) THEN c:= chr(ord(c)-32); {_toupper}
		s1:=c;
		IF indx(protonly,s1) > 0 THEN BEGIN
			isnuc:= false; GOTO 1; {break;}
			END
		ELSE IF indx(primenuc,s1) > 0 THEN pn:= pn+1
		ELSE IF indx(aminos,s1) > 0 THEN pp:= pp+1;
		END;
1:
	IF (isnuc AND (pn > pp)) THEN 
		FOR i:= 0 TO seqlen-1 DO BEGIN
			c:= seq^[i];
			islow:= islower(c);
			IF islow THEN c:= chr(ord(c)-32); {_toupper}
			FOR j:= 1 TO NNB DO 
				IF (c = frombase[j]) THEN BEGIN
					c:= tobase[j];
					IF islow THEN seq^[i]:= chr(ord(c)+32) {_tolower}
					ELSE seq^[i]:= c;
					GOTO 2;{break;}
					END;
2:
			END;
END; {mapbase}


PROCEDURE writeSeq(
                   VAR out     : text;
                       seq     : seqPtr; 
                       seqlen  : longint;
                       inform,                   {format of input seq}
                       outform : integer;        {format for output seq}
                   VAR seqid   : string255);
		   
{ dump sequence to standard output }
VAR 
	spacer, 
	width, 
	tab     : integer;
	numline, 
	endline : boolean;
	idword  : string40;
	endstr  : string10;
	s       : string100;
	i,j,l, 
	l1,line : longint;

BEGIN	 
	spacer:= 0; width:= 50; tab:= 0;
	numline:= false;
	
	i:= 0;
	while (i<length(seqid)) and (seqid[i] = ' ') do i:= i+1; 
{	if i>1 then delete(seqid, 1, i-1);}
	if i>1 then seqid:=substr(seqid,i,length(seqID)-i+1);
	{- sscanf( seqid, "%30s", idword); }
	i:= 0;
	while (i<length(seqid)) and (seqid[i+1] <> ' ') do i:= i+1;
	idword:= substr(seqid, 1, min(30, i));
	endstr:= '';
	endline:= false;
	
	if outform = kUnknown then
	  begin  	{ no header, just sequence }
            if (inform = kIG) then mapbase(seq, seqlen, igbase, iubbase);
            endline:= true;
          end
	else if outform= kGenBank then 
	  begin
            if (inform = kIG) then mapbase(seq, seqlen, igbase, iubbase);
            writeln(out, 'LOCUS       ', idword);
            writeln(out, 'ORIGIN      ', seqid,', ',seqlen:1,' bases.');
            spacer:= 11;
            numline:= true;
            endstr:= '//'; 
            endline:= true;
          END
	else if outform= kNBRF then 
	  begin
            if (inform = kIG) then mapbase(seq, seqlen, igbase, iubbase);
            writeln(out,'>DL;', idword);
            writeln(out,seqid, ', ', seqlen:1,' bases.');
            writeln(out);
            spacer:= 11;
            numline:= true;
            endstr:= '*';
          END
	else if outform=kEMBL then
	  begin
            if (inform = kIG) then mapbase(seq, seqlen, igbase, iubbase);
            writeln(out,'ID   ', idword);
            writeln(out,'SQ   Sequence ',seqid,', ',seqlen,' bases.');
          END
	else if outform=kGCG then
	  begin
            if (inform = kIG) then mapbase(seq, seqlen, igbase, iubbase);
            writeln(out,seqid);
            writeln(out,'    ',idword,'  Length: ', seqlen:1,
	      '  (today)  Check: ', GCGchecksum(seq, seqlen):1,'  ..');
            spacer:= 11;
            numline:= true;
          END
	else if outform=kStrider then
	  begin
            if (inform = kIG) then mapbase(seq, seqlen, igbase, iubbase);
            writeln(out,'; ### from DNA Strider ;-)');
            writeln(out,'; DNA sequence  ',seqid,', ',seqlen:1,' b.p.');
            writeln(out,';');
            endstr:= '//'; 
            endline:= true;
          END
	else if outform=kFitch then
	  begin
            if (inform = kIG) then mapbase(seq, seqlen, igbase, iubbase);
            writeln(out,seqid, ', ',seqlen:1,' bases.');
            spacer:= 4;
            width:= 60;
          END
	else if outform=kIG then 
	  begin
            if (inform <> kIG) then mapbase(seq, seqlen, iubbase, igbase);
            writeln(out,';',seqid, ', ',seqlen:1,' bases.');
            writeln(out,idword);
            endstr:= '1'; { == linear dna }
          END
  	else 
	  begin
            if (inform = kIG) then mapbase(seq, seqlen, igbase, iubbase);
            writeln(out,'>',seqid, ', ',seqlen:1,' bases.');
          END;

	width:= min(width,100);
	l:= 0; l1:= 0; line:= 1;
	for i:= 0 to seqlen-1 do begin
		if (spacer > 0) THEN IF ((l+1) mod spacer = 1) then
			begin l:= l+1; s[l]:= ' ';  
			end;
		l:= l+1; s[l]:= seq^[i]; l1:= l1+1;
		if (l1 = width) or (i = seqlen-1) then begin
			if numline then write(out, (i+1):8, ' ');
			if tab > 0 then write(out,' ':tab);
			s[0]:= chr(l); l:= 0; l1:= 0;
			if (i = seqlen-1) AND NOT endline then 
				writeln(out,s, endstr)
			else 
				writeln(out,s);
			line:= line+1;
			end;
		end;
	if endline then writeln(out, endstr);
end; {writeseq}


{$ifdef Borland}
  begin
{$endif}
      
   end.

