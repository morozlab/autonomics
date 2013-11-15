{$ifdef osf1 }
  [ENVIRONMENT ('seq.pen'),
   INHERIT ('dep_dec.pen','ureadseq.pen')]
  module seq(output);
{$endif }

{$ifdef Dos}
  unit seq;
  INTERFACE
  uses dep_tp,ureadseq;
{$endif }

{$ifdef Windows}
  unit seq;
  INTERFACE
  uses dep_win,ureadseq;
{$endif }


const
      prot               = true;
      dna                = false;
var
      sq_prot_code       : string20;
      sq_legal           : set of char;
      sq_protmode        : boolean;
      read_err : string255;

{$ifdef Borland}

procedure sq_read_file(
                       var buffer;
                       var seqlen   : longint;     {has to be initialized}
                       name         : string255);  {with maximum seqlen}

function sq_read_multifile(
                           var buffer;
                           var seqlen       : longint;
                           var seqid        : string255;
                               name         : string255;
                               choice       : integer;
                           var maxnseq      : integer ): integer;    {returns 0 if successful}

procedure sq_write_file
                       (var buffer;
                            len    : integer;
                            name   : string255);

procedure sq_set_legal(
                       isprot : boolean);

{$endif}

{----------------------------------------------------------------------------}

{$ifdef Borland}
  IMPLEMENTATION
{$endif}


const
      maxbuf                = 32000;
type
      bt            = packed array[1..maxbuf] of char;


procedure sq_read_file(
{$ifdef Borland}                           
                       var buffer;
{$endif}                           
{$ifdef dec}                           
                       var buffer   : packed array[b0..b1:integer] of char;
{$endif}                           
                       var seqlen   : longint;
                       name         : string255);
var
   choice : integer;
   i      : integer;  
   seq    : SeqPtr;   {defined in ureadseq}
   nseq   : integer;
   format : integer;
   seqid  : string255;
   seqerr : integer;

begin
  new(seq);
  choice:=1;                         {read only first sequence}
  format:=0;                         {read any format         }
  seqerr:=readseq(choice,name,seq,seqlen,nseq,format,seqid);
  if seqerr=0 then begin
{$ifdef Borland}                           
    for i:=1 to seqlen do bt(buffer)[i]:=seq^[i-1];
{$endif}                           
{$ifdef dec}                           
    for i:=1 to seqlen do buffer[i]:=seq^[i-1];
{$endif}                           
  end
  else seqlen:=0;
  dispose(seq);
end;

function sq_read_multifile(

{$ifdef Borland}                           
			   var buffer;
{$endif}                           
{$ifdef dec}                           
			   var buffer       : packed array[b0..b1: integer] of char;
{$endif}                           
			   var seqlen       : longint;
                           var seqid        : string255;
                               name         : string255;
                               choice       : integer;
                           var maxnseq      : integer ): integer;    {returns 0 if successful}
var
   i      : integer;
   seq    : SeqPtr;
   format : integer;
   seqerr : integer;
begin
  new(seq);
  format:=0;                         {read any format         }
  seqerr:=readseq(choice,name,seq,seqlen,maxnseq,format,seqid);
  if seqerr=0 then begin
{$ifdef Borland}                           
    for i:=1 to seqlen do bt(buffer)[i]:=seq^[i-1];
{$endif}                           
{$ifdef dec}                           
    for i:=1 to seqlen do buffer[i]:=seq^[i-1];
{$endif}                           
  end
  else seqlen:=0;
  dispose(seq);
  sq_read_multifile := seqerr;
  if choice>maxnseq then sq_read_multifile:=eItemNotFound;
end;

procedure sq_write_file(
{$ifdef Borland}                           
		        var buffer;
{$endif}                           
{$ifdef dec}                           
			var buffer  : packed array[b0..b1: integer] of char;
{$endif}                           
                            len     : integer;
                            name    : string255);
var sfile : text;
    i : integer;
begin
   fopen(sfile,name,'w');
   for i:=1 to len do begin
     if i mod 60 =1 then write(sfile,'     ')
     else if i mod 10 =1 then write(sfile,' ');
{$ifdef Borland}                           
     write(sfile,bt(buffer)[i]);
{$endif}                           
{$ifdef dec}                           
     write(sfile,buffer[i]);
{$endif}                           
     if i mod 60 =0 then writeln(sfile);
   end;
   close(sfile);
end;

procedure sq_set_legal(
                       isprot : boolean);
begin
  if isprot then begin
                sq_legal:=['A'..'Z','-','*'];
                sq_protmode:=true;
              end
              else begin
                sq_legal:=['A','C','G','T','U','N','*','-'];
                sq_protmode:=false;
              end;
end;




{$ifdef Borland}
  begin
    sq_prot_code:= 'ARNDCQEGHILKMFPSTWYV';
    
{$endif    }

{$ifdef dec}
  to begin do begin
    sq_prot_code:= 'ARNDCQEGHILKMFPSTWYV';
  end;
    
{$endif    }

end.
