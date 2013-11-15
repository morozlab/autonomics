{$ifdef osf1 }
  [ENVIRONMENT ('tabutil.pen'),
   INHERIT ('dep_dec.pen')]
  module tabutil;
{$endif }

{$ifdef vms }
  [ENVIRONMENT ('tabutil.pen'),
   INHERIT ('dep_vms.pen')]
  module tabutil;
{$endif }

{$ifdef Dos}
  unit tabutil;
  INTERFACE
  uses dep_tp;
{$endif }

{$ifdef Windows}
  unit tabutil;
  INTERFACE
  uses dep_win;
{$endif }

const
    maxcolumns      = 30;
    maxrows         = 20;

type
    namestring      = string80;
    rowcodetype     = string20; {maxrows}
    tabletype       = array[1..maxcolumns,0..maxrows] of real;
    itemtype        = record
                        ncol,nrow       : integer;
                        colname,rowname : namestring;
                        rowcode         : rowcodetype;
                        refpos          : integer;
                        table           : tabletype
                      end;

{$ifdef Borland}
function readmatrix(
                       var item     : itemtype;   {table item to be read}
                           itemname : namestring; {item identifier}
                           filename : string255;  {where to look}
                    ) : integer;                  {returns 0 if successful}

function scalematrix(
                       var item     : itemtype;   {item to be scaled}
                           scale    : itemtype    {item used fo scaling}
                     ) : integer;                 {returns 0 if successful}

{----------------------------------------------------------------------------}
implementation
{$endif}

function readint(var infile : text;
                 var num    : integer;
                     term   : char)  : integer;   {returns 0 if successful}
var s    : string255;
    c    : char;
begin
  s:='';
  repeat
    read(infile,c);
    if c<>term then s:=s+c;
  until (c=term) or eoln(infile);
  depend_err:='';
  num:=str2int(s);
  if depend_err='' then readint:=0 else readint:=1;
end;

function readreal(var infile : text;
                  var num    : real;
                      term   : char)  : integer;   {returns 0 if successful}
var s    : string255;
    c    : char;
    code : integer;
begin
  s:='';
  repeat
    read(infile,c);
    if c<>term then s:=s+c;
  until (c=term) or eoln(infile);
  depend_err:='';
  num:=str2real(s);
  if depend_err='' then readreal:=0 else readreal:=1;
end;

function readstr(var infile : text;

{$ifdef Borland}                           
                 var strg;
{$endif}                           
{$ifdef dec}                           
                 var strg   : varying[v1] of char;
{$endif}                           
                     term   : char)  : integer;   {returns 0 if successful}

var c    : char;
begin
  strg:='';
  repeat
    read(infile,c);
{$ifdef Borland}                           
    if c<>term then string(strg):=string(strg)+c;
{$endif}                           
{$ifdef dec}                           
    if c<>term then strg:=strg+c;
{$endif}                           
  until (c=term) or eoln(infile);
  readstr:=0;
end;



function readmatrix(
                       var item     : itemtype;   {table item to be read}
                           itemname : namestring; {item identifier}
                           filename : string255  {where to look}
                    ) : integer;              {returns 0 if successful}
var infile   : text;
    termchar : char;
    i,j      : integer;
    err      : integer;
    line     : string255;
begin
  fopen(infile,filename,'r');
  repeat
    readln(infile,line);
    if length(line)=0 then line:=' ';
    while (line[1]<>'[') and not eof(infile) do begin
      readln(infile,line);
      if length(line)=0 then line:=' ';
    end
  until eof(infile) or (indx(line,itemname)=1);
  if eof(infile) then readmatrix:=1 else
  begin
    readln(infile,line);
    if indx(line,'matrix')<>1 then readmatrix:=2 else
    begin
      err:=0;
      err:=err+readint(infile,item.ncol,' ');     readln(infile);
      err:=err+readint(infile,item.nrow,' ');     readln(infile);
      err:=err+readstr(infile,item.colname,' ');  readln(infile);
      err:=err+readstr(infile,item.rowname,' ');  readln(infile);
      err:=err+readstr(infile,item.rowcode,' ');  readln(infile);
      readln(infile,termchar);
      err:=err+readint(infile,item.refpos,' ');   readln(infile);
      repeat
        readln(infile,line);
      until (indx(line,'begin')=1) or eof(infile);
      if eof(infile) then err:=100;
      
      for j:=1 to item.nrow do begin
        for i:=1 to item.ncol do begin
          if not eof(infile) then 
	    err:=err+readreal(infile,item.table[i,j],termchar);
        end;
        if eof(infile) then err:=err+1 else readln(infile);
      end;
      for i:=1 to item.ncol do begin
        item.table[i,0]:=0;
	for j:=1 to item.nrow do begin
	  item.table[i,0]:=item.table[i,0]+item.table[i,j]/item.nrow
	end;
      end;
      
      if not eof(infile) then readln(infile,line);
      if indx(line,'end')<>1 then err:=200;
      if err=0 then readmatrix:=0 else readmatrix:=2+err;
    end;
  end;
  close(infile);
end;



function scalematrix(
                       var item     : itemtype;   {item to be scaled}
                           scale    : itemtype    {item used fo scaling}
                     ) : integer;                 {returns 0 if successful}
var sum   : real;
    i,j   : integer;
    achar : char;
    p2    : integer;

begin
 for i:=1 to item.ncol do begin
   sum:=0;
   for j:=1 to item.nrow do sum:=sum+item.table[i,j];
   for j:=1 to item.nrow do begin
     if item.table[i,j]=0 then item.table[i,j]:=1/sum;   {to be improved !!!}
     item.table[i,j]:=item.table[i,j]/sum;
     achar:=item.rowcode[j];
     p2:=indx(scale.rowcode,achar);                  {should be equal to j }
     item.table[i,j]:=ln(item.table[i,j]/scale.table[1,p2]);
   end;
   item.table[i,0]:=0;
   for j:=1 to item.nrow do begin
     item.table[i,0]:=item.table[i,0]+item.table[i,j]/item.nrow
   end;
 end;
 scalematrix:=0;
end;


{$ifdef Borland}
  begin
{$endif}
      
   end.

