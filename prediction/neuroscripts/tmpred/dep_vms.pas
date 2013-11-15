[ENVIRONMENT('dep_vms.pen')]

module dep_vms;    {implementation dependent parts - DECpascal version}

const pi       = 3.1415926;

const max_longint  = maxint;
      max_shortint = 32768;
      depend_mode  = 'dec';
      c_sep        = '/';  { command line separator}
      d_sep        = ':';  { directory separator   }
                           { put between a logical and a filename}

type string2   = varying[2]   of char;
     string3   = varying[3]   of char;
     string4   = varying[4]   of char;
     string5   = varying[5]   of char;
     string6   = varying[6]   of char;
     string7   = varying[7]   of char;
     string8   = varying[8]   of char;
     string9   = varying[9]   of char;
     string10  = varying[10]  of char;
     string20  = varying[20]  of char;
     string26  = varying[26]  of char;
     string30  = varying[30]  of char;
     string40  = varying[40]  of char;
     string50  = varying[50]  of char;
     string60  = varying[60]  of char;
     string70  = varying[70]  of char;
     string80  = varying[80]  of char;
     string90  = varying[90]  of char;
     string100 = varying[100] of char;
     string128 = varying[128] of char;
     string255 = varying[255] of char;
     longint   = integer;              { 4 byte integer }


const rand_m   = 100000000;
      rand_m1  = 10000;
      rand_b   = 31415821;

var rand_a     : array[0..54] of longint;
    rand_j     : integer;
    depend_err : string255;
    unixmode   : boolean;

procedure lib$get_foreign(var s : varying[s1] of char); extern;

function get_logical(logical : string255) : string255;
begin
  get_logical:=logical;
end;

function upcase(c : char) : char;
begin
  if (ord(c)<=122) and (ord(c)>=97) then upcase:=chr(ord(c)-32)
  else upcase:=c;
end;

function lowcase(c : char) : char;
begin
  if (ord(c)<=90) and (ord(c)>=65) then lowcase:=chr(ord(c)+32)
  else lowcase:=c;
end;
  
function indx(mainstr: string255; substr : string255) : integer;
begin
  if mainstr='' then indx:=0 else indx:=index(mainstr,substr);
end;

function get_command_line : string255;
var cl : string255;
    i  : integer;
begin
  lib$get_foreign(cl);
  for i:=1 to length(cl) do cl[i]:=lowcase(cl[i]);
  get_command_line:=cl;
end;

function fexist(fname:string255) : boolean;
var dummy : text;
begin
  open(dummy,fname,readonly,error:=continue);
  if status(dummy)=0 then begin
    fexist:=true;
    close(dummy);
  end else begin
    fexist:=false;
  end;
end;

procedure ureadln(var f : text; var line:varying[s1] of char);
var done : boolean;
    c    : char;
    p1   : integer;
begin
  readln(f,line);
  if line<>'' then if line[length(line)]=chr(13) 
        then line:=substr(line,1,length(line)-1);
end;

function getunixmode(fname : string255) : boolean;
var line : string255;
    f    : text;
begin
  open(f,fname,readonly);
  reset(f);
  readln(f,line);
  if line[length(line)]=chr(13) then getunixmode:=false
                                else getunixmode:=true;
  close(f);
end;

procedure uwriteln(var f : text; line : string255);
begin
  if unixmode then writeln(f,line)
              else writeln(f,line+chr(13));
end;
 
function getTEMPname(s : string255) : string255;
begin
  getTEMPname:=s+'~';
end;
 
procedure fopen(var f :text; fname: string255; acc : char); {acc = r,w,a }
var line    : string255;
    c       : char;
    uxmode  : boolean;
    tmpname : string255;
    tmpfile : text;
    
begin
  if fname='' then fname:='tt:';
  if (acc='w') or (acc='W') then
    begin
      open(f,fname,new,record_length:=9096);
      rewrite(f);
    end
  else if (acc='a') or (acc='A') then
    begin
      open(f,fname,old,error:=continue);
      extend(f);
    end
  else
    begin
      open(f,fname,readonly);
      reset(f);
      readln(f,line);
      reset(f);
      if line[length(line)]=chr(13) then uxmode:=false else uxmode:=true;
      if not uxmode then begin
        tmpname:=getTEMPname(fname);
        open(tmpfile,tmpname);
        rewrite(tmpfile);
        while not eof(f) do begin
          ureadln(f,line);
          writeln(tmpfile,line);
        end;
        close(tmpfile);
        close(f);
        open(f,tmpname,readonly);
        reset(f);
        if status(f)<>0 then depend_err:='could not convert unix format';
      end;
    end;
end;


function mult(p,q: longint) : longint;
var p1,p0,q1,q0 : longint;
begin
  p1:=p div rand_m1; p0:=p mod rand_m1;
  q1:=q div rand_m1; q0:=q mod rand_m1;
  mult:=(((p0*q1+p1*q0) mod rand_m1) * rand_m1 +p0*q0) mod rand_m;
end;

procedure rand_init(seed : integer);
begin
  rand_a[0]:=seed; rand_j:=0;
  repeat
    rand_j:=rand_j+1;
    rand_a[rand_j]:=(mult(rand_b,rand_a[rand_j-1])+1) mod rand_m
  until rand_j=54;
end;

function rand_int(limit : integer) : integer; {from 0 to limit-1}
begin
  rand_j:=(rand_j+1) mod 55;
  rand_a[rand_j] :=(rand_a[(rand_j+23) mod 55]
                 +  rand_a[(rand_j+54) mod 55]) mod rand_m;
  rand_int:=((rand_a[rand_j] div rand_m1)*limit) div rand_m1;
end;

function rand_real : real;   { from 0(incl) to 1(excl) }
begin
  rand_real:=rand_int(32000) / 32000;
end;


function int2str(i : integer; dig : integer) : string255;
var s : string255;
begin
  writev(s,i:dig);
  int2str:=s;
end;

function real2str(r : real; dig1,dig2 : integer) : string255;
var s : string255;
begin
  writev(s,r:dig1:dig2);
  real2str:=s;
end;

function str2int(s : string80) : integer;
var i : integer;
begin
  readv(s,i,ERROR:=CONTINUE);
  if statusv<>0 then begin
    depend_err:='str2int';
    str2int:=0;
  end else
   str2int:=i;
end;

function str2real(s : string80) : real;
var r : real;
begin
  readv(s,r,ERROR:=CONTINUE);
  if statusv<>0 then begin
    depend_err:='str2real';
    str2real:=0;
  end else
   str2real:=r;
end;

to begin do begin
  rand_init(1);
  depend_err := '';
  unixmode:=true;
end;

end.
