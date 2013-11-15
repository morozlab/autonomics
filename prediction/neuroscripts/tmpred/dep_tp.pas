unit dep_tp;     {implementation dependent parts - turboPascal version}

interface
uses dos;

const pi          = 3.1415926;
     max_longint  = maxlongint;
     max_shortint = maxint;
     depend_mode  = 'tp';
     c_sep        = '/';  { command line separator}
     d_sep        = '\';  { directory separator   }
                          { put between a logical and a filename}

type string2     = string[2];
     string3     = string[3];
     string4     = string[4];
     string5     = string[5];
     string6     = string[6];
     string7     = string[7];
     string8     = string[8];
     string9     = string[9];
     string10    = string[10];
     string20    = string[20];
     string26    = string[26];
     string30    = string[30];
     string40    = string[40];
     string50    = string[50];
     string60    = string[60];
     string70    = string[70];
     string80    = string[80];
     string90    = string[90];
     string100   = string[100];
     string128   = string[128];
     string255   = string[255];

var depend_err : string40;
    unixmode   : boolean; {determines whether uwriteln writes unix or dos}

procedure fopen(var f :text;fname : string255; acc : char);{acc = r,w,a }
function fexist(fname : string255) : boolean;
function indx (mainstring:string255; substring :string255)      : integer;
function substr(var original; pos1,count : integer) : string255;
procedure rand_init(seed : integer);
function rand_int(limit  : integer): integer;     {0 <= rand < limit }
function rand_real : real;                        {0 <= rand <1 }
procedure ureadln(var f : text; var line);      {reads unix and dos style}
function getunixmode(fname : string255 ) : boolean;
procedure uwriteln(var f : text; line : string255); {depends on 'unixmode'}
function int2str(int : integer; dig : integer)      : string255;
function real2str(rl : real; dig1,dig2 : integer)   : string255;
function str2int(var s) : integer;
function str2real(var s) : real;
function get_command_line : string255;  {returns all arguments but prognam}
function get_logical(logical : string255) : string255;
function lowcase (c:char) : char;
{function upcase (c:char) : char;                        (exists anyway)}

implementation

const rand_m   : longint = 100000000;
      rand_m1  : longint = 10000;
      rand_b   : longint = 31415821;

var rand_a    : array[0..54] of longint;
    rand_j    : integer;
    filecount : integer;

function lowcase (c:char) : char;
begin
  if (ord(c)<=90) and (ord(c)>=65) then lowcase:=chr(ord(c)+32)
  else lowcase:=c;
end;

function get_logical(logical : string255) : string255;
begin
  get_logical:=getenv(logical);
end;

function get_command_line : string255;
var cl : string255;
    i  : integer;
begin
  cl:='';
  if paramcount>0 then
    for i:=1 to paramcount do cl:=cl+paramstr(i);
  get_command_line := cl;
end;

function fexist;
var f : text;
begin
  assign(f,fname);
  {$I-}
  reset(f);
  {$I+}
  if ioresult=0 then fexist:=true else fexist:=false;
end;

function getTEMPname(s : string255) : string255;
begin
  s[length(s)]:='~';
  getTEMPname:=s;
end;

procedure fopen(var f :text; fname: string255; acc : char); {acc = r,w,a }
var tmpfile : text;
    tmpname : string255;
    line    : string255;
    c       : char;
    buf1,
    buf2    : packed array[1..6*1024] of char;
    label l_exit;
begin
  if (acc='w') or (acc='W') then
    begin
      assign(f,fname);
      rewrite(f);
    end
  else if (acc='a') or (acc='A') then
    begin
      assign(f,fname);
      append(f);
    end
  else
    begin
      assign(f,fname);
      reset(f);
      repeat read(f,c) until eoln(f) or (c=chr(10));
      reset(f);
      if c=chr(10) then uxmode:=true else uxmode:=false;
      if uxmode then begin
        {$I-}
        settextbuf(f,buf1);
        tmpname:=getTEMPname(fname);
        assign(tmpfile,tmpname);
        rewrite(tmpfile);
        settextbuf(f,buf2);
        while not eof(f) do begin
          ureadln(f,line);
          if ioresult<>0 then goto l_exit;
          writeln(tmpfile,line);
        end;
        close(tmpfile);
        close(f);
        assign(f,tmpname);
        reset(f);
        {$I+}
      l_exit:
        if ioresult<>0 then depend_err:='could not convert unix format';
      end;
    end;
end;


function indx;
begin
  if mainstring='' then indx:=0 else
    indx:=pos(string(substring),string(mainstring));
end;

function substr(var original;pos1,count : integer) : string255;
begin
  substr:=copy(string(original),pos1,count);
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

procedure ureadln(var f : text; var line);
type string1 = string[1];
var done : boolean;
    c    : char;
    lin  : string255;
begin
  done:=false;
  string1(line):='';
  repeat
    read(f,c);
    if c=#10 then done:=true
             else if c<>#13 then string(line):=string(line)+c;
    if eoln(f) then begin
      readln(f);
      done:=true;
    end;
  until done;
end;

function getunixmode(fname : string255) : boolean;
var c : char;
    f : text;
begin
  assign(f,fname);
  reset(f);
  repeat
    read(f,c);
  until (c=chr(10)) or eoln(f);
  if c=chr(10) then getunixmode:=true
               else getunixmode:=false;
  close(f);
end;

procedure uwriteln(var f : text; line : string255);
begin
  if unixmode then write(f,line,chr(10))
              else writeln(f,line);
end;

function int2str(int : integer; dig : integer) : string255;
var s : string255;
begin
  str(int:dig,s);
  int2str:=s;
end;

function real2str(rl : real; dig1,dig2 : integer) : string255;
var s : string255;
begin
  str(rl:dig1:dig2,s);
  real2str:=s;
end;

function str2real(var s) : real;
var r   : real;
    err : integer;
begin
  val(string(s),r,err);
  if err<>0 then begin
    depend_err:='str2real';
    str2real:=0;
  end
  else str2real:=r;
end;

function str2int(var s) : integer;
var i   : integer;
    err : integer;
begin
  val(string(s),i,err);
  if err<>0 then begin
    depend_err:='str2real';
    str2int:=0;
  end
  else str2int:=i;
end;

begin
  depend_err:='';
  filecount:=100;
  rand_init(1);
  unixmode:=false;
end.
