program preprocess(input,output);



{ preprocessing program to evaluate conditional compilation              }

{ in turbo-pascal style. Only simple (non-nested) commands are accepted. }

{ commands removed are padded with empty comment lines in order to       }

{ preserve consistent line numbering                                     } 



var line,

    tag     : varying[255] of char;

    p1,p2   : integer;

    buffer  : varying[255] of char;

    infile  : text;

    outfile : text;

    inname  : varying[255] of char;

    outname : varying[255] of char;

    ndef    : integer;

    define  : array[1..10] of varying[255] of char;

    i       : integer;

    found   : boolean;

    copymode: boolean;

    label l_exit;

   



procedure lib$get_foreign(var cl : varying[s1] of char); extern;



begin

  ndef:=2;

  define[1]:='dec';

  define[2]:='vms';

  lib$get_foreign(buffer);

  writeln(buffer);

  if length(buffer)>0 then begin

    inname:=buffer;

    p1:=length(inname);

    while (inname[p1]<>'.') and (p1>1) do p1:=p1-1;

    if p1=1 then begin

      outname:=inname+'.p';

      inname:=inname+'.pas';

    end else begin

      outname:=substr(inname,1,p1)+'p';

    end;

  end else begin

   writeln('call this program either:  pp_vms myfile.pas');

   writeln('or:                        pp_vms myfile');

   writeln('in order to get myfile.p as processed output');

   goto l_exit;

  end;

  open(infile,inname,readonly);

  reset(infile);

  open(outfile,outname);

  rewrite(outfile);

  copymode:=true;

  while not eof(infile) do begin

    readln(infile,line);

    if index(line,'{$ifdef')>0 then begin

      p1:=index(line,'{$ifdef');

      p2:=index(line,'}');

      tag:=substr(line,p1+7,p2-p1-7);

      found:=false;

      for i:=1 to ndef do if (index(tag,define[i])>0) then found:=true;

      if found=false then copymode:=false else copymode:=true;

      line:='{}';

    end 

    else if index(line,'{$endif')>0 then begin

      copymode:=true;

      line:='{}';

    end;

    if copymode then writeln(outfile,line) else writeln(outfile,'{}');

  end;

  l_exit:

end.

