program preprocess(input,output);
var line,
    newline : varying[255] of char;
    p1,p2   : integer;
    buffer  : varying[255] of char;
    infile  : text;
    outfile : text;
    inname  : varying[255] of char;
    outname : varying[255] of char;

begin
  if argc=2 then begin
    argv(1,buffer);
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
    inname:='INPUT';
    outname:='OUTPUT';
  end;
  open(infile,inname,readonly);
  reset(infile);
  open(outfile,outname,error:=continue);
  rewrite(outfile,error:=continue);
  while not eof(infile) do begin
    readln(infile,line);
    if index(line,'{$ifdef')>0 then begin
      p1:=index(line,'{$ifdef');
      p2:=index(line,'}');
      if p1>2 then newline:=substr(line,1,p1-1)
              else newline:='';
      newline:=newline+'#ifdef';
      newline:=newline+substr(line,p1+7,p2-p1-7);
    end 
    else if index(line,'{$endif')>0 then begin
      p1:=index(line,'{$endif');
      p2:=index(line,'}');
      if p1>2 then newline:=substr(line,1,p1-1)
              else newline:='';
      newline:=newline+'#endif';
      newline:=newline+substr(line,p1+7,p2-p1-7);
    end
    else newline:=line;

    writeln(outfile,newline);
  end;
end.
