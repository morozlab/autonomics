{$ifdef osf1}
  [INHERIT ('dep_dec.pen','seq.pen','ureadseq.pen','tabutil.pen')]
  program tmpred(input,output);
{$endif}

{$ifdef vms}
  [INHERIT ('dep_vms.pen','seq.pen','ureadseq.pen','tabutil.pen')]
  program tmpred(input,output);
{$endif}

{$ifdef Dos}
  program tmpred(input,output);
  uses dep_tp,ureadseq,seq,tabutil;
{$endif}

{$ifdef Windows}
  program tmpred(input,output);
  uses wincrt;dep_win,ureadseq,seq,tabutil;
{$endif}

const
    maxseqlen  = 6000;
    maxhelix   = 150;

type
    seqtype         = packed array[1..maxseqlen] of char;
    profiletype     = array[1..maxseqlen] of integer;
    featuretype     = record
                        position           : integer;
                        score              : integer;
                      end;
    helixtype       = record
                        center,nterm,cterm : featuretype;
                        sh_cterm,sh_nterm  : featuretype;
                        score              : integer;
                        nt_in              : boolean;
                      end;
    helixarray      = array[1..maxhelix div 2] of helixtype;
    cltagtype       = string40;
    cllinetype      = string255;

var
    tseq                            : seqtype;
    seqlen                          : longint;
    seqname,tablename               : string255;

    predname                        : string255;
    predfile                        : text;


    io_center,io_nt,io_ct,
    oi_center,oi_nt,oi_ct,aa_freq   : itemtype;

    iom,ion,ioc,
    oim,oin,oic                     : profiletype;
    io_score,oi_score               : profiletype;

    minwidth, maxwidth              : integer;

    io_helix,oi_helix               : helixarray;
    io_count,oi_count               : integer;

    thres1                          : integer;      {low orientational significance level}
    thres2                          : integer;      {high orientational significance level}
    thres3                          : integer;      {TM-existence significance level}
    thres4                          : integer;      {average orientation significance level}

    errcount                        : integer;

   { *** command line parameters ***  }
   interactflag : boolean; {user interaction ? }
   clinflag     : boolean; {sequence input filename specified ? }
   cloutflag    : boolean; {main output filename specified ? }
   clout2flag   : boolean; {auxiliary(2) output filename specified ? }
   clout3flag   : boolean; {auxiliary(3) output filename specified ? }
   clparflag    : boolean; {parameter (matrix) filename specified ? }
   clminflag    : boolean; {minimal helix length specified ?}
   clmaxflag    : boolean; {maximal helix length specified ?}
   htmlflag     : boolean; {HTML output required ?}
   inname,
   outname,
   outname2,
   outname3,
   parname      : string255;
   minhel,
   maxhel       : integer;

{$ifdef Borland} 
function min(i,j : integer) : integer;
begin if i<j then min:=i else min:=j end;

function max(i,j : integer) : integer;
begin if i>j then max:=i else max:=j end;
{$endif}

function findmax(var prf        : profiletype;
                     start,stop : integer;
                 var position   : integer) : integer;
var m,i : integer;
begin
  m:=prf[start];position:=start;
  for i:=start to stop do if prf[i]>m then begin
    m:=prf[i];
    position:=i;
  end;
  findmax:=m;
end;



{******** Command Line section *****************}

function get_cl_filename(cl : cllinetype; tag : cltagtype) : string255;
var p1,p2,p3 : integer;
    cl2      : cllinetype;
begin
  p1:=indx(cl,tag);
  cl2:=substr(cl,p1+length(tag),length(cl)-(p1+length(tag))+1);
  p2:=indx(cl2,c_sep); if p2=0 then p2:=length(cl2)+1;
  p3:=indx(cl2,' '); if p3=0 then p3:=length(cl2)+1;
  get_cl_filename:=substr(cl2,1,min(p2,p3)-1);
end;

function get_cl_int(cl : cllinetype; tag : cltagtype) : integer;
var p1,p2    : integer;
    cl2      : cllinetype;
    dummy    : string255;
begin
  p1:=indx(cl,tag);
  cl2:=substr(cl,p1+length(tag),length(cl)-(p1+length(tag))+1);
  p2:=indx(cl2,' '); if p2=0 then p2:=length(cl2)+1;
  dummy:=substr(cl2,1,p2-1);
  get_cl_int:=str2int(dummy);
end;


procedure process_command_line;
var cl   : cllinetype;
    cs   : char;
    incr : integer;
begin
  cl:=get_command_line;
  if indx(cl,c_sep+'def')>0 then interactflag:=false else interactflag:=true;
  if indx(cl,c_sep+'html')>0 then htmlflag:=true else htmlflag:=false;
  
  if indx(cl,c_sep+'in=')>0 then begin
    clinflag:=true;
    inname:=get_cl_filename(cl,c_sep+'in=');
  end else clinflag:=false;
  
  if indx(cl,c_sep+'out=')>0 then begin
    cloutflag:=true;
    outname:=get_cl_filename(cl,c_sep+'out=');
  end else cloutflag:=false;
  
  if indx(cl,c_sep+'out2=')>0 then begin
    clout2flag:=true;
    outname2:=get_cl_filename(cl,c_sep+'out2=');
  end else clout2flag:=false;
  
  if indx(cl,c_sep+'out3=')>0 then begin
    clout3flag:=true;
    outname3:=get_cl_filename(cl,c_sep+'out3=');
  end else clout3flag:=false;
  
  if indx(cl,c_sep+'par=')>0 then begin
    clparflag:=true;
    parname:=get_cl_filename(cl,c_sep+'par=');
  end else clparflag:=false;
  
  if indx(cl,c_sep+'min=')>0 then begin
    clminflag:=true;
    minhel:=get_cl_int(cl,c_sep+'min=');   
    if (minhel<11) or (minhel>40) then clminflag:=false;
  end else clminflag:=false;
  
  if indx(cl,c_sep+'max=')>0 then begin
    clmaxflag:=true;
    maxhel:=get_cl_int(cl,c_sep+'max=');   
    if (maxhel<11) or (maxhel>40) then clmaxflag:=false;
  end else clmaxflag:=false;

end;

{******** end of Command Line section *****************}


function getmatrices : integer;
begin
  if readmatrix(aa_freq,'[aa_freq]',tablename)<>0 then writeln('error in aa_freq');
  if readmatrix(io_center,'[io_center]',tablename)<>0 then writeln('error in io_center');
  if readmatrix(oi_center,'[oi_center]',tablename)<>0 then writeln('error in oi_center');
  if readmatrix(io_nt,'[io_nt]',tablename)<>0 then writeln('error in io_nt');
  if readmatrix(oi_nt,'[oi_nt]',tablename)<>0 then writeln('error in oi_nt');
  if readmatrix(io_ct,'[io_ct]',tablename)<>0 then writeln('error in io_ct');
  if readmatrix(oi_ct,'[oi_ct]',tablename)<>0 then writeln('error in oi_ct');
  if scalematrix(io_center,aa_freq)<>0 then writeln('error in io_center');
  if scalematrix(oi_center,aa_freq)<>0 then writeln('error in oi_center');
  if scalematrix(io_nt,aa_freq)<>0 then writeln('error in io_nt');
  if scalematrix(oi_nt,aa_freq)<>0 then writeln('error in oi_nt');
  if scalematrix(io_ct,aa_freq)<>0 then writeln('error in io_ct');
  if scalematrix(oi_ct,aa_freq)<>0 then writeln('error in oi_ct');
  getmatrices:=0;
end;

procedure make_profile(item : itemtype; var profile : profiletype);
var k,i : integer;
    p   : integer;
    m   : real;
begin
  for k:=1 to seqlen do begin
    m:=0;
    for i:=1 to item.ncol do begin
      p:=k+i-item.refpos;
      if (p>=1) and (p<=seqlen) then
        m:=m+item.table[i,indx(item.rowcode,tseq[p])];
    end;
    profile[k]:=round(m*100);    {for memories sake}
  end;
end;

{------------------------------------------------------------------------------}

procedure do_part1;
begin
  make_profile(io_center,iom);
  make_profile(io_nt,ion);
  make_profile(io_ct,ioc);
  make_profile(oi_center,oim);
  make_profile(oi_nt,oin);
  make_profile(oi_ct,oic);
end;

procedure output_part1;
var instring : string255;
    outfile  : text;
    tc       : char;
    i        : integer;
begin
  {tc:=chr(9);}
  tc:=',';
  if clout2flag then instring:=outname2
                else instring:='';
  if interactflag and not clout2flag then begin
    write('filename for auxiliary table output (<RETURN> if none) : ');
    readln(instring);
  end;
  if instring<>'' then begin
    fopen(outfile,instring,'w');
    writeln(outfile,'"TM-PREDICT auxiliary table output for sequence : ',seqname,'"');
    write(outfile,'"seq"':5,tc);
    write(outfile,'"i-o center"':12,tc,'"i-o Nter"':12,tc,'"i-o Cter"':12,tc);
    writeln(outfile,'"o-i center"':12,tc,'"o-i Nter"':12,tc,'"o-i Cter"':12);
    for i:=1 to seqlen do begin
      write(outfile,tseq[i]:5,tc,iom[i]:12,tc,ion[i]:12,tc,ioc[i]:12,tc);
      writeln(outfile,oim[i]:12,tc,oin[i]:12,tc,oic[i]:12);
    end;
    close(outfile);
  end;
end;

{------------------------------------------------------------------------------}


procedure make_curve(        le     : longint;
                      var    m,n,c  : profiletype;
                          minw,maxw : integer;
                      var score     : profiletype);
var
    i           : integer;
    minhw,maxhw : integer;
    dummy       : integer;
begin
  minhw:=minw div 2;
  maxhw:=(maxw+1) div 2;
  for i:=1 to le do begin
    if i<=minhw then score[i]:=0 else
    if i>(le-minhw) then score[i]:=0 else
    begin
      score[i]:=m[i]
                +findmax(n,max(i-maxhw,1),i-minhw,dummy)
                +findmax(c,i+minhw,min(i+maxhw,le),dummy);
    end;
  end;
end;


procedure do_part2;
begin
  make_curve(seqlen,iom,ion,ioc,minwidth,maxwidth,io_score);
  make_curve(seqlen,oim,oin,oic,minwidth,maxwidth,oi_score);
end;

procedure output_part2;
var
   instring : string255;
   outfile  : text;
   i        : integer;
   tc       : char;

begin
  {tc:=chr(9);}
  tc:=',';
  if clout3flag then instring:=outname3
                else instring:='';
  if interactflag and not clout3flag then begin
    write('filename for prediction graphics values (<RETURN> if none) : ');
    readln(instring);
  end;
  if instring<>'' then begin
    fopen(outfile,instring,'w');
    writeln(outfile,'# TM-PREDICT prediction graphics output for : ',seqname);
    writeln(outfile,'# minimal TM segment length : ',minwidth);
    writeln(outfile,'# maximal TM segment length : ',maxwidth);
    writeln(outfile,'# seq':3,tc,'i-o':10,tc,'o-i':10);
    for i:=1 to seqlen do begin
      writeln(outfile,tseq[i]:3,tc,io_score[i]:10,tc,oi_score[i]:10);
    end;
    close(outfile);
  end;
end;

{------------------------------------------------------------------------------}

function find_helix( var seq        : seqtype;
                         le         : longint;
                         minw,maxw  : integer;
                     var start      : integer;
                     var helix      : helixtype;
                     var s,m,n,c    : profiletype ) : boolean;
var
  minhw,maxhw : integer;
  found,done  : boolean;
  dummy       : integer;
  i,j         : integer;
begin
  minhw:=minw div 2;
  maxhw:=(maxw+1) div 2;
  i:=start;
  if i<=minhw then i:=minhw+1; {*** neu ***}
  found:=false;
  while (i<=le-minhw) and (not found) do begin
    if (s[i]=findmax(s,i-minhw,i+maxhw,dummy)) and (s[i]>0) then begin
      found:=true;
      helix.center.position:=i;
      helix.center.score:=m[i]; {determine center of TM-segment}

      helix.nterm.score:=findmax(n,max(i-maxhw,1),i-minhw,helix.nterm.position);  {optimal..}
      helix.cterm.score:=findmax(c,i+minhw,min(i+maxhw,le),helix.cterm.position); {..termini}

      j:=i-minhw;    {determine nearest N-terminus}
      done:=false;
      while (j-1>=1) and (j-1>=i-maxhw) and (not done) do begin
        if (n[j-1]>n[j]) then j:=j-1 else done:=true;
      end;
      helix.sh_nterm.score:=n[j];
      helix.sh_nterm.position:=j;

      j:=i+minhw;   {determine nearest C-terminus}
      done:=false;
      while (j+1<=le) and (j+1<=i+maxhw) and (not done) do begin
        if (c[j+1]>c[j]) then j:=j+1 else done:=true;
      end;
      helix.sh_cterm.score:=c[j];
      helix.sh_cterm.position:=j;
    end;
    i:=i+1;
  end;
  if found then begin
    start:=helix.sh_cterm.position+1;
    helix.score:=helix.center.score+helix.nterm.score+helix.cterm.score;
    find_helix:=true;
  end else begin
    start:=le;
    find_helix:=false;
  end;
end;


procedure do_part3;
var start : integer;
    helix : helixtype;
    found : boolean;
begin
  start:=1;
  io_count:=0;
  repeat
    found:=find_helix(tseq,seqlen,minwidth,maxwidth,start,helix,io_score,iom,ion,ioc);
    if found then begin
      io_count:=io_count+1;
      helix.nt_in:=true;
      io_helix[io_count]:=helix;
    end;
  until (start>=seqlen) or (io_count>(maxhelix div 2));
  start:=1;
  oi_count:=0;
  repeat
    found:=find_helix(tseq,seqlen,minwidth,maxwidth,start,helix,oi_score,oim,oin,oic);
    if found then begin
      oi_count:=oi_count+1;
      helix.nt_in:=false;
      oi_helix[oi_count]:=helix;
    end;
  until (start>=seqlen) or (oi_count>(maxhelix div 2));
end;

procedure output_part3(var f : text);
var i : integer;
begin
  if htmlflag then begin
    write(f,'Sequence: ',tseq[1],tseq[2],tseq[3],'...',
            tseq[seqlen-2],tseq[seqlen-1],tseq[seqlen]);
    writeln(f,',   length: ',seqlen:7,'<BR>');
    writeln(f,'Prediction parameters: TM-helix length between ',minwidth:2,
            ' and ',maxwidth:2);
    writeln(f,'<HR>');
    writeln(f,' ');
    writeln(f,'<H2> 1.) Possible transmembrane helices </H2>');
    writeln(f,'The sequence positions in brackets denominate the core region.<BR>');
    writeln(f,'Only scores above ',thres3:4,' are considered significant.<BR>');
    writeln(f,'<PRE>');
  end else begin
    writeln(f,'TMpred prediction output for : ',seqname); 
    writeln(f,' ');
    write(f,'Sequence: ',tseq[1],tseq[2],tseq[3],'...',
          tseq[seqlen-2],tseq[seqlen-1],tseq[seqlen]);
    writeln(f,'   length: ',seqlen:7);
    writeln(f,'Prediction parameters: TM-helix length between ',minwidth:2,
            ' and ',maxwidth:2);
    writeln(f,' ');
    writeln(f,' ');
    writeln(f,'1.) Possible transmembrane helices');
    writeln(f,'==================================');
    writeln(f,'The sequence positions in brackets denominate the core region.');
    writeln(f,'Only scores above ',thres3:4,' are considered significant.');
    writeln(f,' ');
  end;
  writeln(f,'Inside to outside helices : ',io_count:3,' found');
  if io_count>0 then begin
    writeln(f,'from ':11, 'to  ':11,'score':7,'center':7);
    for i:=1 to io_count do with io_helix[i] do begin
      writeln(f,nterm.position:4,' (',sh_nterm.position:4,')',
                cterm.position:4,' (',sh_cterm.position:4,')',
                score:7,center.position:7);
    end;
  end;
  writeln(f);
  writeln(f,'Outside to inside helices : ',oi_count:3,' found');
  if oi_count>0 then begin
    writeln(f,'from ':11, 'to  ':11,'score':7,'center':7);
    for i:=1 to oi_count do with oi_helix[i] do begin
      writeln(f,nterm.position:4,' (',sh_nterm.position:4,')',
                cterm.position:4,' (',sh_cterm.position:4,')',
                score:7,center.position:7);
    end;
  end;
  writeln(f,'');
  if htmlflag then writeln(f,'</PRE> <HR>')
              else writeln(f,'');
end;

procedure out4(var f : text; io,oi : integer);

var io_sign,oi_sign : string2;
    brack1,brack2   : char;

begin
  if (oi=0) then io_sign:='++'
  else if (io=0) then oi_sign:='++'
  else if io_helix[io].score-oi_helix[oi].score>thres2 then
  begin
    io_sign:='++';
    oi_sign:='  ';
  end
  else if io_helix[io].score-oi_helix[oi].score>thres1 then
  begin
    io_sign:=' +';
    oi_sign:='  ';
  end
  else if io_helix[io].score-oi_helix[oi].score<(-thres2) then
  begin
    io_sign:='  ';
    oi_sign:='++';
  end
  else if io_helix[io].score-oi_helix[oi].score<(-thres1) then
  begin
    io_sign:='  ';
    oi_sign:=' +';
  end
  else begin
    io_sign:='  ';
    oi_sign:='  ';
  end;

  if io=0 then begin
    write(f,' ':26,'|')
  end
  else with io_helix[io] do begin
    if score<thres3
      then begin brack1:='('; brack2:=')'; end
      else begin brack1:=' '; brack2:=' '; end;
    write(f,brack1,' ',nterm.position:4,'-',cterm.position:4,' ');
    write(f,'(',cterm.position-nterm.position+1:2,')');
    write(f,' ',score:4,' ',io_sign:2,' ',brack2,' |');
  end;
  if oi<>0 then with oi_helix[oi] do begin
    if score<thres3
      then begin brack1:='('; brack2:=')'; end
      else begin brack1:=' '; brack2:=' '; end;
    write(f,brack1,' ',nterm.position:4,'-',cterm.position:4,' ');
    write(f,'(',cterm.position-nterm.position+1:2,')');
    writeln(f,' ',score:4,' ',oi_sign:2,' ',brack2,' ');
  end else begin
    writeln(f)
  end;
end;

procedure output_part4(var f : text);
var iocount,oicount : integer;
begin
  iocount:=1; oicount:=1;
  if htmlflag then begin
    writeln(f,'<H2> 2.) Table of correspondences </H2>');
    writeln(f,'Here is shown, which of the inside->outside helices correspond');
    writeln(f,'to which of the outside->inside helices.<BR>');
    writeln(f,'<BLOCKQUOTE>');
    writeln(f,'Helices shown in brackets are considered insignificant.<BR>');
    writeln(f,'A "+"-symbol indicates a preference of this orientation.<BR>');
    writeln(f,'A "++"-symbol indicates a strong preference of this orientation.<BR>');
    writeln(f,'</BLOCKQUOTE>');
    writeln(f,'<PRE>');
  end else begin
    writeln(f,' ');
    writeln(f,'2.) Table of correspondences');
    writeln(f,'============================');
    writeln(f,'Here is shown, which of the inside->outside helices correspond');
    writeln(f,'to which of the outside->inside helices.');
    writeln(f,'  Helices shown in brackets are considered insignificant.');
    writeln(f,'  A "+"  symbol indicates a preference of this orientation.');
    writeln(f,'  A "++" symbol indicates a strong preference of this orientation.');
  end;
  writeln(f,' ');
  writeln(f,'inside->outside ':27,'|',' outside->inside');

  repeat
    if (iocount>io_count) and (oicount>oi_count) then
    begin
    end
    else if (iocount>io_count) and (oicount<=oi_count) then
    begin
      out4(f,0,oicount);
      oicount:=oicount+1;
    end
    else if (oicount>oi_count) and (iocount<=io_count) then
    begin
      out4(f,iocount,0);
      iocount:=iocount+1;
    end
    else if io_helix[iocount].center.position<=oi_helix[oicount].center.position then
    begin
      if oi_helix[oicount].sh_nterm.position<=io_helix[iocount].sh_cterm.position
         then begin
           out4(f,iocount,oicount);
           iocount:=iocount+1;
           oicount:=oicount+1;
         end else begin
           out4(f,iocount,0);
           iocount:=iocount+1;
         end;
    end else
    begin
      if io_helix[iocount].sh_nterm.position<=oi_helix[oicount].sh_cterm.position
         then begin
           out4(f,iocount,oicount);
           iocount:=iocount+1;
           oicount:=oicount+1;
         end else begin
           out4(f,0,oicount);
           oicount:=oicount+1;
         end;
    end;
  until (iocount>io_count) and (oicount>oi_count);
  writeln(f,'');
  if htmlflag then writeln(f,'</PRE><HR>')
              else writeln(f,'');
end;

procedure part5_findnext(    helix : helixarray;           {the array to search in}
                         var ct    : integer;              {where to start searching}
                             maxct : integer;              {no of TM-segs in array}
                             ntpos : integer;              {earliest pos. of N-term}
                             thres : integer);             {threshold for ignoring}
                                                           {**returns ct=0 if not found**}
begin
  while ( (helix[ct].score<thres) or (helix[ct].sh_nterm.position<ntpos) ) and (ct<maxct)
    do ct:=ct+1;
  if (helix[ct].score<thres) or (helix[ct].sh_nterm.position<ntpos)
    then ct:=0
end;


procedure part5_model(     ori_io   : boolean;       {is N-term inside?}
                       var model    : helixarray;    {put model here}
                       var model_ct : integer);      {number of segments in model}
var
    io_ct,oi_ct  : integer;
    nt_pos       : integer;
    done         : boolean;
begin
  io_ct:=1; oi_ct:=1;
  nt_pos:=1;
  done:=false;
  model_ct:=0;

  repeat
    if ori_io then begin
      part5_findnext(io_helix,io_ct,io_count,nt_pos,thres3);
      if io_ct=0
      then done:=true
      else begin
        model_ct:=model_ct+1;
        model[model_ct]:=io_helix[io_ct];
        nt_pos:=model[model_ct].sh_cterm.position+1;
      end;
    end else begin
      part5_findnext(oi_helix,oi_ct,oi_count,nt_pos,thres3);
      if oi_ct=0
      then done:=true
      else begin
        model_ct:=model_ct+1;
        model[model_ct]:=oi_helix[oi_ct];
        nt_pos:=model[model_ct].sh_cterm.position+1;
      end;
    end;
    ori_io:=not ori_io;
  until done;
end;

procedure out_part5_model( var f       : text;
                               model   : helixarray;
                               count   : integer;
                               totscore: integer);
var i : integer;
begin
  writeln(f,count:2,' strong transmembrane helices, total score : ',totscore);
    writeln(f,'#':2,' ','from':4,' ','to':4,' ','length ','score':5,' ','orientation');
  for i:=1 to count do with model[i] do begin
    write(f,i:2,' ',nterm.position:4,' ',cterm.position:4,' ');
    write(f,'(',cterm.position-nterm.position+1:2,')   ',score:5,' ');
    if nt_in then writeln(f,'i-o') else writeln(f,'o-i');
  end;
  writeln(f);
end;

procedure output_part5(var f : text);
var in_model,out_model         : helixarray;
    in_mod_count,out_mod_count : integer;
    in_mod_score,out_mod_score : integer;
    i                          : integer;
begin
  if htmlflag then begin
    writeln(f,'<H2> 3.) Suggested models for transmembrane topology</H2>');
    writeln(f,'These suggestions are purely speculative and should be used with');
    writeln(f,'<B>extreme caution</B> since they are based on the assumption that');
    writeln(f,'all transmembrane helices have been found.<BR>');
    writeln(f,'In most cases, the Correspondence Table shown above or the');
    writeln(f,'prediction plot that is also created should be used for the');
    writeln(f,'topology assignment of unknown proteins.<BR>');
    writeln(f,'<PRE>');
  end else begin
    writeln(f,'3.) Suggested models for transmembrane topology');
    writeln(f,'===============================================');
    writeln(f,'These suggestions are purely speculative and should be used with');
    writeln(f,'EXTREME CAUTION since they are based on the assumption that');
    writeln(f,'all transmembrane helices have been found.');
    writeln(f,'In most cases, the Correspondence Table shown above or the');
    writeln(f,'prediction plot that is also created should be used for the');
    writeln(f,'topology assignment of unknown proteins.');
  end;
  writeln(f,'');
  part5_model(true,in_model,in_mod_count);
  in_mod_score:=0;
  for i:=1 to in_mod_count do in_mod_score:=in_mod_score+in_model[i].score;
  part5_model(false,out_model,out_mod_count);
  out_mod_score:=0;
  for i:=1 to out_mod_count do out_mod_score:=out_mod_score+out_model[i].score;
  writeln(f,'2 possible models considered, only significant TM-segments used');
  writeln(f);
  if in_mod_count<>out_mod_count then begin
    writeln(f,'*** the models differ in the number of TM-helices ! ***');
    writeln(f);
  end;
  if in_mod_score>out_mod_score then begin
    if in_mod_count=0 then begin
      writeln(f,'!!! probably no transmembrane protein - no possible model found !!!');
    end
    else begin
      if (in_mod_score-out_mod_score/in_mod_count)>thres4   {to be improved}
      then  writeln(f,'-----> STRONGLY prefered model: N-terminus inside')
      else  writeln(f,'-----> slightly prefered model: N-terminus inside');
      out_part5_model(f,in_model,in_mod_count,in_mod_score);
      writeln(f,'------> alternative model');
      out_part5_model(f,out_model,out_mod_count,out_mod_score);
    end;
  end else begin
    if out_mod_count=0 then begin
      writeln(f,'!!! probably no transmembrane protein - no possible model found !!!');
    end
    else begin
      if (out_mod_score-in_mod_score/out_mod_count)>thres4   {to be improved}
      then  writeln(f,'-----> STRONGLY prefered model: N-terminus outside')
      else  writeln(f,'-----> slightly prefered model: N-terminus outside');
      out_part5_model(f,out_model,out_mod_count,out_mod_score);
      writeln(f,'------> alternative model');
      out_part5_model(f,in_model,in_mod_count,in_mod_score);
    end;
  end;
  if htmlflag then writeln(f,'</PRE><HR>');
end;

{------------------------------------------------------------------------------}


function getsequence : integer;
var i : integer;
begin
  seqlen:=maxseqlen;
  sq_read_file(tseq,seqlen,seqname);
  if (seqlen=0) or (seqlen>=maxseqlen) then getsequence:=1 
  else begin 
    getsequence:=0;
    for i:=1 to seqlen do begin
      tseq[i]:=upcase(tseq[i]);
{***   if index('ACDEFGHIKLMNPQRSTVWY',tseq[i])=0 then begin   }
{***     tseq[i]:='A';                                         }
{***     writeln('Invalid amino acid replaced by Alanine');    }
{***   end;                                                    } 
    end;
  end;
end;

begin
  {screensize.y:=80;      }
  thres1:=80;
  thres2:=200;
  thres3:=500;
  thres4:=80;
  process_command_line;
  if interactflag then begin
    writeln('Prediction of transmembrane helices by the matrix method');
    writeln('different probabilities for inside-to-outside and outside-');
    writeln('to-inside helices are calculated');
    writeln('');
    writeln('');
  end;
  if clparflag then tablename:=parname
               else tablename:='matrix.tab';
  if not fexist(tablename) then begin
    tablename:=get_logical('TMPREDDIR')+d_sep+'matrix.tab';
  end;
  if interactflag then begin
    while (not fexist(tablename)) and (errcount<3) do begin
      errcount:=errcount+1;
      write('matrix file:',tablename,' does not exist, enter valid name : ');
      readln(tablename);
    end;
  end;
  
  if clinflag then seqname:=inname
              else seqname:='tmpred.in';
  if interactflag and not clinflag then begin
    write('Name of Sequence (in arbitrary format) : ');
    readln(seqname);
  end;
  errcount:=0;
  if interactflag then while (not fexist(seqname)) and (errcount<3) do begin
    errcount:=errcount+1;
    write('sequence file does not exist, try again : ');
    readln(seqname);
  end else begin
    if not fexist(seqname) then writeln('sequence file not found');
  end;
  
  if getsequence<>0 then writeln(chr(7),'error in the sequence');
  if interactflag then writeln('please wait a while ....');
  if getmatrices<>0 then writeln(chr(7),'error in the matrices');
  do_part1;                                        {the auxiliary table}
  output_part1;
  
  if interactflag then begin
    if clminflag then minwidth:=minhel else begin
      write('minimal length of transmembrane sequence : ');
      readln(minwidth);
    end;
    if clmaxflag then maxwidth:=maxhel else begin
      write('maximal length of transmembrane sequence : ');
      readln(maxwidth);
    end;
  end else begin
    minwidth:=17;
    maxwidth:=35;
    if clminflag then minwidth:=minhel;
    if clmaxflag then maxwidth:=maxhel;
  end;
  
  do_part2;                                        {the predicition graphics}
  output_part2;
  
  do_part3; {the actual prediction}
  if cloutflag then predname:=outname
               else predname:='tmpred.out';
  if interactflag and not cloutflag then begin
    write('file for predicion output :');
    readln(predname);
  end;
  
  fopen(predfile,predname,'w');
  output_part3(predfile);
  output_part4(predfile);
  output_part5(predfile);
  close(predfile);
{$ifdef Windows}
  donewincrt; }
{$endif}
end.

