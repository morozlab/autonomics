{Batch mode TM-predict}

{$ifdef osf1}
  [INHERIT ('dep_dec.pen','seq.pen','ureadseq.pen','tabutil.pen')]
  program btmpred(input,output);
{$endif}

{$ifdef vms}
  [INHERIT ('dep_vms.pen','seq.pen','ureadseq.pen','tabutil.pen')]
  program btmpred(input,output);
{$endif}

{$ifdef Dos}
  program btmpred(input,output);
  uses dep_tp,ureadseq,seq,tabutil;
{$endif}

{$ifdef Windows}
  program btmpred(input,output);
  uses wincrt;dep_win,ureadseq,seq,tabutil;
{$endif}


const
{$ifdef Windows}
    maxseqlen  = 1100;
{$endif}
{$ifdef Dos}
    maxseqlen  = 2100;
{$endif}
{$ifdef dec}
    maxseqlen  = 30000;
{$endif}

maxhelix   = 50;

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
    assocarray      = array[1..maxhelix div 2] of integer;

var
    tseq                        : seqtype;
    seqlen                          : longint;
    seqname,tablename               : string255;
    outname1,outname2               : string255;
    seqid                           : string255;

    io_center,io_nt,io_ct,
    oi_center,oi_nt,oi_ct,aa_freq   : itemtype;

    iom,ion,ioc,
    oim,oin,oic                     : profiletype;
    io_score,oi_score               : profiletype;

    minwidth, maxwidth              : integer;

    io_helix,oi_helix               : helixarray;
    io_count,oi_count               : integer;
    io_assoc,oi_assoc               : assocarray;

    thres1                          : integer;      {low orientational significance level}
    thres2                          : integer;      {high orientational significance level}
    thres3                          : integer;      {TM-existence significance level}
    thres4                          : integer;      {average orientation significance level}

    errcount                        : integer;
    no_seq,max_no_seq               : integer;
    seqerr                          : integer;
    outdev1,outdev2                 : text;

    i                               : integer;


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
        m:=m+item.table[i,indx(item.rowcode,upcase(tseq[p]))];
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
  tc:=chr(9);
  write('filename for auxiliary table output (<RETURN> if none) : ');
  readln(instring);
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
  tc:=chr(9);
  write('filename for prediction graphics values (<RETURN> if none) : ');
  readln(instring);
  if instring<>'' then begin
    fopen(outfile,instring,'w');
    writeln(outfile,'TM-PREDICT prediction graphics output for : ',seqname);
    writeln(outfile,'minimal TM segment length : ',minwidth);
    writeln(outfile,'maximal TM segment length : ',maxwidth);
    writeln(outfile,'seq':3,tc,'i-o':10,tc,'o-i':10);
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
    null_helix : helixtype;
begin
  start:=1;
  io_count:=0;
  null_helix.score:=-maxint;
  repeat
    found:=find_helix(tseq,seqlen,minwidth,maxwidth,start,helix,io_score,iom,ion,ioc);
    io_count:=io_count+1;
    if found then begin
      helix.nt_in:=true;
      io_helix[io_count]:=helix;
    end else begin
      io_helix[io_count]:=null_helix;
    end;
  until (start>=seqlen) or (io_count>(maxhelix div 2));
  start:=1;
  oi_count:=0;
  repeat
    found:=find_helix(tseq,seqlen,minwidth,maxwidth,start,helix,oi_score,oim,oin,oic);
    oi_count:=oi_count+1;
    if found then begin
      helix.nt_in:=false;
      oi_helix[oi_count]:=helix;
    end else begin
      oi_helix[oi_count]:=null_helix;
    end;
  until (start>=seqlen) or (oi_count>(maxhelix div 2));
end;

procedure output_part3;
var i : integer;
begin
  writeln('Inside to outside helices : ',io_count:3,' found');
  if io_count>0 then begin
    writeln('from':16, 'to':16,'score':10,'center':10);
    for i:=1 to io_count do with io_helix[i] do begin
      writeln(nterm.position:7,'(',sh_nterm.position:7,')',
              cterm.position:7,'(',sh_cterm.position:7,')',
              score:10,center.position:10);
    end;
  end;
  writeln;
  writeln('Outside to inside helices : ',oi_count:3,' found');
  if oi_count>0 then begin
    writeln('from':16, 'to':16,'score':10,'center':10);
    for i:=1 to oi_count do with oi_helix[i] do begin
      writeln(nterm.position:7,'(',sh_nterm.position:7,')',
              cterm.position:7,'(',sh_cterm.position:7,')',
              score:10,center.position:10);
    end;
  end;
  writeln;
end;

procedure out4(io,oi : integer);

var io_sign,oi_sign : string2;
    brack1,brack2   : char;

begin
  if (oi=0) then io_sign:='++'
  else if (io=0) then oi_sign:='++'
  else if io_helix[io].score-oi_helix[oi].score>thres2 then
  begin
    io_sign:='++';
    oi_sign:='--';
  end
  else if io_helix[io].score-oi_helix[oi].score>thres1 then
  begin
    io_sign:=' +';
    oi_sign:=' -';
  end
  else if io_helix[io].score-oi_helix[oi].score<(-thres2) then
  begin
    io_sign:='--';
    oi_sign:='++';
  end
  else if io_helix[io].score-oi_helix[oi].score<(-thres1) then
  begin
    io_sign:=' -';
    oi_sign:=' +';
  end
  else begin
    io_sign:='  ';
    oi_sign:='  ';
  end;
  if io=0 then begin
    write(' ':25,'|')
  end
  else with io_helix[io] do begin
    if score<thres3
      then begin brack1:='('; brack2:=')'; end
      else begin brack1:=' '; brack2:=' '; end;
    write(brack1,' ',nterm.position:3,'-',cterm.position:3,' ');
    write('(',cterm.position-nterm.position+1:2,')');
    write(' ',score:4,' ',io_sign:2,' ',brack2,' |');
  end;
  if oi<>0 then with oi_helix[oi] do begin
    if score<thres3
      then begin brack1:='('; brack2:=')'; end
      else begin brack1:=' '; brack2:=' '; end;
    write(brack1,' ',nterm.position:3,'-',cterm.position:3,' ');
    write('(',cterm.position-nterm.position+1:2,')');
    writeln(' ',score:4,' ',oi_sign:2,' ',brack2,' ');
  end else begin
    writeln
  end;
end;

procedure do_part4;
var iocount,oicount : integer;
begin
  iocount:=1; oicount:=1;
  repeat
    if (iocount>io_count) and (oicount>oi_count) then begin end
    else if (iocount>io_count) and (oicount<=oi_count) then
    begin
      oi_assoc[oicount]:=0;
      oicount:=oicount+1;
    end else if (oicount>oi_count) and (iocount<=io_count) then
    begin
      io_assoc[iocount]:=0;
      iocount:=iocount+1;
    end else if io_helix[iocount].center.position<=oi_helix[oicount].center.position then
    begin
      if oi_helix[oicount].sh_nterm.position<=io_helix[iocount].sh_cterm.position
         then begin
           io_assoc[iocount]:=oicount;
           oi_assoc[oicount]:=iocount;
           iocount:=iocount+1;
           oicount:=oicount+1;
         end else begin
           io_assoc[iocount]:=0;
           iocount:=iocount+1;
         end;
    end else
    begin
      if io_helix[iocount].sh_nterm.position<=oi_helix[oicount].sh_cterm.position
         then begin
           io_assoc[iocount]:=oicount;
           oi_assoc[oicount]:=iocount;
           iocount:=iocount+1;
           oicount:=oicount+1;
         end else begin
           oi_assoc[oicount]:=0;
           oicount:=oicount+1;
         end;
    end;
  until (iocount>io_count) and (oicount>oi_count);
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
                       var model_ct : integer;       {number of segments in model}
                       var plus_ct  : integer;       {number of correct significants}
                       var minus_ct : integer);      {number of false significants}
var
    io_ct,oi_ct  : integer;
    nt_pos       : integer;
    done         : boolean;
begin
  io_ct:=1; oi_ct:=1;
  nt_pos:=1;
  done:=false;
  model_ct:=0;
  plus_ct:=0;
  minus_ct:=0;

  repeat
    if ori_io then begin
      part5_findnext(io_helix,io_ct,io_count,nt_pos,thres3);
      if io_ct=0
      then done:=true
      else begin
        model_ct:=model_ct+1;
        model[model_ct]:=io_helix[io_ct];
        if io_assoc[io_ct]<>0 then begin
          if (oi_helix[io_assoc[io_ct]].score-io_helix[io_ct].score) >thres2
              then minus_ct:=minus_ct+1
          else if (io_helix[io_ct].score-oi_helix[io_assoc[io_ct]].score) > thres2
              then plus_ct:=plus_ct+1;
        end;
        nt_pos:=model[model_ct].sh_cterm.position+1;
      end;
    end else begin
      part5_findnext(oi_helix,oi_ct,oi_count,nt_pos,thres3);
      if oi_ct=0
      then done:=true
      else begin
        model_ct:=model_ct+1;
        model[model_ct]:=oi_helix[oi_ct];
        if oi_assoc[oi_ct]<>0 then begin
          if (io_helix[oi_assoc[oi_ct]].score-oi_helix[oi_ct].score) >thres2
              then minus_ct:=minus_ct+1
          else if (oi_helix[oi_ct].score-io_helix[oi_assoc[oi_ct]].score) > thres2
              then plus_ct:=plus_ct+1;
        end;
        nt_pos:=model[model_ct].sh_cterm.position+1;
      end;
    end;
    ori_io:=not ori_io;
  until done;
end;

(*procedure out_part5_model(model : helixarray; count: integer; totscore : integer);
var i : integer;
begin
  writeln(outdev,count:2,' strong transmembrane helices, total score : ',totscore);
    writeln(outdev,'#':2,' ','from':4,' ','to':4,' ','length ','score':5,' ','orientation');
  for i:=1 to count do with model[i] do begin
    write(outdev,i:2,' ',nterm.position:4,' ',cterm.position:4,' ');
    write(outdev,'(',cterm.position-nterm.position+1:2,')   ',score:5,' ');
    if nt_in then writeln(outdev,'i-o') else writeln(outdev,'o-i');
  end;
  writeln;
end; *)

(*procedure output_part5;
var in_model,out_model         : helixarray;
    in_mod_count,out_mod_count : integer;
    in_mod_score,out_mod_score : integer;
    i                          : integer;
begin
  part5_model(true,in_model,in_mod_count);
  in_mod_score:=0;
  for i:=1 to in_mod_count do in_mod_score:=in_mod_score+in_model[i].score;
  part5_model(false,out_model,out_mod_count);
  out_mod_score:=0;
  for i:=1 to out_mod_count do out_mod_score:=out_mod_score+out_model[i].score;
  writeln(outdev,'Prediction for sequence : ',seqid,' length=',seqlen:4);
  writeln(outdev,'2 possible models considered, only significant TM-segments used');
  writeln;
  if in_mod_count<>out_mod_count then begin
    writeln(outdev,'*** the models differ in the number of TM-helices ! ***');
    writeln;
  end;
  if in_mod_score>out_mod_score then begin
    if in_mod_count=0 then begin
      writeln(outdev,'!!! probably no transmembrane protein - no possible model found !!!');
    end
    else begin
      if (in_mod_score-out_mod_score/in_mod_count)>thres4 
      then  writeln(outdev,'-----> STRONGLY prefered model: N-terminus inside')
      else  writeln(outdev,'-----> slightly prefered model: N-terminus inside');
      out_part5_model(in_model,in_mod_count,in_mod_score);
      writeln(outdev,'------> alternative model');
      out_part5_model(out_model,out_mod_count,out_mod_score);
    end;
  end else begin
    if out_mod_count=0 then begin
      writeln(outdev,'!!! probably no transmembrane protein - no possible model found !!!');
    end
    else begin
      if (out_mod_score-in_mod_score/out_mod_count)>thres4
      then  writeln(outdev,'-----> STRONGLY prefered model: N-terminus outside')
      else  writeln(outdev,'-----> slightly prefered model: N-terminus outside');
      out_part5_model(out_model,out_mod_count,out_mod_score);
      writeln(outdev,'------> alternative model');
      out_part5_model(in_model,in_mod_count,in_mod_score);
    end;
  end;
end;       *)

procedure batch_output_part5;
var in_model,out_model         : helixarray;
    in_mod_count,out_mod_count : integer;
    in_mod_score,out_mod_score : integer;
    i                          : integer;
    in_plus_ct,in_minus_ct,
    out_plus_ct,out_minus_ct   : integer;
begin
  part5_model(true,in_model,in_mod_count,in_plus_ct,in_minus_ct);
  in_mod_score:=0;
  for i:=1 to in_mod_count do in_mod_score:=in_mod_score+in_model[i].score;

  part5_model(false,out_model,out_mod_count,out_plus_ct,out_minus_ct);
  out_mod_score:=0;
  for i:=1 to out_mod_count do out_mod_score:=out_mod_score+out_model[i].score;

  if in_mod_score>out_mod_score then begin
    write(outdev1,'"',seqid,'",',seqlen:4,',',in_mod_count:2);
    if in_mod_count=0 then begin
      writeln(outdev1,',"","",0,0');
    end else begin
      if (in_mod_score-out_mod_score/in_mod_count)>thres4  
      then write(outdev1,',"i","s"') else write(outdev1,',"i","w"');
      writeln(outdev1,',',in_plus_ct:2,',',in_minus_ct:2);
      for i:=1 to in_mod_count do begin
        write(outdev2,'"',seqid,'",',i:2,',',in_model[i].score:4);
        write(outdev2,',',in_model[i].nterm.position:4,',',in_model[i].cterm.position:4);
        if in_model[i].nt_in then writeln(outdev2,',"io"')
                             else writeln(outdev2,',"oi"');
      end;
    end;
  end
  else begin
    write(outdev1,'"',seqid,'",',seqlen:4,',',out_mod_count:2);
    if out_mod_count=0 then begin
      writeln(outdev1,',"","",0,0');
    end else begin
      if (out_mod_score-in_mod_score/out_mod_count)>thres4 
      then write(outdev1,',"o","s"') else write(outdev1,',"o","w"');
      writeln(outdev1,',',out_plus_ct:2,',',out_minus_ct:2);
      for i:=1 to out_mod_count do begin
        write(outdev2,'"',seqid,'",',i:2,',',out_model[i].score:4);
        write(outdev2,',',out_model[i].nterm.position:4,',',out_model[i].cterm.position:4);
        if out_model[i].nt_in then writeln(outdev2,',"io"')
                              else writeln(outdev2,',"oi"');
      end;
    end;
  end;
end;

{------------------------------------------------------------------------------}

begin
{$ifdef Windows}
  screensize.y:=80;
{$endif}
  thres1:=80;
  thres2:=200;
  thres3:=1000;
  thres4:=80;
  writeln('Prediction of transmembrane helices by the matrix method');
  writeln('different probabilities for inside-to-outside and outside-');
  writeln('to-inside helices are calculated');
  writeln;
  tablename:='matrix.tab';
  if not fexist(tablename) then begin
    tablename:=get_logical('TMPREDDIR')+d_sep+'matrix.tab';
  end;
  while (not fexist(tablename)) and (errcount<3) do begin
    errcount:=errcount+1;
    write('matrix file:',tablename,' does not exist, enter valid name : ');
    readln(tablename);
  end;
  write('Name of Sequence (in arbitrary format) : ');
  readln(seqname);
  errcount:=0;
  while (not fexist(seqname)) and (errcount<3) do begin
    errcount:=errcount+1;
    write('file does not exist, try again : ');
    readln(seqname);
  end;
  write('minimal length of transmembrane sequence : ');
  readln(minwidth);
  write('maximal length of transmembrane sequence : ');
  readln(maxwidth);
  if getmatrices<>0 then writeln(chr(7),'error in the matrices');

  write('name of first output file (all entries) : ');
  readln(outname1);

{$ifdef Borland}  
  if outname1<>'' then begin
    fopen(outdev1,outname1,'w');
  end else begin
    assigncrt(outdev1);
    rewrite(outdev1);
  end;
{$endif}
{$ifdef dec}  
    fopen(outdev1,outname1,'w');
{$endif}

  write('name of second output file (all helices) : ');
  readln(outname2);
{$ifdef Borland}  
  if outname2<>'' then begin
    fopen(outdev2,outname2,'w');
  end else begin
    assigncrt(outdev2);
    rewrite(outdev2);
  end;
{$endif}
{$ifdef dec}  
    fopen(outdev2,outname2,'w');
{$endif}

  no_seq:=0;
  seqlen:=maxseqlen;
  seqerr:=sq_read_multifile(tseq,seqlen,seqid,seqname,no_seq,max_no_seq);
  for no_seq:=1 to max_no_seq do begin
    seqlen:=maxseqlen;
    seqerr:=sq_read_multifile(tseq,seqlen,seqid,seqname,no_seq,max_no_seq);
    if (seqerr<>0) then writeln(chr(7),'error in the sequence No:',no_seq,' name : ',seqid)
    else begin
      for i:=1 to length(seqid) do seqid[i]:=upcase(seqid[i]);
      writeln(seqid);
      do_part1;                                        {the auxiliary table}
      { output_part1;  }
      do_part2;                                        {the predicition graphics}
      {output_part2; }
      do_part3;                                        {the actual prediction}
      {output_part3;} 
      do_part4;
      batch_output_part5;
    end;
  end;
  close(outdev1);
  close(outdev2);
{$ifdef Windows}
  donewincrt;
{$endif}
end.

