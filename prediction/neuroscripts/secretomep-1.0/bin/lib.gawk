#!/usr/freeware/bin/gawk -f

function fileexists( testfile,   l) {
   l=(getline xx<testfile);
   if (l!=-1) return 1;
   return 0;
}

# This functions make a sums vector of the values in bins defined by
#  the parameters;
# 
#  bsz;       Bin size (assumed odd, or increased by 1!)
#  ndb;       Number of bins to be distributed (minimum 3)
#  len;       Length of sequence
#
#  setXval;   Values of sites.
#  setXpos;   Positions of sites.
#  cntX;      Number of sites.
#
#  setXvalA;  Resulting sum vector. [1..ndb]
#
function get_dist_bin(bsz,ndb,len,setXval,setXpos,cntX,setXvalA,      i,j,p1,p2,c,l,d) {

   if (bsz%2==0) bsz=bsz+1;

   c=int(bsz/2.0)+1;

#   print "bsz:"bsz" ndb:"ndb" len:"len;

   p1[1]=1;
   p2[1]=bsz;
   setXvalA[1]=0;

   p1[ndb]=len-bsz+1;
   p2[ndb]=len;
   setXvalA[ndb]=0;

   l=(len-2*c);
   d=(p1[ndb]-p2[1]+2*c-2)/(ndb-1);
   if (d<0) d=0;
#   print "l: "l" d: "d" c:"c;
   
#   print "1: p1: "p1[1]"   p2:"p2[1];

   for (i=2;i<ndb;i++) {
      setXvalA[i]=0;
      pc=int(c+(d*(i-1)));
      p1[i]=pc-c+1;
      p2[i]=pc+c-1;
#      print i": p1: "p1[i]"   p2:"p2[i];
   }

#   print ndb": p1: "p1[ndb]"   p2:"p2[ndb];

   for (i=1;i<cntX;i++) {
      for (j=1;j<=ndb;j++) {
	 if (setXpos[i]>=p1[j] && setXpos[i]<=p2[j]) {
	    setXvalA[j]+=setXval[i];
	 }
      }
   }
}

# This functions make a sums vector of the values in bins defined by
#  the parameters;
#
#  sbb, nsbb; Static bins beginning, and number of sbb. [1..nsbb]
#             a list of positions that defines the bins relative to the
#             beginning of the sequence.
#  sbe,nsbe ; Static bins ending, and number of sbe. [1..nsbe]
#             a list of positions that defines the bins relative to the
#             end of the sequence.
#  ndb;       Number of dynamic bins not covered by the beginning and end.
#
#  len;       Length of sequence.
#  setXval;   Values of sites.
#  setXpos;   Positions of sites.
#  cntX;      Number of sites.
#
#  setXvalA;  Resulting sum vector. [1..nsbb+nsbe+ndb]
#
function get_n_bin(sbb,nsbb,sbe,nsbe,ndb,len,setXval,setXpos,cntX,setXvalA,   nsbei,ndbi,i,j,p1,p2) {

   nsbei=nsbb+ndb;
   ndbi=nsbb;
   sbe[nsbe+1]=0;
   nsbe++;

   for (i=1;i<=nsbb+nsbe+ndb;i++)
     setXvalA[i]=0;
   
   for (i=1;i<cntX;i++) {

      p1=0;
      # Check begin bins
      for (j=1;j<=nsbb;j++) {
	 p2=sbb[j];
	 if (setXpos[i]>p1 && setXpos[i]<=p2) {
	    setXvalA[j]+=setXval[i];
#	    print "p: "setXpos[i]" bb:"j" v:"setXval[i];
	 }
	 p1=p2;
      }

      # Check end bins
      p1=len-sbe[1];
      for (j=2;j<=nsbe;j++) {	 
	 p2=len-sbe[j];
	 if (setXpos[i]>p1 && setXpos[i]<=p2) {
	    setXvalA[j+nsbei-1]+=setXval[i];
#	    print "p: "setXpos[i]" eb:"j-1" v:"setXval[i];
	 }
	 p1=p2;
      }

      if (nsbb==0)
	p1=0;
      else
	p1=sbb[nsbb];
      if (nsbe==0)
	p2=len;
      else
	p2=len-sbe[1];
      # Check dynamic middle bins
      if (setXpos[i]>p1 && setXpos[i]<=p2) {

         # Normalize the dynamic values by dividing with binsize
	 b=int((ndb*(setXpos[i]-p1-1))/(p2-p1));
	 setXvalA[ndbi+b+1]+= (setXval[i]/ ((p2-p1)/ndb));
#	 print "p: "setXpos[i]" mb:"b+1" v:"setXval[i];
      }
   }  
   
   
   return;
}

function showarray(numbin,setXmval,cntXm,dobool,fname,   i,fn2) {

   newfilename=fname"/type";
   l=(getline xx<newfilename);

   if (l!=-1) {
      fn2=fname"/feature.out";
      
      for (i=1;i<=numbin;i++) {
	 printf "%.3f ",setXmval[i] >> fn2 ;
      }
      # 0 or 1 -> assign from cntXm
      if (dobool==1) {
	 if (cntXm>1)
	   printf "1 " >> fn2;
	 else
	   printf "0 " >> fn2;
      }
      # 2 or 3 -> assign from dobool (0,1)
      if (dobool==2) 
	printf "0 " >> fn2;
      if (dobool==3)
	printf "1 " >> fn2;
   }
}

function endarray(fname,           fn2) {

   # make a newline on all arrays   
   if (fileexists(fname"/type")) {
      fn2=fname"/feature.out";
      printf "\n" >> fn2;
   }
   if (fileexists(fname"_5d/type")) {
      fn2=fname"_5d/feature.out";
      printf "\n" >> fn2;
   }
   if(fileexists(fname"_10d/type")) {
      fn2=fname"_10d/feature.out";
      printf "\n" >> fn2;
   }
   if(fileexists(fname"_101w5/type")) {
      fn2=fname"_101w5/feature.out";
      printf "\n" >> fn2;
   }
   if(fileexists(fname"_51w10/type")) {
      fn2=fname"_51w10/feature.out";
      printf "\n" >> fn2;
   }
   if(fileexists(fname"_5w100/type")) {
      fn2=fname"_5w100/feature.out";
      printf "\n" >> fn2;
   }
}

function beginarray(ident,fname,  fn2) {
   # make a newline on all arrays   
   if (fileexists(fname"/type")) {
      fn2=fname"/feature.out";
      printf "%s ",ident >> fn2;
   }
   if(fileexists(fname"_5d/type")) {
      fn2=fname"_5d/feature.out";
      printf "%s ",ident >> fn2;
   }
   if(fileexists(fname"_10d/type")) {
      fn2=fname"_10d/feature.out";
      printf "%s ",ident >> fn2;
   }
   if(fileexists(fname"_101w5/type")) {
      fn2=fname"_101w5/feature.out";
      printf "%s ",ident >> fn2;
   }
   if (fileexists(fname"_51w10/type")) {
      fn2=fname"_51w10/feature.out";
      printf "%s ",ident >> fn2;
   }
   if (fileexists(fname"_5w100/type")) {
      fn2=fname"_5w100/feature.out";
      printf "%s ",ident >> fn2;
   }
}

function clearit(arr) {
   foreach (i in arr)
     delete arr[i];
}

function makearrays(len,setXmval,setXmpos,cntXm,dobool,fname,  fn,startbin,numstartbin,endbin,numendbin,numdbin,numbin) {

   # Example feature vector; 

   if (fileexists(fname"/type")) {

      startbin[1]=25;
      numstartbin=1;
      
      endbin[1]=25;
      numendbin=1;
      
      numdbin=3;
      numbin= numstartbin+numendbin+numdbin;

      get_n_bin(startbin,numstartbin,endbin,numendbin,numdbin,len,setXmval,setXmpos,cntXm,setXmvalA1);
      fn=fname;
      showarray(numbin,setXmvalA1,cntXm,dobool,fn);
   }

#   clearit(setXmvalA1);

   # Example feature vector  5 dynamic; 

   if (fileexists(fname"_5d/type")) {
      numstartbin=0;
      numendbin=0;
      
      numdbin=5;
      numbin= numstartbin+numendbin+numdbin;
      
      fn=fname"_5d";
      get_n_bin(startbin,numstartbin,endbin,numendbin,numdbin,len,setXmval,setXmpos,cntXm,setXmvalA5d);
      showarray(numbin,setXmvalA5d,cntXm,dobool,fn);
   }
#   clearit(setXmvalA5d);
   
   # Example feature vector  10 dynamic; 
   if (fileexists(fname"_10d/type")) {
      numstartbin=0;
      numendbin=0;
      
      numdbin=10;
      numbin= numstartbin+numendbin+numdbin;
      fn=fname"_10d";
      get_n_bin(startbin,numstartbin,endbin,numendbin,numdbin,len,setXmval,setXmpos,cntXm,setXmvalA10d);
      showarray(numbin,setXmvalA10d,cntXm,dobool,fn);
   }
#   clearit(setXmvalA10d);

   # Example feature vector  101w 5 dynamic; 
   if (fileexists(fname"_101w5/type")) {
      binsize=101;
      numdbin=5;
      numbin= numstartbin+numendbin+numdbin;
      fn=fname"_101w5";
      get_dist_bin(binsize,numdbin,len,setXmval,setXmpos,cntXm,setXmvalA101w5);
      showarray(numbin,setXmvalA101w5,cntXm,dobool,fn);
   }
#   clearit(setXmvalA101w5);

   # Example feature vector  51w 10 dynamic; 
   if (fileexists(fname"_51w10/type")) {
      binsize=51;
      numdbin=10;
      numbin= numstartbin+numendbin+numdbin;
      fn=fname"_51w10";
      get_dist_bin(binsize,numdbin,len,setXmval,setXmpos,cntXm,setXmvalA51w10);
      showarray(numbin,setXmvalA51w10,cntXm,dobool,fn);
   }
#   clearit(setXmvalA51w10);

   # Example feature vector  5w 100 dynamic; 
   if (fileexists(fname"_5w100/type")) {
      binsize=5;
      numdbin=100;
      numbin= numstartbin+numendbin+numdbin;
      fn=fname"_5w100";
      get_dist_bin(binsize,numdbin,len,setXmval,setXmpos,cntXm,setXmvalA5w100);
      showarray(numbin,setXmvalA5w100,cntXm,dobool,fn);
   }
#   clearit(setXmvalA5w100);
#   clearit(setXmval);
}

