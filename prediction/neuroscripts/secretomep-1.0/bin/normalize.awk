#!/usr/freeware/bin/gawk -f
BEGIN{

   fname=feature"/feature.stat";
   getline < fname;
   getline < fname;
   col=2;
   while ( getline < fname > 0 ) {
      ident[col]=$2;
      avg[col]=0+$3;
      dev[col]=0+$4;
#    print "avg:"avg[col]"  dev:"dev[col] >> "/dev/stderr";
      if (dev[col]+0<=0)
	dev[col]=1;
      norm[col]=$7;
      col++;
   }
   close(fname);
   
   w=0;
   fname=feature"/feature.out";
   while ( getline < fname > 0 ) {
      printf $1" ";
      for (i=2;i<=NF;i++) {
	 if (col<=i && w<4) {
	    print "Warning, encountered column "i" ("$1")" >> "/dev/stderr";
	    w++;
	 }else{
	    if ($i=="x") {
	       $i=0;
	    }else if (norm[i]=="log") {
	       if ($i<0 && w[i]<4) {
		  printf("warning: norm type 'log' with value %s (%s,%s) \n",$i,ident[i],$1) >> "/dev/stderr";
		  w[i]++;
	       }
	       if ($i<0) $i=0;

	       $i=log(1+$i);
	    } else if (norm[i]=="bool") {
	       if ($i!=1 && $i!=0 && w[i]<4) {
		  printf("warning: norm type 'bool' with value %s (%s,%s)\n",$i,ident[i],$1) >> "/dev/stderr";
		  w[i]++;
	       }
	    }else if (norm[i]=="pct") {
	       if (($i<0 || $i>100) && w[i]<4) {
		  printf("warning: norm type 'pct' with value %s (%s,%s)\n",$i,ident[i],$1) >> "/dev/stderr";
		  w[i]++;
	       }
	    }else if (norm[i]=="freq") {
	       if (($i<0 || $i>1) && w[i]<4) {
		  printf("warning: norm type 'freq' with value %s (%s,%s)\n",$i,ident[i],$1) >> "/dev/stderr";
		  w[i]++;
	       }
	    }else if (norm[i]=="0norm") {
	       if ($i<0 && w[i]<4) {
		  printf("warning: norm type '0norm' with value %s (%s,%s)\n",$i,ident[i],$1) >> "/dev/stderr";
		  w[i]++;
	       }
	    }
	    
	    printf "%.3f ",($i-avg[i])/dev[i];
	 }
      }
      printf "\n";
      
   }
   close(fname);
}
