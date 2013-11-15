#!/usr/freeware/bin/gawk -f
##
##
## merge.gawk <template> <feature> <feature> .. <feature>
##
## feature; name of feature directory.
##
##

BEGIN{
   
   template=ARGV[1];
   delete ARGV[1];
   datacnt=0;
   for (i=2;i<ARGC;i++) {
      feature[i]=ARGV[i];
      filename[i]=template".Features/"ARGV[i]"/feature.norm.out";
      
      delete ARGV[i];
 
      while (getline < filename[i] > 0 ) {
#	 printf "." >> "/dev/stderr";	
	 data[i,$1]=substr($0,length($1)+2,length($0));
	 idset[$1]=1;
	 datacnt++;
      }

#      printf "\n" >> "/dev/stderr";	
      close(filename[i]);
   }

   if (datacnt==0) {
      printf "Error: all ids are missing all features\n" >> "/dev/stderr";
      exit;
   }


   count=0;
   for (id in idset) {
# REMOVE THE FOLLOWING LINE IF YOU DON'T WANT THE ID ALSO!
#      printf("%-11s ",id);
#      printf "." >> "/dev/stderr";
      assignment="";
      ok=1;
      for (f=2;f<ARGC;f++) {
	 if (data[f,id]=="") 
	   ok=0;
	 assignment=assignment" "data[f,id];
      }

      if (ok) {
         printf("%s 1 %s\n",assignment,id);
         count++;
      }else{
         printf "Error: id "id" is missing features ( " >> "/dev/stderr";
	 for (f=2;f<ARGC;f++) {
	    if (data[f,id]=="") { 
	       printf feature[f]" " >> "/dev/stderr";
	    }
	 }
	 printf ")\n" >> "/dev/stderr"; 
      }
   }
#   print "## Written to howlin file set "setsel" : "count >> "/dev/stderr";
}




