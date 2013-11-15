#!/usr/freeware/bin/gawk -f

BEGIN {

  # do average of entries in files
  filecnt=ARGC;
  for (i=1;i<filecnt;i++) {
      while (getline < ARGV[i] > 0) {
            idset[$2]=1;
            score[$2,i]=$1;
      }
      close(ARGV[i]);
      delete ARGV[i];
  }

  for (id in idset) {
      ok=1;
      sum=0;
      cnt=0;
      for (i=1;i<filecnt;i++) {
          if (score[id,i]=="") ok=0;
          sum+=score[id,i];
          cnt++;
      }
      if (ok && cnt>0) {
         printf("%.6f %s\n",sum/cnt,id);
      }
  }
}

