#! /usr/freeware/bin/gawk -f

BEGIN {

  # passed from the outside: tmpdir is $PFUNTMP, datadir $PFUNDAT

  k=1;
  fn=datadir"/secretome.prob"
  while ( getline < fn > 0 )
     seca[k++]=$1;

  fname=tmpdir"/secretome.pscore";
  while (getline < fname > 0) {
     idset[$2]=1;
     secretomepscore[$2,1]=$1;
     isdefsecretome[$2,1]=1;
  }

  printf("############# SecretomeP 1.0 predictions #############\n\n");

  for (id in idset) {
  printf(">%s\n\n",id);

  printf("\n# Non-classical secretion		   Prob     Odds\n");
  issec=-1;
  if (isdefsecretome[id,1] && secretomepscore[id,1]>0) {
     if (secretomepscore[id,1]>2*seca[1])
  	issec=1;
     else
  	issec=0;
  }
  else {
     issec=-1; # fail
  }

  if ( isdefsecretome[id,1] ) {
     printf("  Secreted 			%s%5.3f %8.3f\n",
  	    (issec==1)?" => ":"    ",
  	    secretomepscore[id,1],
  	    ( secretomepscore[id,1]/seca[1] ));
  }
  else {
     printf("  Secreted 			%s%s\n",
  	    (issec==1)?" => ":"    ",
  	    "Failed");
  }

  printf("\n//\n\n");
  }

}
