
/*****************************************************************************/
/***  (seg.c)                                                              ***/
/*** #include precomputed ln(fact(n)) array in "lnfac.h"                   ***/
/*****************************************************************************/

/*--------------------------------------------------------------(includes)---*/

#include "genwin.h"
#include "lnfac.h"

/*---------------------------------------------------------------(defines)---*/

#define LENCORRLIM 120

/*---------------------------------------------------------------(structs)---*/

struct Segment
  {
   int begin;
   int end;
   struct Segment *next;
  };

/*---------------------------------------------------------------(globals)---*/

int window = 12;
int downset, upset;
double locut = 2.2;
double hicut = 2.5;

int hilenmin = 0;
int overlaps = FALSE;
int hionly = FALSE;
int loonly = FALSE;
int entinfo = TRUE;
int singleseq = FALSE;
int prettyseq = FALSE;
int prettytree = TRUE;
int charline = 60;
int maxtrim = 100;

double getprob(), lnperm(), lnass();

/*------------------------------------------------------------------(main)---*/

main(argc, argv)
  int argc;
  char *argv[];

  {struct Database *db;
   struct Sequence *seq;
   struct Segment *segs;
   int ctime;

   genwininit();
/* readlencorr(); */                        /* #include lencorr file */
   getparams(argc, argv);

   if ((db=opendbase(argv[1]))==NULL)
     {
      fprintf(stderr, "Error opening file %s\n", argv[1]);
      exit(1);
     }

   for (seq=firstseq(db); seq!=NULL; seq=nextseq(db))
     {
      segs = (struct Segment *) NULL;
      segseq(seq, &segs, 0);
      mergesegs(seq, segs);

      if (singleseq) singreport(seq, segs);
      else if (prettyseq) prettyreport(seq, segs);
      else if (prettytree) pretreereport(seq, segs);
      else report(seq, segs);

      freesegs(segs);
      closeseq(seq);
     }

   closedbase(db);
   exit(0);
  }

/*-------------------------------------------------------------(getparams)---*/

getparams(argc, argv)
  int argc;
  char *argv[];

  {int i;
   int nargc;
   char **nargv;
   extern char *optarg;
   extern int optind;
   int c;

   if (argc<2)
     {
      usage();
      exit(1);
     }

   for (i=2; argc>i && argv[i][0]!='-'; i++)
     {
      if (i==2)
        {
         window = atoi(argv[2]);
        }
      else if (i==3)
        {
         locut = atof(argv[3]);
        }
      else if (i==4)
        {
         hicut = atof(argv[4]);
        }
     }

   if (locut>hicut)
     {
      fprintf(stderr, "Warning: trigger entropy greater than extension entropy\n");
      hicut = locut;
     }

   downset = (window+1)/2 - 1;
   upset = window - downset;

   if (i==argc) return;

   nargc = argc-i+1;
   nargv = argv+(i-1);
   while ((c=getopt(nargc, nargv, "m:olhaxpqc:nt:"))!=-1)
     {
      switch(c)
        {
         case 'm':
            hilenmin = atoi(optarg);
            break;
         case 'o':
            overlaps = TRUE;
            hilenmin = 0;
            break;
         case 'l':
            loonly = TRUE;
            singleseq = FALSE;
            prettyseq = FALSE;
            prettytree = FALSE;
            break;
         case 'h':
            hionly = TRUE;
            singleseq = FALSE;
            prettyseq = FALSE;
            prettytree = FALSE;
            break;
         case 'a':
            hionly = FALSE;
            loonly = FALSE;
            singleseq = FALSE;
            prettyseq = FALSE;
            prettytree = FALSE;
            break;
         case 'x':
            singleseq = TRUE;
            prettyseq = FALSE;
            prettytree = FALSE;
            hionly = TRUE;
            loonly = FALSE;
            break;
         case 'p':
            prettytree = TRUE;
            prettyseq = FALSE;
            singleseq = FALSE;
            hionly = FALSE;
            loonly = FALSE;
            break;
         case 'q':
            prettyseq = TRUE;
            prettytree = FALSE;
            singleseq = FALSE;
            hionly = FALSE;
            loonly = FALSE;
            break;
         case 'c':
            charline = atoi(optarg);
            break;
         case 'n':
            entinfo = FALSE;
            break;
         case 't':
            maxtrim = atoi(optarg);
            break;
         case '?':
            fprintf(stderr, "Unknown option.\n");
            usage();
            exit(1);
            break;
        }
     }   

   return;
  }

/*---------------------------------------------------------------(segment)---*/

segseq(seq, segs, offset)
  struct Sequence *seq;
  struct Segment **segs;
  int offset;

  {struct Segment *seg, *leftsegs;
   struct Sequence *leftseq;
   int first, last, lowlim;
   int loi, hii, i;
   int leftend, rightend, lend, rend;
   double *H, *seqent();

   H = seqent(seq);
   if (H==NULL) return;

   first = downset;
   last = seq->length - upset;
   lowlim = first;

   for (i=first; i<=last; i++)
     {
      if (H[i]<=locut && H[i]!=-1)
        {
         loi = findlo(i, lowlim, H);
         hii = findhi(i, last, H);

         leftend = loi - downset;
         rightend = hii + upset - 1;

         trim(openwin(seq, leftend, rightend-leftend+1), &leftend, &rightend);

         if (i+upset-1<leftend)   /* check for trigger window in left trim */
           {
            lend = loi - downset;
            rend = leftend - 1;

            leftseq = openwin(seq, lend, rend-lend+1);
            leftsegs = (struct Segment *) NULL;
            segseq(leftseq, &leftsegs, offset+lend);
            if (leftsegs!=NULL)
              {
               if (*segs==NULL) *segs = leftsegs;
               else appendseg(*segs, leftsegs);
              }
            closewin(leftseq);

/*          trim(openwin(seq, lend, rend-lend+1), &lend, &rend);
            seg = (struct Segment *) malloc(sizeof(struct Segment));
            seg->begin = lend;
            seg->end = rend;
            seg->next = (struct Segment *) NULL;
            if (segs==NULL) segs = seg;
            else appendseg(segs, seg);  */
           }

         seg = (struct Segment *) malloc(sizeof(struct Segment));
         seg->begin = leftend + offset;
         seg->end = rightend + offset;
         seg->next = (struct Segment *) NULL;

         if (*segs==NULL) *segs = seg;
         else appendseg(*segs, seg);

         i = min(hii, rightend+downset);
         lowlim = i + 1;
/*       i = hii;     this ignores the trimmed residues... */
        }
     }

   free(H);
   return;
  }

/*----------------------------------------------------------------(seqent)---*/

double *seqent(seq)
  struct Sequence *seq;

  {struct Sequence *win;
   double *H;
   int i, first, last;

   if (window>seq->length)
     {
      return((double *) NULL);
     }

   H = (double *) malloc(seq->length*sizeof(double));

   for (i=0; i<seq->length; i++)
     {
      H[i] = -1.;
     }

   win = openwin(seq, 0, window);
   enton(win);

   first = downset;
   last = seq->length - upset;

   for (i=first; i<=last; i++)
     {
      if (seq->punctuation && hasdash(win))
        {H[i] = -1;
         shiftwin(win, 1);
         continue;}
      H[i] = win->entropy;
      shiftwin(win, 1);
     }

   closewin(win);
   return(H);
  }

/*---------------------------------------------------------------(hasdash)---*/

hasdash(win)
  struct Sequence *win;

  {int i, len;

   len = win->length;

   for (i=0; i<len; i++)
     {
      if (win->seq[i]=='-')
        {
         return(TRUE);
        }
     }

   return(FALSE);
  }

/*----------------------------------------------------------------(findlo)---*/

findlo(i, limit, H)
  int i, limit;
  double *H;

  {int j;

   for (j=i; j>=limit; j--)
     {
      if (H[j]==-1) break;
      if (H[j]>hicut) break;
     }

   return(j+1);
  }

/*----------------------------------------------------------------(findhi)---*/

findhi(i, limit, H)
  int i, limit;
  double *H;

  {int j;

   for (j=i; j<=limit; j++)
     {
      if (H[j]==-1) break;
      if (H[j]>hicut) break;
     }

   return(j-1);
  }

/*------------------------------------------------------------------(trim)---*/

trim(seq, leftend, rightend)
  struct Sequence *seq;
  int *leftend, *rightend;

  {struct Sequence *win;
   double prob, minprob;
   int shift, len, i;
   int lend, rend;
   int minlen;

/* fprintf(stderr, "%d %d\n", *leftend, *rightend);  */

   lend = 0;
   rend = seq->length - 1;
   minlen = 1;
   if ((seq->length-maxtrim)>minlen) minlen = seq->length-maxtrim;

   minprob = 1.;
   for (len=seq->length; len>minlen; len--)
     {
      win = openwin(seq, 0, len);
      stateon(win);
      i = 0;

      shift = TRUE;
      while (shift)
        {
         prob = getprob(win->state, len);
         if (prob<minprob)
           {
            minprob = prob;
            lend = i;
            rend = len + i - 1;
           }
         shift = shiftwin(win, 1);
         i++;
        }
      closewin(win);
     }

/* fprintf(stderr, "%d-%d ", *leftend, *rightend);  */

   *leftend = *leftend + lend;
   *rightend = *rightend - (seq->length - rend - 1);

/* fprintf(stderr, "%d-%d\n", *leftend, *rightend);  */

   closewin(seq);
   return;
  }

/*---------------------------------------------------------------(getprob)---*/

double getprob(sv, total)
  int *sv;
  int total;

  {double ans, totseq;

   totseq = ((double) total) * log((double) 20.);

   ans = lnass(sv) + lnperm(sv, total) - totseq;

   return(ans);
  }

/*----------------------------------------------------------------(lnperm)---*/

double lnperm(sv, tot)
  int *sv;
   int tot;

  {double ans;
   int i;

   ans = lnfac[tot];

   for (i=0; sv[i]!=0; i++) 
     {
      ans -= lnfac[sv[i]];
     }

   return(ans);
  }

/*-----------------------------------------------------------------(lnass)---*/

double lnass(sv)
  int *sv;

  {double ans;
   int class, total, i;

   ans = lnfac[20];
   if (sv[0]==0) return(ans);
   if (sv[19]==1) return(0.);

   total = 20;
   class = 1;
   for (i=1;; i++)
     {
      if (sv[i]==sv[i-1])
        {
         class++;
         continue;
        }
      else
        {
         total -= class;
         ans -= lnfac[class];
         if (sv[i]==0)
           {
            ans -= lnfac[total];
            break;
           }
         else
           {
            class = 1;
            continue;
           }
        }
     }

   return(ans);
  }

/*-------------------------------------------------------------(mergesegs)---*/

mergesegs(seq, segs)
  struct Sequence *seq;
  struct Segment *segs;

  {struct Segment *seg, *nextseg;
   int len;

   if (overlaps) return;
   if (segs==NULL) return;

   if (segs->begin<hilenmin) segs->begin = 0;

   seg = segs;
   nextseg = seg->next;

   while (nextseg!=NULL)
     {
      if (seg->end>=nextseg->begin)               /* overlapping segments */
        {
         seg->end = nextseg->end;
         seg->next = nextseg->next;
         free(nextseg);
         nextseg = seg->next;
         continue;
        }
      len = nextseg->begin - seg->end - 1;
      if (len<hilenmin)                            /* short hient segment */
        {
         seg->end = nextseg->end;
         seg->next = nextseg->next;
         free(nextseg);
         nextseg = seg->next;
         continue;
        }
      seg = nextseg;
      nextseg = seg->next;
     }

   len = seq->length - seg->end - 1;
   if (len<hilenmin) seg->end = seq->length - 1;

   return;
  }

/*----------------------------------------------------------------(report)---*/

report(seq, segs)
  struct Sequence *seq;
  struct Segment *segs;

  {struct Sequence *subseq;
   struct Segment *seg, *nextseg;
   static int hi = 1;
   static int lo = 0;

   if (segs==NULL)
     {
      enton(seq);
      seqout(seq, hi, 1, seq->length);
/*    fputc('\n', stdout);   -- for spacing after each sequence */
      return;
     }

   if (segs->begin>0)
     {
      subseq = openwin(seq, 0, segs->begin);
      enton(subseq);
      seqout(subseq, hi, 1, segs->begin);
      closewin(subseq);
     }

   for (seg=segs; seg!=NULL; seg=seg->next)
     {
      subseq = openwin(seq, seg->begin, seg->end-seg->begin+1);
      enton(subseq);
      seqout(subseq, lo, seg->begin+1, seg->end+1);
      closewin(subseq);

      if (seg->next==NULL)
        {
         break;
        }

      nextseg = seg->next;
      
      if (nextseg->begin<=seg->end)
        {
         fprintf(stderr, "Overlapping segments: %s\n", seq->id);
         continue;
        }

      if (nextseg->begin==seg->end+1)
        {
         continue;
        }

      subseq = openwin(seq, seg->end+1, nextseg->begin-seg->end-1);
      enton(subseq);
      seqout(subseq, hi, seg->end+2, nextseg->begin);
      closeseq(subseq);
     }

   if (seg->end+1==seq->length)
     {
/*    fputc('\n', stdout);   -- for spacing after each sequence */
      return;
     }

   subseq = openwin(seq, seg->end+1, seq->length-seg->end-1);
   enton(subseq);
   seqout(subseq, hi, seg->end+2, seq->length);
   closeseq(subseq);

/* fputc('\n', stdout);   -- for spacing after each sequence */
   return;
  }

/*------------------------------------------------------------(singreport)---*/

singreport(seq, segs)
  struct Sequence *seq;
  struct Segment *segs;

  {char *proseq;
   struct Segment *seg;
   int begin, end, i, ctr;

   proseq = seq->seq;
   upper(proseq);

   for (seg=segs; seg!=NULL; seg=seg->next)
     {
      begin = seg->begin;
      end = seg->end;

      for (i=begin; i<=end; i++)
        {
         proseq[i] = 'x';
        }
     }

   fprintf(stdout, "%s\n", seq->header);

   for (i=0, ctr=0; proseq[i]!='\0'; i++, ctr++)
     {
      if (ctr==charline)
        {
         fputc('\n', stdout);
         ctr = 0;
        }
      fputc(proseq[i], stdout);
     }

   fputc('\n', stdout);
   fputc('\n', stdout);
  }

/*----------------------------------------------------------(prettyreport)---*/

prettyreport(seq, segs)
  struct Sequence *seq;
  struct Segment *segs;

  {char *proseq;
   char format[10];
   struct Segment *seg;
   int begin, end, i, ctr;
   int leftspace;

   leftspace = (int) ceil(log10((double)seq->length));
   sprintf(format, "%%%dd ", leftspace);

   proseq = seq->seq;
   upper(proseq);

   for (seg=segs; seg!=NULL; seg=seg->next)
     {
      begin = seg->begin;
      end = seg->end;

      for (i=begin; i<=end; i++)
        {
         if (isalpha(proseq[i])) proseq[i] = tolower(proseq[i]);
        }
     }

   fprintf(stdout, "%s\n", seq->header);

   for (i=0; i<=leftspace; i++)
     {
      fputc(' ', stdout);
     }
   for (i=0, ctr=1; i<charline; i++, ctr++)
     {
      if (ctr==10)
        {
         fputc('.', stdout);
         ctr = 0;
        }
      else
        {
         fputc(' ', stdout);
        }
     }
   fputc('\n', stdout);
   fprintf(stdout, format, 1);

   for (i=0, ctr=0; proseq[i]!='\0'; i++, ctr++)
     {
      if (ctr==charline)
        {
         fputc('\n', stdout);
         fprintf(stdout, format, i+1);
         ctr = 0;
        }
      fputc(proseq[i], stdout);
     }

   fprintf(stdout, " %d\n", seq->length);
   fputc('\n', stdout);
  }

/*---------------------------------------------------------(pretreereport)---*/

pretreereport(seq, segs)
  struct Sequence *seq;
  struct Segment *segs;

  {struct Sequence *win;
   struct Segment *seg;
   char buffer[100], leftfmt[20], rightfmt[20];
   char *curseq;
   int left, right, len;
   int current, nextloent;
   int cline;

   cline = charline / 2;
   
   fprintf(stdout, "%s\n\n", seq->header);
   sprintf(leftfmt, "%%%ds", cline);
   sprintf(rightfmt, "%%-%ds", cline);

   current = 0;

   for (seg=segs; ; seg=seg->next)
     {
      if (seg==NULL) nextloent = seq->length;
      else nextloent = seg->begin;

      if (current < nextloent)
        {
         left = current;
         right = nextloent - 1;
         len = right - left + 1;
         win = openwin(seq, left, len);
         curseq = win->seq;
         upper(curseq);

         space(cline);
         fprintf(stdout, " %4d-%-4d ", left+1, right+1);
         
         while (len>0)
           {
            if (len<=cline)
              {
               fprintf(stdout, rightfmt, curseq);
               fprintf(stdout, "\n");
               break;
              }
            else
              {
               strncpy(buffer, curseq, cline);
               buffer[cline] = '\0';
               fprintf(stdout, rightfmt, buffer);
               fprintf(stdout, "\n");
               space(cline+11);
               curseq += cline;
               len -= cline;
              }
           }

         closewin(win);
        }

      if (seg==NULL) break;

      left = seg->begin;
      right = seg->end;
      len = right - left + 1;
      win = openwin(seq, left, len);
      curseq = win->seq;
      lower(curseq);

      strncpy(buffer, curseq, cline);
      buffer[cline] = '\0';
      fprintf(stdout, leftfmt, buffer);
      fprintf(stdout, " %4d-%-4d ", left+1, right+1);
      fprintf(stdout, "\n");

      len -= cline;
      if (len>0) curseq += cline;

      while (len>0)
        {
         strncpy(buffer, curseq, cline);
         buffer[cline] = '\0';
         fprintf(stdout, leftfmt, buffer);
         fprintf(stdout, "\n");
         len -= cline;
         if (len>0) curseq += cline;
        }

      closewin(win);
      current = right + 1;
     }

   fprintf(stdout, "\n");
  }

/*-----------------------------------------------------------------(space)---*/

space(len)
  int len;

  {int i;

   for (i=0; i<len; i++) fputc(' ', stdout);
  }

/*----------------------------------------------------------------(seqout)---*/

seqout(seq, hilo, begin, end)
  struct Sequence *seq;
  int hilo;
  int begin, end;

#define HDRLEN 60

  {char *proseq, *id, *header;
   char outbuf[HDRLEN+1];
   static int hi = 1;
   static int lo = 0;
   int i, ctr, iend;

   if (hionly && hilo==lo) return;
   if (loonly && hilo==hi) return;

   proseq = seq->seq;
   id = seq->id;
   if (id==NULL) id = seq->parent->id;
   header = seq->header;
   if (header==NULL) header = seq->parent->header;

   iend = findchar(header, ' ');
   if (iend!=-1) header = header+iend;

   if (entinfo)
     {
      fprintf(stdout, ">%s(%d-%d)", id, begin, end);
/*    if (iend!=-1 && strlen(header)<=HDRLEN) fprintf(stdout, "%s", header);
      else if (iend!=-1) for (i=0; i<HDRLEN; i++) fputc(header[i], stdout); */
      fprintf(stdout, " entropy=%4.2f (%d/%4.2f/%4.2f)\n",
           seq->entropy, window, locut, hicut);
     }
   else
     {
      fprintf(stdout, ">%s(%d-%d)", id, begin, end);
      if (iend!=-1)   /* fprintf(stdout, "%s\n", header); */
        {
         strncpy(outbuf, header, HDRLEN);
         outbuf[HDRLEN] = '\0';
         fprintf(stdout, "%s\n", outbuf);
        }
      else fputc('\n', stdout);
     }
   
   if (hilo==lo)
     {
      lower(proseq);
     }
   else if (hilo==hi && seq->length>=hilenmin)
     {
      upper(proseq);
     }
   else
     {
      lower(proseq);
     }

   for (i=0, ctr=0; proseq[i]!='\0'; i++, ctr++)
     {
      if (ctr==charline)
        {
         fputc('\n', stdout);
         ctr = 0;
        }
      fputc(proseq[i], stdout);
     }

   fputc('\n', stdout);
   fputc('\n', stdout);
  }

/*-------------------------------------------------------------(appendseg)---*/

appendseg(segs, seg)
  struct Segment *segs, *seg;

  {struct Segment *temp;

   temp = segs;
   while (1)
     {
      if (temp->next==NULL)
        {
         temp->next = seg;
         break;
        }
      else
        {
         temp = temp->next;
        }
     }

   return;
  }

/*--------------------------------------------------------------(freesegs)---*/

freesegs(segs)
  struct Segment *segs;

  {struct Segment *temp;

   while (segs!=NULL)
     {
      temp = segs->next;
      free(segs);
      segs = temp;
     }
  }

/*-----------------------------------------------------------------(usage)---*/

usage()

  {
   fprintf(stderr, "\
Usage: seg <file> <window> <locut> <hicut> <options>\n\
         <file>   - filename containing fasta-formatted sequence(s) \n\
         <window> - OPTIONAL window size (default 12) \n\
         <locut>  - OPTIONAL low (trigger) entropy (default 2.2) \n\
         <hicut>  - OPTIONAL high (extension) entropy (default 2.5) \n\
         <options> \n\
            -m <size> minimum length for a high-entropy segment (default 0) \n\
                shorter segments are merged with adjacent low-entropy segments \n\
            -o  show overlapping low-entropy segments (default merge) \n\
            -t <maxtrim> maximum segment trimming (default 100) \n\
            -h  show only high-entropy segments (fasta format) \n\
            -l  show only low-entropy segments (fasta format) \n\
            -a  show all segments (fasta format) \n\
            -n  do not add entropy information to the header line \n\
            -x  each input sequence is represented by a single output sequence\n\
                with low-entropy regions replaced by strings of 'X'-es.\n\
            -c <chars> number of sequence characters/line (default 60)\n\
            -p  prettyprint each segmented sequence (tree format) \n\
            -q  prettyprint each segmented sequence (block format) \n");
   exit(1);
  }


/*---------------------------------------------------------------------------*/
