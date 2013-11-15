
/*****************************************************************************/
/***   (genwin.c)                                                          ***/
/*****************************************************************************/

/*--------------------------------------------------------------(includes)---*/

#include "genwin.h"

/*---------------------------------------------------------------(defines)---*/

#ifndef EOF
#define EOF '\377'
#endif

#define STRSIZE 100

/*----------------------------------------------------------------(protos)---*/

struct Sequence *readentry();

/*---------------------------------------------------------------(globals)---*/

char *blastdbs[] =
  {"bba", "bbn", "embl", "gbupdate", "genbank", "genpept", "gpupdate",
   "nr", "nrdb", "nrdb.shuf", "pir", "pseq", "swissprot", "tfdaa"};

int nblastdbs = 14;

char *blastdir = "/net/cruncher/usr/ncbi/db/fasta/";
char *indexdir = "/net/cruncher/usr/ncbi/db/index/";

int nabets;
struct Alphabet **abets;
int ntvecs;
struct TransVector **tvecs;
int nsvecs;
struct ScoreVector **svecs;
int nsmats;
struct ScoreMatrix **smats;

int aaindex[128];
char aachar[20];

struct strlist
  {
   char string[STRSIZE];
   struct strlist *next;
  } *str, *curstr;

/*---------------------------------------------------------------(tmalloc)---*/

#define TESTMAX 1000
void *tmalloc();
int record_ptrs[TESTMAX] = {0,0,0,0};
int rptr = 0;

/*------------------------------------------------------------(genwininit)---*/

genwininit()

  {int i;
   char c;

   for (i=0; i<128; i++)
     {
      c = (char) i;
      if (c=='a' || c=='A') {aaindex[i] = 0; aachar[0] = c;}
      else if (c=='c' || c=='C') {aaindex[i] = 1; aachar[1] = c;}
      else if (c=='d' || c=='D') {aaindex[i] = 2; aachar[2] = c;}
      else if (c=='e' || c=='E') {aaindex[i] = 3; aachar[3] = c;}
      else if (c=='f' || c=='F') {aaindex[i] = 4; aachar[4] = c;}
      else if (c=='g' || c=='G') {aaindex[i] = 5; aachar[5] = c;}
      else if (c=='h' || c=='H') {aaindex[i] = 6; aachar[6] = c;}
      else if (c=='i' || c=='I') {aaindex[i] = 7; aachar[7] = c;}
      else if (c=='k' || c=='K') {aaindex[i] = 8; aachar[8] = c;}
      else if (c=='l' || c=='L') {aaindex[i] = 9; aachar[9] = c;}
      else if (c=='m' || c=='M') {aaindex[i] = 10; aachar[10] = c;}
      else if (c=='n' || c=='N') {aaindex[i] = 11; aachar[11] = c;}
      else if (c=='p' || c=='P') {aaindex[i] = 12; aachar[12] = c;}
      else if (c=='q' || c=='Q') {aaindex[i] = 13; aachar[13] = c;}
      else if (c=='r' || c=='R') {aaindex[i] = 14; aachar[14] = c;}
      else if (c=='s' || c=='S') {aaindex[i] = 15; aachar[15] = c;}
      else if (c=='t' || c=='T') {aaindex[i] = 16; aachar[16] = c;}
      else if (c=='v' || c=='V') {aaindex[i] = 17; aachar[17] = c;}
      else if (c=='w' || c=='W') {aaindex[i] = 18; aachar[18] = c;}
      else if (c=='y' || c=='Y') {aaindex[i] = 19; aachar[19] = c;}
      else {aaindex[i] = 20;}
     }

   return;
  }
        
/*-------------------------------------------------------------(opendbase)---*/

extern struct Database *opendbase(name)
  char *name;

  {struct Database *dbase;

   dbase = (struct Database *) malloc(sizeof(struct Database));

   if (blastdb(name))
     {
      dbase->filename = (char *) malloc(strlen(blastdir)+strlen(name)+1);
      dbase->indexname = (char *) malloc(strlen(indexdir)+strlen(name)+1);
      strcpy(dbase->filename, blastdir);
      strcat(dbase->filename, name);
      strcpy(dbase->indexname, indexdir);
      strcat(dbase->indexname, name);
     }
   else
     {
      dbase->filename = (char *) malloc(strlen(name)+1);
      dbase->indexname = (char *) malloc(strlen(name)+1);
      strcpy(dbase->filename, name);
      strcpy(dbase->indexname, name);
     }

   if (strcmp(dbase->filename, "-")==0)
     {
      dbase->fp = stdin;
     }
   else if ((dbase->fp=fopen(dbase->filename, "r"))==NULL)
     {
      free(dbase->filename);
      free(dbase->indexname);
      free(dbase);
      return((struct Database *) NULL);
     }

   dbase->db = dbm_open(dbase->indexname, O_RDONLY);

   dbase->filepos = 0L;

   return(dbase);
  }

/*---------------------------------------------------------------(blastdb)---*/

int blastdb(name)
  char *name;

  {int i;

   for (i=0; i<nblastdbs; i++)
     {
      if (strcmp(name, blastdbs[i])==0) {return(TRUE);}
     }

   return(FALSE);
  }

/*------------------------------------------------------------(closedbase)---*/

extern closedbase(dbase)
  struct Database *dbase;

  {
   fclose(dbase->fp);
   if (dbase->db!=NULL) dbm_close(dbase->db);

   free(dbase->filename);
   free(dbase->indexname);
   free(dbase);

   return;
  }

/*---------------------------------------------------------------(openseq)---*/

extern struct Sequence *openseq(dbase, id)
  struct Database *dbase;
  char *id;

  {datum key, content;
   struct Sequence *seq;
   char *buffer;
   int derror;

   if (dbase->db==NULL)
     {fprintf(stderr, "No index files for database %s.\n", dbase->filename);
      exit(1);}

   key.dptr = (char *) malloc(strlen(id)+1);
   strcpy(key.dptr, id);
   upper(key.dptr);
   key.dsize = strlen(id);
   content = dbm_fetch(dbase->db, key);
   if (content.dptr==NULL)
     {fprintf(stderr, "Entry %s not found.\n", id);
      return((struct Sequence *) NULL);}

   buffer = (char *) malloc(content.dsize+1);
   strncpy(buffer, content.dptr, content.dsize);
   buffer[content.dsize] = '\0';

/* content.dptr[content.dsize] = '\0'; bug in ndbm? *//* Nope. dope in ncbi */

   dbase->filepos = strtol(content.dptr, NULL, 10);
   if (fseek(dbase->fp, dbase->filepos, 0)==-1)
     {fprintf(stderr, "Error fseek-ing position in file.\n");
      exit(1);}

   free(buffer);
   free(key.dptr);

   seq = readentry(dbase);

   return(seq);
  }

/*--------------------------------------------------------------(firstseq)---*/

extern struct Sequence *firstseq(dbase)
  struct Database *dbase;

  {
   if (dbase->filepos!=0L)
     {
      dbase->filepos = 0L;
      if (fseek(dbase->fp, dbase->filepos, 0)!=0)
        {fprintf(stderr, "Error positioning file %s for firstseq.\n",
                           dbase->filename);
         exit(1);}
     }

   return(readentry(dbase));
  }

/*---------------------------------------------------------------(nextseq)---*/

extern struct Sequence *nextseq(dbase)
  struct Database *dbase;

  {
   return(readentry(dbase));
  }

/*--------------------------------------------------------------(closeseq)---*/

extern closeseq(seq)
  struct Sequence *seq;

  {
   if (seq==NULL) return;

   if (seq->id!=NULL)          free(seq->id);
   if (seq->name!=NULL)        free(seq->name);
   if (seq->organism!=NULL)    free(seq->organism);
   if (seq->header!=NULL)      free(seq->header);
   if (seq->state!=NULL)       free(seq->state);
   if (seq->composition!=NULL) free(seq->composition);

   free(seq->seq);

   free(seq->config);

   free(seq);
   return;
  }

/*---------------------------------------------------------------(openwin)---*/

extern struct Sequence *openwin(parent, start, length)
  struct Sequence *parent;
  int start, length;

  {struct Sequence *win;
   int i;

   if (start<0 || length<0 || start+length>parent->length)
     {
      return((struct Sequence *) NULL);
     }

   win = (struct Sequence *) malloc(sizeof(struct Sequence));

/*---                                          ---[set links, up and down]---*/

   win->parent = parent;
   if (parent->root==NULL)
     {win->root = parent;}
   else
     {win->root = parent->root;}
   win->children = (struct Sequence **) NULL;

/* parent->children = ***foo***                   ---[not yet implemented]---*/

   win->id = (char *) NULL;
   win->name = (char *) NULL;
   win->organism = (char *) NULL;
   win->header = (char *) NULL;

/*---                          ---[install the local copy of the sequence]---*/

   win->start = start;
   win->length = length;
   win->seq = (char *) malloc(sizeof(char)*length + 1);
   memcpy(win->seq, (parent->seq)+start, length);
   win->seq[length] = '\0';

/*---                          ---[setup window implementation parameters]---*/

   win->config = (struct Configuration *) malloc(sizeof(struct Configuration));

   win->config->iseq = win->seq;
   win->config->ilength = win->length;
   win->config->printper = 60;

/*---                                                 ---[set local flags]---*/

   win->rubberwin = FALSE;
   win->floatwin = FALSE;
   win->punctuation = FALSE;
   win->bogus = FALSE;
   for (i=0; i<length; i++)
     {
      if (bogus(win->seq[i]))
        {
         win->bogus = TRUE;
         break;
        }
     }

/*---                                   ---[initially unconfiguerd window]---*/

   win->entropy = -2.;
   win->state = (int *) NULL;
   win->composition = (int *) NULL;
   win->classvec = (char *) NULL;
   win->scorevec = (double *) NULL;

   return(win);
  }

/*---------------------------------------------------------------(nextwin)---*/

extern struct Sequence *nextwin(win, shift)
  struct Sequence *win;
  int shift;

  {
   if ((win->start+shift)<0 ||
       (win->start+win->length+shift)>win->parent->length)
     {
      return((struct Sequence *) NULL);
     }
   else
     {
      return(openwin(win->parent, win->start+shift, win->length));
     }
  }

/*--------------------------------------------------------------(shiftwin)---*/

extern int shiftwin(win, shift)
  struct Sequence *win;
  int shift;

  {int i, j;

/*---                                  ---[check boundaries of new window]---*/

   if (!win->rubberwin)                                   /*---[rigidwins]---*/
     {
      if ((win->start + shift)<0 ||
          (win->start + win->length + shift) > win->parent->length)
        {
         return(FALSE);
        }
     }
   else                                                  /*---[rubberwins]---*/
     {
      if ((win->start + win->length + shift)<=0 ||
          (win->start + shift) > win->parent->length)
        {
         return(FALSE);
        }
     }

/*---                    ---[half-update configured fixed function fields]---*/

   if (win->composition!=NULL)
     {
      if (shift>0)
        {
         for (i=0; i<shift && i<win->length; i++)
           {
            if (!isalpha(win->seq[i]) || bogus(win->seq[i]))
              {
               continue;
              }
            if (win->state!=NULL)
              {
               decrementsv(win->state, win->composition[aaindex[win->seq[i]]]);
              }
            win->composition[aaindex[win->seq[i]]]--;
           }
        }
      else if (shift<0)
        {
         for (i=0; i<shift && i<win->length; i++)
           {
            j = win->length - i - 1;
            if (!isalpha(win->seq[j]) || bogus(win->seq[j]))
              {
               continue;
              }
            if (win->state!=NULL)
              {
               decrementsv(win->state, win->composition[aaindex[win->seq[j]]]);
              }
            win->composition[aaindex[win->seq[i]]]--;
           }
        }
     }
   

/*---                                       ---[shift in the new sequence]---*/

   if (!win->rubberwin)                                   /*---[rigidwins]---*/
     {
      win->start += shift;
      memcpy(win->seq, (win->parent->seq) + win->start, win->length);
     }
   else                                                  /*---[rubberwins]---*/
     {
      fprintf(stderr, "Rubber windows not yet implemented.\n");
      exit(1);
     }

/*---                                                    ---[update flags]---*/

   win->bogus = FALSE;
   for (i=0; i < win->length; i++)
     {
      if (bogus(win->seq[i]))
        {
         win->bogus = TRUE;
        }
     }

/*---                          ---[half-update configured fixed functions]---*/

   if (win->composition!=NULL)
     {
      if (shift>0)
        {
         for (i=0; i<shift && i<win->length; i++)
           {
            j = win->length - i - 1;
            if (!isalpha(win->seq[j]) || bogus(win->seq[j]))
              {
               continue;
              }
            if (win->state!=NULL)
              {
               incrementsv(win->state, win->composition[aaindex[win->seq[j]]]);
              }
            win->composition[aaindex[win->seq[j]]]++;
           }
        }
      else if (shift>0)
        {
         for (i=0; i<shift && i<win->length; i++)
           {
            if (!isalpha(win->seq[i]) || bogus(win->seq[i]))
              {
               continue;
              }
            if (win->state!=NULL)
              {
               incrementsv(win->state, win->composition[aaindex[win->seq[i]]]);
              }
            win->composition[aaindex[win->seq[i]]]++;
           }
        }
     }

                /*---[#ifdef ENTROPY, then table-lookup or theTree lookup]---*/

   if (win->entropy>-2.)
     {
      if (win->state==NULL) {stateon(win);}
      win->entropy = entropy(win->state);
     }

/*---                             ---[update configured generic functions]---*/
/*                                          [not yet implemented]------------*/

/*---                                      ---[update floating subwindows]---*/
/*                                             ---[not yet implemented]------*/

   return(TRUE);
  }

/*--------------------------------------------------------------(closewin)---*/

extern closewin(win)
  struct Sequence *win;

  {
   if (win==NULL) return;

   free(win->config->iseq);
   free(win->config);

   if (win->state!=NULL)       free(win->state);
   if (win->composition!=NULL) free(win->composition);
   if (win->classvec!=NULL)    free(win->classvec);
   if (win->scorevec!=NULL)    free(win->scorevec);

   free(win);
   return;
  }

/*----------------------------------------------------------------(compon)---*/

extern compon(win)
  struct Sequence *win;

  {int i, aa;

   win->composition = (int *) malloc(20*sizeof(int));

   for (aa=0; aa<20; aa++)
     {
      win->composition[aa] = 0;
     }

   for (i=0; i<win->length; i++)
     {
      if (!isalpha(win->seq[i]) || bogus(win->seq[i]))
        {
         continue;
        }
      win->composition[aaindex[win->seq[i]]]++;
     }

   return;
  }

/*---------------------------------------------------------------(stateon)---*/

extern stateon(win)
  struct Sequence *win;

  {int aa, i, j;

   if (win->composition==NULL) {compon(win);}

   win->state = (int *) malloc(21*sizeof(int));

   for (aa=0; aa<21; aa++)
     {
      win->state[aa] = 0;
     }

   for (aa=0; aa<20; aa++)
     {
      if (win->composition[aa]==0) {continue;}

      for (i=0; i<20; i++)
        {
         if (win->composition[aa]<win->state[i]) {continue;}

         for (j=19; j>i; j--)
           {
            win->state[j] = win->state[j-1];
           }

         win->state[i] = win->composition[aa];
         break;
        }
     }

   return;
  }

/*-----------------------------------------------------------------(enton)---*/

extern enton(win)
  struct Sequence *win;

  {
   if (win->state==NULL) {stateon(win);}

   win->entropy = entropy(win->state);

   return;
  }

/*---------------------------------------------------------------(entropy)---*/

extern double entropy(sv)
  int *sv;

  {double ent;
   int i, total;

   total = 0;
   for (i=0; sv[i]!=0; i++)
     {
      total += sv[i];
     }
   if (total==0) return(0.);

   ent = 0.0;
   for (i=0; sv[i]!=0; i++)
     {
      ent += ((double)sv[i])*log(((double)sv[i])/(double)total)/log(2.);
     }

   ent = fabs(ent/(double)total);

   return(ent);
  }

/*-----------------------------------------------------------(decrementsv)---*/

decrementsv(sv, class)
  int *sv;
  int class;

  {int i;

   for (i=0; sv[i]!=0; i++)
     {
      if (sv[i]==class && sv[i+1]<class)
        {
         sv[i]--;
         break;
        }
     }
  }

/*-----------------------------------------------------------(incrementsv)---*/

incrementsv(sv, class)
  int *sv;
  int class;

  {int i;

   for (i=0; ; i++)
     {
      if (sv[i]==class)
        {
         sv[i]++;
         break;
        }
     }
  }

/*---------------------------------------------------------------(openmat)---*/

extern struct Matrix *openmat(parent, start, period, length)
  struct Sequence *parent;
  int start, period, length;

  {struct Matrix *mat;
   int i, j, istart, ilength;
   
   if (start<0 || start>parent->length || period<=0 || length<=0 ||
       (start+(period*(length-1)))>=parent->length)
     {
      return((struct Matrix *) NULL);
     }

   mat = (struct Matrix *) malloc(sizeof(struct Matrix));

   mat->parent = parent;

   mat->start = start;
   mat->period = period;
   mat->length = length;

   mat->seq = (char **) malloc(length*sizeof(char *));
   for (i=0, istart=start; i<length; i++, istart+=period)
     {
      mat->seq[i] = (char *) malloc(period*sizeof(char)+1);
      
      ilength = min(period, (parent->length)-istart);
      memcpy(mat->seq[i], (parent->seq)+istart, ilength);
      mat->seq[i][ilength] = '\0';
      if (ilength<period)
        {
         for (j=ilength; j<=period; j++)
           {
            mat->seq[i][j] = '\0';
           }
        }
     }

   return(mat);
  }

/*--------------------------------------------------------------(closemat)---*/

extern closemat(mat)
  struct Matrix *mat;

  {int i;

   for (i=0; i<mat->length; i++)
     {
      free(mat->seq[i]);
     }

   free(mat->seq);

   free(mat);

   return;
  }

/*-------------------------------------------------------------(readentry)---*/

struct Sequence *readentry(dbase)
  struct Database *dbase;

  {struct Sequence *seq;
   char c;

   seq = (struct Sequence *) malloc(sizeof(struct Sequence));

   seq->db = dbase;

/*---                                    ---[backpointers null at the top]---*/

   seq->parent = (struct Sequence *) NULL;
   seq->root = (struct Sequence *) NULL;
   seq->children = (struct Sequence **) NULL;

/*---                                                       ---[set flags]---*/

   seq->rubberwin = FALSE;
   seq->floatwin = FALSE;

/*---                                                  ---[read from file]---*/

   if (!readhdr(seq))
     {
      return((struct Sequence *) NULL);
     }
   while (1)  /*---[skip multiple headers]---*/
     {
      c = fgetc(dbase->fp);
      if (isalpha(c))
        {
         ungetc(c, dbase->fp);
         break;
        }
      while ((c=fgetc(dbase->fp))!='\n');
     }
   readseq(seq);

/*---                                   ---[set implementation parameters]---*/

   seq->config = (struct Configuration *) malloc(sizeof(struct Configuration));

   seq->config->iseq = (char *) seq->seq;
   seq->config->ilength = seq->length;
   seq->config->printper = 60;

/*---                                          ---[initially unconfigured]---*/

   seq->entropy = -2.;
   seq->state = (int *) NULL;
   seq->composition = (int *) NULL;
   seq->classvec = (char *) NULL;
   seq->scorevec = (double *) NULL;

   return(seq);
  }

/*---------------------------------------------------------------(readhdr)---*/

readhdr(seq)
  struct Sequence *seq;

  {FILE *fp;
   char *bptr, *curpos, c;
   int i, itotal;
   int idend, namend, orgend;

   fp = seq->db->fp;

   if ((c=fgetc(fp))==EOF)
     {
      free(seq);
      return(FALSE);
     }
   
   while (isspace(c))
     {
      c = fgetc(fp);
     }

   if (c!='>')
     {fprintf(stderr, "Error reading fasta format - '>' not found.\n");
      exit(1);}
   ungetc(c, fp);
/*                                               ---[read the header line]---*/
   str = (struct strlist *) malloc (sizeof(struct strlist));
   str->next = NULL;
   curstr = str;

   for (i=0,itotal=0,c=fgetc(fp); ; c=fgetc(fp))
     {
      if (c=='\n') break;

      if (i==STRSIZE-1)
        {curstr->string[i] = '\0';
         curstr->next = (struct strlist *) malloc (sizeof(struct strlist));
         curstr = curstr->next;
         curstr->next = NULL;
         i = 0;}

      curstr->string[i] = c;
      itotal++;
      i++;
     }

   curstr->string[i] = '\0';
   seq->header = (char *) malloc (itotal+2);
   seq->header[0] = '\0';

   for (curstr=str, curpos=seq->header; curstr!=NULL;)
     {
      if (curstr->next==NULL)
        {memccpy(curpos, curstr->string, '\0', STRSIZE);}
      else
        {memccpy(curpos, curstr->string, '\0', STRSIZE-1);}

      str = curstr;
      curstr = curstr->next;
      free (str);

      if (curstr!=NULL) {curpos = curpos+STRSIZE-1;}
     }

   bptr = (seq->header)+1;
   seq->name = (char *) NULL;
   seq->organism = (char *) NULL;
/*                                                   ---[parse out the id]---*/
   idend = findchar(bptr, ' ');
   if (idend==-1) {idend = findchar(bptr, '\n');}
   if (idend==-1) {idend = findchar(bptr, '\0');}
   if (idend==-1)
     {fprintf(stderr, "Error parsing header line - id.\n");
      fputs(seq->header, fp);
      exit(1);}   

   seq->id = (char *) malloc((idend+1)*sizeof(char));
   memcpy(seq->id, bptr, idend);
   seq->id[idend] = '\0';

   if (bptr[idend]=='\n' || bptr[idend]=='\0') {return(TRUE);}

/*                                         ---[parse out the protein name]---*/
   bptr = bptr + idend + 1;
   while (bptr[0]==' ') {bptr++;}

   namend = findchar(bptr, '-');
   if (namend==-1) {namend = findchar(bptr, '\n');}
   if (namend==-1) {namend = findchar(bptr, '\0');}
   if (namend==-1)
     {fprintf(stderr, "Error parsing header line - name.\n");
      fputs(seq->header, fp);
      return(TRUE);}

   seq->name = (char *) malloc((namend+1)*sizeof(char));
   memcpy(seq->name, bptr, namend);
   seq->name[namend] = '\0';

   if (bptr[namend]=='\n' || bptr[namend]=='\0') {return(TRUE);}

/*                                                 ---[parse out organism]---*/
   bptr = bptr + namend + 1;
   while (bptr[0]==' ') {bptr++;}

   orgend = findchar(bptr, '|');
   if (orgend==-1) {orgend = findchar(bptr, '#');}
   if (orgend==-1) {orgend = findchar(bptr, '\n');}
   if (orgend==-1) {orgend = findchar(bptr, '\0');}
   if (orgend==-1)
     {fprintf(stderr, "Error parsing header line - organism.\n");
      fputs(seq->header, fp);
      return(TRUE);}

   seq->organism = (char *) malloc((orgend+1)*sizeof(char));
   memcpy(seq->organism, bptr, orgend);
   seq->organism[orgend] = '\0';

/*                                    ---[skip over multiple header lines]---*/
   while (TRUE)
     {
      c = fgetc(fp);
      if (c=='>')
        {
         skipline(fp);
        }
      else
        {
         ungetc(c,fp);
         break;
        }
     }

   return(TRUE);
  }

/*--------------------------------------------------------------(skipline)---*/

skipline(fp)
  FILE *fp;

  {char c;

   while ((c=fgetc(fp))!='\n' && c!=EOF)
     {}

   return;
  }

/*--------------------------------------------------------------(findchar)---*/

extern int findchar(str, chr)
  char *str;
  char chr;

  {int i;

   for (i=0; ; i++)
     {
      if (str[i]==chr)
        {
         return(i);
        }
      if (str[i]=='\0')
        {
         return(-1);
        }
     }
   }

/*---------------------------------------------------------------(readseq)---*/

readseq(seq)
  struct Sequence *seq;

  {FILE *fp;
   int i, itotal;
   char c;
   char *curpos;

   fp = seq->db->fp;

   seq->bogus = FALSE;
   seq->punctuation = FALSE;   

   str = (struct strlist *) malloc (sizeof(struct strlist));
   str->next = NULL;
   curstr = str;

   for (i=0,itotal=0,c=fgetc(fp); ; c=fgetc(fp))
     {
      if (!isalpha(c))
        {
         if (c=='>')
           {ungetc(c, fp);
            break;}
         if (c==EOF)
           {break;}
         if (c=='-')
           {seq->punctuation = TRUE;}
         else if (c=='*')
           {seq->punctuation = TRUE;}
         else
           {seq->punctuation = TRUE;
            continue;}
        }

      if (i==STRSIZE-1)
        {curstr->string[i] = '\0';
         curstr->next = (struct strlist *) malloc (sizeof(struct strlist));
         curstr = curstr->next;
         curstr->next = NULL;
         i = 0;}

      if (bogus(c)) {seq->bogus = TRUE;}

      curstr->string[i] = c;
      itotal++;
      i++;
     }

   curstr->string[i] = '\0';
   seq->seq = (char *) malloc (itotal+2);
   seq->seq[0] = '\0';

   for (curstr=str, curpos=seq->seq; curstr!=NULL;)
     {
      if (curstr->next==NULL)
        {memccpy(curpos, curstr->string, '\0', STRSIZE);}
      else
        {memccpy(curpos, curstr->string, '\0', STRSIZE-1);}

      str = curstr;
      curstr = curstr->next;
      free (str);

      if (curstr!=NULL) {curpos = curpos+STRSIZE-1;}
     }

   seq->length = strlen(seq->seq);

   return;
  }
/*-----------------------------------------------------------------(upper)---*/

extern upper(string)
  char *string;

  {int i;

   for (i=0; string[i]!='\0'; i++)
     {
      if (islower(string[i]))
        {
         string[i] = toupper(string[i]);
        }
     }
  }

/*-----------------------------------------------------------------(lower)---*/

extern lower(string)
  char *string;
  
  {int i;

   for (i=0; string[i]!='\0'; i++)
     {
      if (isupper(string[i]))
        {
         string[i] = tolower(string[i]);
        }
     }
  }

/*-------------------------------------------------------------------(min)---*/

int min(a, b)
  int a, b;

  {
   if (a<b) {return(a);}
   else {return(b);}
  }

/*-------------------------------------------------------------------(max)---*/

int max(a, b)
  int a, b;

  {
   if (a<b) {return(b);}
   else {return(a);}
  }

/*-----------------------------------------------------------------(bogus)---*/
/*                                         ---[try a macro implementation]---*/

int bogus(c)
  char c;

  {
   if (aaindex[c]>=20) {return(TRUE);}
   else {return(FALSE);}
  }

/*                                                                           */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------(tmalloc)---*/

void *tmalloc(size)
  size_t size;

  {void *ptr;

   ptr = (void *) malloc(size);

   if (rptr>TESTMAX)
     {
      exit(2);
     }

   record_ptrs[rptr] = (int) ptr;
   rptr++;

   return(ptr);
  }

/*-----------------------------------------------------------------(tfree)---*/

tfree(ptr)
  void *ptr;

  {int i;

   for (i=0; i<rptr; i++)
     {
      if (record_ptrs[i]==(int)ptr)
        {
         record_ptrs[i] = 0;
         break;
        }
      }

   free(ptr);
  }

/*---------------------------------------------------------------------------*/
