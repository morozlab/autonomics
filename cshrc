set history = 5000
set savehist = (5000 merge)
# set histfile = ~/.tcsh_history
set ignoreeof
set prompt="`hostname`{$LOGNAME}\!: "
set notify
setenv PERL5LIB /srv/data2/pwilliams/pantherScore1.03/lib
source ~/.aliases
set noclobber
setenv EDITOR emacs

# acis

setenv PATH .:/srv/data2/pwilliams/hmmer-2.3.2/bin/bin:/srv/data2/pwilliams/ncbi-blast-2.2.24+/bin:/usr/local/bin:/home/pwilliams/scripts:/home/pwilliams/bin:/home/pwilliams/bin/bin:/usr/sbin:/srv/data2/pwilliams/MUMmer3.23:$PATH

setenv CUTADAPT_PATH /usr/local/bin
setenv MIRA_PATH /usr/bin
setenv PFAM_LOCAL_PATH /srv/data/pfam/PfamScan
setenv TRINITY_PATH /srv/data2/pwilliams/trinityrnaseq_r20131110
setenv KHMER_PYTHON_PATH /srv/data2/software/Khmer/python
setenv KHMER_PATH /srv/data2/software/Khmer/scripts
setenv PFAM_HPC_DATA_PATH /scratch/lfs/moroz/pfam/
setenv HPC_USER_NAME plw1080
setenv REMOTE_DIR /scratch/lfs/moroz/
if ($?PYTHONPATH) then
  setenv PYTHONPATH $KHMER_PYTHON_PATH/:$PYTHONPATH
else
  setenv PYTHONPATH $KHMER_PYTHON_PATH
endif
setenv PATH_TO_AUTONOMICS /home/pwilliams
setenv AUTONOMICS_PATH $PATH_TO_AUTONOMICS/autonomics
setenv PYTHONPATH $PATH_TO_AUTONOMICS/:$AUTONOMICS_PATH/:$PYTHONPATH
setenv PROJECT_PATH /srv/data2/pipeline
setenv PATH $AUTONOMICS_PATH/bin:$AUTONOMICS_PATH/prediction:$PATH
if ($?CLASSPATH) then
  setenv CLASSPATH $AUTONOMICS_PATH/java/bin:$AUTONOMICS_PATH/java/mysql-connector-java-5.0.8-bin.jar:$CLASSPATH
else
  setenv CLASSPATH $AUTONOMICS_PATH/java/bin:$AUTONOMICS_PATH/java/mysql-connector-java-5.0.8-bin.jar
endif
setenv USING_NB 0
setenv RUNNING_AUTONOMICS_PIPELINE 1
