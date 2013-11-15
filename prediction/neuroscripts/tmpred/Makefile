PC  = fpc
PFLAGS = -Dosf1 -Ddec -C all -O5 -math_library fast

all: tmpred btmpred

pp_osf1: pp_osf1.pas
	$(PC) -C bounds -o pp_osf1 pp_osf1.pas

dep_dec.o: dep_dec.pas 
	$(PC) $(PFLAGS) -g0 -c dep_dec.pas

ureadseq.o: pp_osf1 ureadseq.pas dep_dec.o
	pp_osf1 ureadseq.pas
	$(PC) $(PFLAGS)  -c ureadseq.p

seq.o: pp_osf1 seq.pas ureadseq.o dep_dec.o
	pp_osf1 seq.pas
	$(PC) $(PFLAGS)  -c seq.p

tabutil.o: pp_osf1 tabutil.pas dep_dec.o
	pp_osf1 tabutil.pas
	$(PC) $(PFLAGS)  -c tabutil.p

test: pp_osf1 test.pas seq.o ureadseq.o dep_dec.o
	pp_osf1 test.pas
	$(PC) $(PFLAGS) -o test test.p seq.o ureadseq.o dep_dec.o

tmpred: pp_osf1 tmpred.pas tabutil.o seq.o ureadseq.o dep_dec.o
	pp_osf1 tmpred.pas
	$(PC) $(PFLAGS) -o tmpred tmpred.p tabutil.o seq.o ureadseq.o dep_dec.o

btmpred: pp_osf1 btmpred.pas tabutil.o seq.o ureadseq.o dep_dec.o
	pp_osf1 btmpred.pas
	$(PC) $(PFLAGS) -o btmpred btmpred.p tabutil.o seq.o ureadseq.o dep_dec.o

