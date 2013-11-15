
all : genwin seg

genwin : genwin.c genwin.h
	cc -c genwin.o genwin.c -lm -ldbm

seg : seg.c lnfac.h genwin.h genwin.o
	cc -o seg seg.c genwin.o -lm -ldbm
	cp seg /home/federhen/bin.sun4

