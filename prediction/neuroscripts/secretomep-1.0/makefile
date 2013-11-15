
# makefile for SecretomeP 1.0

PFUN	= /usr/cbs/bio/src/secretomep-1.0

PARSERS	= parsers/allprot/Process.server \
	  parsers/prop/Process.server \
	  parsers/psort/Process.server \
	  parsers/seg/Process.server \
	  parsers/tmhmm/Process.server

BINDIR	= /usr/cbs/bio/bin
MANDIR	= /usr/cbs/bio/man/man1


# development and debugging ===================================================

all:	$(PARSERS)

parsers/allprot/Process.server: parsers/allprot/Process.server.in bin/lib.gawk
	cat bin/lib.gawk > parsers/allprot/Process.server
	cat parsers/allprot/Process.server.in >> parsers/allprot/Process.server

parsers/seg/Process.server: parsers/seg/Process.server.in bin/lib.gawk
	cat bin/lib.gawk > parsers/seg/Process.server
	cat parsers/seg/Process.server.in >> parsers/seg/Process.server

parsers/tmhmm/Process.server: parsers/tmhmm/Process.server.in bin/lib.gawk
	cat bin/lib.gawk > parsers/tmhmm/Process.server
	cat parsers/tmhmm/Process.server.in >> parsers/tmhmm/Process.server

parsers/psort/Process.server: parsers/psort/Process.server.in bin/lib.gawk
	cat bin/lib.gawk > parsers/psort/Process.server
	cat parsers/psort/Process.server.in >> parsers/psort/Process.server

parsers/prop/Process.server: parsers/prop/Process.server.in bin/lib.gawk
	cat bin/lib.gawk > parsers/prop/Process.server
	cat parsers/prop/Process.server.in >> parsers/prop/Process.server


# installation ================================================================
install:	secretomep
		cp secretomep $(BINDIR)
		chmod 755 $(BINDIR)/secretomep

install.man:	secretomep.1
		cp secretomep.1 $(MANDIR)
		chmod 444 $(MANDIR)/secretomep.1

# other =======================================================================
clean:
		-rm -rf $(PFUN)/tmp/*


# end of makefile =============================================================
