prefix=/usr
bindir=$(prefix)/bin
datadir=$(prefix)/share

FFLAGS=-O3 -ffast-math -march=native -g -fbacktrace

build: bin/aqfindstar

bin/aqfindstar: bin/%: imutil.o %.f90
	mkdir -p bin
	gfortran -o $@ -lcfitsio $^

imutil.o: imutil.f90
	gfortran -c -o $@ $^

installdir:
	install -dvZ $(DESTDIR)$(datadir)/aquila
	install -dvZ $(DESTDIR)$(bindir)

install: build installdir
	install -pvZ -m 644 python/aqutils.py $(DESTDIR)$(datadir)/aquila
	install -pvZ -m 755 python/aqraw2fit.py $(DESTDIR)$(datadir)/aquila
	ln -srf $(DESTDIR)$(datadir)/aquila/aqraw2fit.py $(DESTDIR)$(bindir)/aqraw2fit
	install -pvZ -m 755 python/aqcfaflat.py $(DESTDIR)$(datadir)/aquila
	ln -srf $(DESTDIR)$(datadir)/aquila/aqcfaflat.py $(DESTDIR)$(bindir)/aqcfaflat
	install -pvZ -m 755 python/aqenum.sh $(DESTDIR)$(datadir)/aquila
	ln -srf $(DESTDIR)$(datadir)/aquila/aqenum.sh $(DESTDIR)$(bindir)/aqenum
	install -pvZ -m 755 python/aqfngen.py $(DESTDIR)$(datadir)/aquila
	ln -srf $(DESTDIR)$(datadir)/aquila/aqfngen.py $(DESTDIR)$(bindir)/aqfngen
	install -pvZ -m 755 python/aqcfa2rggb.py $(DESTDIR)$(datadir)/aquila
	ln -srf $(DESTDIR)$(datadir)/aquila/aqcfa2rggb.py $(DESTDIR)$(bindir)/aqcfa2rggb
	install -pvZ -m 755 python/aqrggb2grey.py $(DESTDIR)$(datadir)/aquila
	ln -srf $(DESTDIR)$(datadir)/aquila/aqrggb2grey.py $(DESTDIR)$(bindir)/aqrggb2grey
	install -pvZ -m 755 bin/aqfindstar $(DESTDIR)$(bindir)

clean:
	rm -f *.o *.mod
	rm -rf bin

distclean: clean