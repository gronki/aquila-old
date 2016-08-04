
BINARIES=bin/test-align bin/test-wavele

all: $(BINARIES)

$(BINARIES): bin/%: %.f90 ../imutil.o
	mkdir -p bin
	mkdir -p temp
	gfortran -o $@ -I.. -L.. -lcfitsio $^


../imutil.o: ../imutil.f90
	cd ..; make


clean:
	rm -f *.o
	rm -rf bin temp
	mkdir -p temp
