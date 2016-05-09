include ../config.h

LOCALFLAGS = $(FREEFORMAT_FLAG)

SRC = add_xray.F90 add_pdb.F90
OBJ = add_xray.o add_pdb.o

LIBS = ../lib/nxtsec.o

all: add_xray$(SFX) add_pdb$(SFX)

add_xray$(SFX): add_xray.o $(LIBS)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ add_xray.o \
		$(LDFLAGS) $(AMBERLDFLAGS) $(LIBS)

add_pdb$(SFX): add_pdb.o $(LIBS)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ add_pdb.o \
		$(LDFLAGS) $(AMBERLDFLAGS) $(LIBS)

clean:
	/bin/rm -f $(OBJ) add_xray$(SFX) add_pdb$(SFX) 

install: all
	/bin/mv add_xray$(SFX) add_pdb$(SFX) $(BINDIR)

uninstall:
	-rm -f $(BINDIR)/add_xray$(SFX)
	-rm -f $(BINDIR)/add_pdb$(SFX)

../lib/nxtsec.o: ../lib/nxtsec.F
	$(MAKE) -C ../lib nxtsec.o

add_xray.o:   add_xray.F90
	$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ add_xray.F90

add_pdb.o:   add_pdb.F90
	$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ add_pdb.F90

