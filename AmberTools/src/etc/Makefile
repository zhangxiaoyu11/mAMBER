#
#************************************************************************
#                             AMBER                                    **
#                                                                      **
#                Copyright (c) 1986, 1991, 1995                        **
#             Regents of the University of California                  **
#                       All Rights Reserved.                           **
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
#

.SUFFIXES : .c .F .F90 .o

SHELL=/bin/sh
NSRC= 	nucgen.F90 ngfil.F90
NOBJ=	nucgen.o ngfil.o
ASRC=   ambpdb.F90 forFIRSTtype.F90  forFIRSTbond.F90 \
                   forFIRSThbond.F90 forFIRSTteth.F90 \
		   forFIRSToutp.F90  forFIRSTsugar.F90 calcpka.F90
AOBJ=   ambpdb.o forFIRSTtype.o forFIRSTbond.o forFIRSThbond.o \
                 forFIRSTteth.o forFIRSToutp.o forFIRSTsugar.o
include ../config.h

PROGS= ambpdb$(SFX) elsize$(SFX) molsurf$(SFX) add_pdb$(SFX) resp$(SFX) calcpka$(SFX)

PERLPROGS= lmodprmtop process_mdout.perl process_minout.perl mdout2pymbar.pl

SHELLPROGS= charmmlipid2amber.x sgldinfo.sh sgldwt.sh

all:		$(PROGS) python

nukit$(SFX):	nukit.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o nukit nukit.o $(LDFLAGS) $(AMBERLDFLAGS)

nucgen$(SFX):	nucgen_modules.o $(NOBJ) syslib
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o nucgen$(SFX) \
		nucgen_modules.o $(NOBJ) ../lib/sys2.a $(LDFLAGS) $(AMBERLDFLAGS)

ambpdb$(SFX): 	$(AOBJ) syslib
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o ambpdb$(SFX) \
		$(AOBJ) ../lib/sys2.a $(LDFLAGS) $(AMBERLDFLAGS)

add_pdb$(SFX): add_pdb.o syslib
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o add_pdb$(SFX) \
		add_pdb.o ../lib/sys2.a $(LDFLAGS) $(AMBERLDFLAGS)

calcpka$(SFX):
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o calcpka$(SFX) calcpka.F90

elsize$(SFX): elsize.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) -o elsize$(SFX) elsize.o $(LM)

molsurf$(SFX):  molsurf.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) -o molsurf$(SFX) molsurf.o $(LM)
molsurf.o:  molsurf.h

resp$(SFX):	lapack.o resp.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
		lapack.o resp.o -o resp$(SFX)

syslib::
	cd ../lib && $(MAKE) sys2.a

install:	all
	/bin/mv $(PROGS) $(BINDIR)
	/bin/cp $(PERLPROGS) $(BINDIR)
	/bin/cp $(SHELLPROGS) $(BINDIR)

parallel:
	$(FC) $(FPPFLAGS) -o numprocs numprocs.F90
	/bin/cp numprocs ../../test
	/bin/cp numprocs ../../../test

python:
	./py_setup.sh $(BINDIR) $(PYTHON)

clean:
	/bin/rm -f *.o 90 *.mod  $(PROGS)

uninstall:
	-(for file in $(PROGS) $(PERLPROGS) $(SHELLPROGS) ; do \
	/bin/rm -f $(BINDIR)/$${file} ; \
	done)
	-/bin/rm -fr $(BINDIR)/softcore_setup.py $(BINDIR)/cpinutil.py $(BINDIR)/cpinutils

#  default rules for Fortran and C compilation:

.F90.o: $<
	$(FC) $(FPPFLAGS) -c $(FREEFORMAT_FLAG) $(FNOOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

# lapack.F still uses fixed-format F77
.F.o:   $<
	$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

.c.o:
	$(CC) -c $(CPPFLAGS) $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

include depend
