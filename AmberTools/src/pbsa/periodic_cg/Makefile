#
#************************************************************************
#                              AMBER                                   **
#                                                                      **
#               Copyright (c) 1986, 1991, 1995, 1997, 1999             **
#                Regents of the University of California               **
#                       All Rights Reserved.                           **
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
.SUFFIXES : .F90
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
#
# Makefile for UCSF Amber 8   PBSA
#
include ../config.h

# Sander Fortran source files are free format
LOCALFLAGS = $(FREEFORMAT_FLAG)

SHELL=/bin/sh
VPATH=

SRC=	\
	constants.F90 pbsa.F90 runmd.F90 runmin.F90 mdread.F90 rdparm.F90 getcor.F90 \
	force.F90 is_init.F90 sa_driver.F90 np_force.F90 pb_force.F90 \
        pb_exmol.F90 pb_dielectrics.F90 pb_mpfrc.F90 pb_direct.F90 pb_list.F90 \
	egb.F90 ene.F90 decomp.F90 icosasurf.F90 findmask.F90 printe.F90 locmem.F90 \
	fastwt.F90 echoin.F90 parallel.F90 shake.F90 mdwrit.F90 minrit.F90 \
	set.F90 setmm.F90 dynlib.F90 mdfil.F90 nmlsrc.F90 \
	rstack.F90 istack.F90 rfree.F90 rgroup.F90 \
	random.F90 amopen.F90 debug.F90 degcnt.F90 \
	mexit.F90 new_time.F90 matinv.F90 assert.F90 veclib.F90 constantph.F90

OBJ=	\
	constants.o pbsa.o runmd.o runmin.o mdread.o rdparm.o getcor.o \
	force.o is_init.o sa_driver.o np_force.o pb_force.o \
        pb_exmol.o pb_dielectrics.o pb_mpfrc.o pb_direct.o pb_list.o \
	egb.o ene.o decomp.o icosasurf.o findmask.o printe.o locmem.o \
	fastwt.o echoin.o parallel.o shake.o mdwrit.o minrit.o \
	set.o setmm.o dynlib.o mdfil.o nmlsrc.o \
	rstack.o istack.o rfree.o rgroup.o \
	random.o amopen.o debug.o degcnt.o \
	mexit.o new_time.o matinv.o assert.o veclib.o constantph.o

all::    pbsa$(SFX)

install:	pbsa$(SFX)
	/bin/mv pbsa$(SFX) ../../exe

#---------------------------------------------------------------------------
pbsa$(SFX):	$(OBJ) syslib nxtsec lapack blas
	$(LOAD) -o pbsa$(SFX) $(OBJ) ../lib/sys.a ../lib/nxtsec.o \
		$(LOADLIB) ../lapack/lapack.a ../blas/blas.a

#-----------LIBS

syslib:
	cd ../lib; make sys.a

lapack:  
	cd ../lapack; make

blas:  
	cd ../blas; make

nxtsec:  
	cd ../lib; make nxtsec.o

#-----------OBJ L3 : compile at highest optimization level

runmd.o:	runmd.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

runmin.o:	runmin.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

force.o:	force.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

ene.o:		ene.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

sa_driver.o:	sa_driver.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

np_force.o:	np_force.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

pb_force.o:	pb_force.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

pb_exmol.o:	pb_exmol.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

pb_dielectrics.o:	pb_dielectrics.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

pb_mpfrc.o:	pb_mpfrc.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

pb_direct.o:	pb_direct.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

pb_list.o:	pb_list.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

is_init.o:	is_init.F90 
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

egb.o:	egb.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

icosasurf.o:	icosasurf.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

fastwt.o:	fastwt.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

indexn.o:     	indexn.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

shake.o:	shake.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

rstack.o:	rstack.F90 rstack.h
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

istack.o:	istack.F90 istack.h
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

debug.o:	debug.F90 
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

new_time.o:	new_time.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

veclib.o:	veclib.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

random.o:	random.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

constantph.o:	constantph.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

clean:
	-/bin/rm pbsa$(SFX) 
	-/bin/rm *.o
	-/bin/rm 
	-/bin/rm *.mod
	-/bin/rm *.d

#----------HEADER DEPENDENCIES (constructed from ./makedepend)

amopen.o: \
                    copyright.h
assert.o: \
                    copyright.h\
                    assert.h
constantph.o: \
                    copyright.h\
                    dprec.h\
                    dynph.h\
                    random.h\
                    md.h\
                    files.h\
                    constants.h\
                    extra.h
debug.o: \
                    copyright.h\
                    dprec.h\
                    flocntrl.h\
                    debug.h\
                    md.h
decomp.o: \
                    copyright.h\
                    assert.h\
                    dprec.h\
                    memory.h
degcnt.o: \
                    copyright.h\
                    dprec.h
dynlib.o: \
                    copyright.h\
                    dprec.h\
                    box.h\
                    files.h\
                    md.h\
                    extra.h\
                    parallel.h
echoin.o: \
                    copyright.h
egb.o: \
                    copyright.h\
                    assert.h\
                    dprec.h\
                    icosasurf.o\
                    decomp.o\
                    parallel.h\
                    md.h\
                    parms.h\
                    def_time.h\
                    constants.h\
                    gbsa.h
ene.o: \
                    copyright.h\
                    dprec.h\
                    decomp.o\
                    parallel.h\
                    md.h\
                    box.h\
                    parms.h\
                    flocntrl.h
fastwt.o: \
                    copyright.h\
                    dprec.h\
                    parallel.h\
                    md.h\
                    constants.h
findmask.o: \
                    dprec.h
force.o: \
                    copyright.h\
                    dprec.h\
                    egb.o\
                    pb_force.o\
                    np_force.o\
                    constants.h\
                    def_time.h\
                    parallel.h\
                    md.h\
                    box.h\
                    memory.h\
                    parms.h\
                    files.h\
                    extra.h\
                    flocntrl.h
getcor.o: \
                    dprec.h\
                    copyright.h\
                    files.h
icosasurf.o: \
                    copyright.h\
                    assert.h\
                    dprec.h\
                    decomp.o
istack.o: \
                    copyright.h\
                    memory.h\
                    istack.h
locmem.o: \
                    copyright.h\
                    dprec.h\
                    box.h\
                    memory.h\
                    md.h\
                    dynph.h\
                    parallel.h
matinv.o: \
                    copyright.h\
                    dprec.h
mdfil.o: \
                    copyright.h\
                    files.h\
                    parallel.h
mdread.o: \
                    copyright.h\
                    dprec.h\
                    box.h\
                    constants.h\
                    def_time.h\
                    extra.h\
                    files.h\
                    md.h\
                    memory.h\
                    parms.h\
                    dynph.h\
                    random.h\
                    decomp.o\
                    findmask.o
mdwrit.o: \
                    copyright.h\
                    dprec.h\
                    files.h
mexit.o: \
                    copyright.h\
                    parallel.h
minrit.o: \
                    copyright.h\
                    dprec.h\
                    md.h\
                    files.h\
                    box.h
new_time.o: \
                    copyright.h\
                    dprec.h\
                    def_time.h\
                    new_time.h\
                    extra.h\
                    istack.h\
                    rstack.h\
                    parallel.h
nmlsrc.o: \
                    copyright.h
parallel.o: \
                    copyright.h\
                    dprec.h\
                    parallel.h\
                    extra.h\
                    md.h\
                    memory.h\
                    parms.h\
                    box.h\
                    files.h\
                    nonper.h\
                    flocntrl.h\
                    debug.h\
                    rstack.h\
                    istack.h\
                    new_time.h
sa_driver.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    constants.o\
                    is_def.h
np_force.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    is_def.h\
                    def_time.h\
                    pb_md.h\
                    constants.o\
		    sa_driver.o
pb_force.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    is_def.h\
                    def_time.h\
                    md.h\
                    pb_md.h\
                    box.h\
                    constants.o\
		    pb_exmol.o \
		    pb_dielectrics.o \
		    sa_driver.o 
pb_exmol.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    is_def.h\
                    constants.o
                    
pb_dielectrics.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    is_def.h\
                    constants.o
                    
pb_mpfrc.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    constants.o\
                    is_def.h
pb_direct.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    constants.o\
                    is_def.h
pb_list.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    is_def.h\
                    constants.o\
                    pb_force.o
is_init.o: \
                    copyright.h\
                    is_copyright.h\
                    dprec.h\
                    is_def.h\
                    files.h\
                    md.h\
                    pb_md.h\
                    parms.h\
                    constants.o\
                    sa_driver.o\
                    np_force.o\
                    pb_force.o
printe.o: \
                    copyright.h\
                    dprec.h\
                    decomp.o\
                    box.h\
                    extra.h\
                    files.h\
                    md.h\
                    memory.h\
                    constants.h
random.o: \
                    copyright.h\
                    dprec.h\
                    random.h
rdparm.o: \
                    copyright.h\
                    dprec.h\
                    md.h\
                    memory.h\
                    parms.h\
                    files.h\
                    box.h
rfree.o: \
                    copyright.h\
                    dprec.h
rgroup.o: \
                    copyright.h\
                    dprec.h
rstack.o: \
                    copyright.h\
                    memory.h\
                    rstack.h
runmd.o: \
                    copyright.h\
                    dprec.h\
                    assert.h\
                    decomp.o\
                    fastwt.o\
                    parallel.h\
                    files.h\
                    md.h\
                    box.h\
                    memory.h\
                    extra.h\
                    def_time.h\
                    constants.h\
                    pb_md.h\
                    random.h\
                    parms.h
runmin.o: \
                    copyright.h\
                    assert.h\
                    dprec.h\
                    fastwt.o\
                    parallel.h\
                    constants.h\
                    md.h\
                    box.h\
                    files.h\
                    memory.h\
                    extra.h\
                    pb_md.h
pbsa.o: \
                    copyright.h\
                    dprec.h\
                    assert.h\
                    files.h\
                    parallel.h\
                    egb.o\
                    decomp.o\
                    fastwt.o\
                    memory.h\
                    box.h\
                    md.h\
                    parms.h\
                    extra.h\
                    nonper.h\
                    def_time.h
set.o: \
                    copyright.h\
                    dprec.h\
                    files.h\
                    md.h\
                    parms.h\
                    memory.h\
                    parallel.h\
                    extra.h
setmm.o: \
                    copyright.h\
                    dprec.h\
                    memory.h\
                    box.h
shake.o: \
                    copyright.h\
                    dprec.h\
                    parallel.h\
                    md.h\
                    extra.h
veclib.o: \
                    copyright.h\
                    dprec.h\
                    constants.h
