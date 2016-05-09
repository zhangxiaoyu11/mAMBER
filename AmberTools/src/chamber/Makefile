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
SHELL=/bin/sh
CHAMBERSRC=	main.F90 psfprm.F90 psf_strings.F90 write_prmtop.F90 gbradii.F90 sort.F90 mkexcl.F90 cmap.F90 molnt.F90
CHAMBEROBJ=	main.o psfprm.o psf_strings.o write_prmtop.o gbradii.o sort.o mkexcl.o cmap.o molnt.o

include ../config.h

# Chamber Fortran source files are free format

LOCALFLAGS = $(FREEFORMAT_FLAG)

PROGS = chamber$(SFX)

all:	install

chamber$(SFX):	$(CHAMBEROBJ)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o chamber$(SFX) $(CHAMBEROBJ) \
		$(LDFLAGS) $(AMBERLDFLAGS)

install:
	$(MAKE) $(PROGS)
	mv $(PROGS) $(BINDIR)

clean:
	-rm -f *.o  *.mod $(PROGS)

depend:
	./makedepend > depend

uninstall:
	-for file in $(PROGS) ; do \
	rm -f $(BINDIR)/$${file} ; \
	done

#----------HEADER DEPENDENCIES (constructed from ./makedepend > depend)
include depend
# DO NOT DELETE
