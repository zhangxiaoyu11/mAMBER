
include ../config.h

yes:  install

no:
	@echo "Skipping installation of mdgx"

noparallel:
	@echo "Skipping installation of mdgx"

install: $(BINDIR)/mdgx$(SFX) 

parallel: $(BINDIR)/mdgx.MPI$(SFX)

uninstall:
	-/bin/rm -f $(BINDIR)/mdgx$(SFX)
	-/bin/rm -f $(BINDIR)/mdgx.MPI$(SFX)

MDGX_CFILES = \
	    CompFrc.c \
	    BSpline.c \
	    Grid.c \
	    Random.c \
	    mdgxVector.c \
	    Matrix.c \
	    ChargeMap.c \
	    pmeRecip.c \
	    mleRecip.c \
	    pmeDirect.c \
	    CrdManip.c \
	    CellManip.c \
	    Topology.c \
	    Trajectory.c \
	    SpecialMath.c \
	    Nonbonded.c \
	    Bonded.c \
	    Parse.c \
	    Command.c \
	    Constraints.c \
	    Thermostats.c \
	    Barostats.c \
	    Integrator.c \
	    Timings.c \
	    Manual.c \
	    VirtualSites.c \
	    Buckingham.c \
	    ThermoDyn.c \
	    ChargeFit.c \
	    ParamOut.c \
	    ParamRead.c \
	    ParamFit.c \
	    Debug.c \
	    BroadcastCommand.c \
	    MPIMap.c \
	    MPITypeCast.c \
	    ptrajmask.c \
	    ../sff/AmberNetcdf.c \
	    mdgx.c

MDGX_OBJS = \
	    CompFrc.o \
	    BSpline.o \
	    Grid.o \
	    Random.o \
	    mdgxVector.o \
	    Matrix.o \
	    ChargeMap.o \
	    pmeRecip.o \
	    mleRecip.o \
	    pmeDirect.o \
	    CrdManip.o \
	    CellManip.o \
	    Topology.o \
	    Trajectory.o \
	    SpecialMath.o \
	    Nonbonded.o \
	    Bonded.o \
	    Parse.o \
	    Command.o \
	    Constraints.o \
	    Thermostats.o \
	    Barostats.o \
	    Integrator.o \
	    Timings.o \
	    Manual.o \
	    VirtualSites.o \
	    Buckingham.o \
	    ThermoDyn.o \
	    ChargeFit.o \
	    ParamOut.o \
	    ParamRead.o \
	    ParamFit.o \
	    Debug.o \
	    BroadcastCommand.o \
	    MPIMap.o \
	    MPITypeCast.o \
	    ptrajmask.o \
	    ../sff/AmberNetcdf.o \
	    mdgx.o

MDGX_HEADERS = \
	       CompFrc.h \
	       BSpline.h \
	       Grid.h \
	       Random.h \
	       mdgxVector.h \
	       Matrix.h \
	       ChargeMap.h \
	       pmeRecip.h \
	       mleRecip.h \
	       pmeDirect.h \
	       CrdManip.h \
	       CellManip.h \
	       Topology.h \
	       Trajectory.h \
	       SpecialMath.h \
	       Nonbonded.h \
	       Bonded.h \
	       Parse.h \
	       Command.h \
	       Constraints.h \
	       Thermostats.h \
	       Barostats.h \
	       Integrator.h \
	       Timings.h \
	       Manual.h \
	       VirtualSites.h \
	       Buckingham.h \
	       ThermoDyn.h \
	       ChargeFit.h \
	       ParamOut.h \
	       ParamRead.h \
	       ParamFit.h \
	       Debug.h \
	       BroadcastCommand.h \
	       MPIMap.h \
	       MPITypeCast.h \
	       ptrajmask.h \
	       ../sff/AmberNetcdf.h 

$(BINDIR)/mdgx$(SFX) : $(NETCDF) $(LIBDIR)/libfftw3.a $(MDGX_OBJS)
	$(CC) $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) \
	-o $@ $(MDGX_OBJS) $(LIBDIR)/libfftw3.a $(NETCDFLIB) $(LM)

$(BINDIR)/mdgx.MPI$(SFX) : $(NETCDF) $(LIBDIR)/libfftw3.a $(MDGX_OBJS)
	$(CC) $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) \
	-o $@ $(MDGX_OBJS) $(LIBDIR)/libfftw3.a $(NETCDFLIB) $(LM)

$(LIBDIR)/libfftw3.a:
	cd ../fftw-3.3 && $(MAKE) && $(MAKE) -j 1 install

$(INCDIR)/netcdf.mod:
	cd ../netcdf/src && $(MAKE) install

.c.o:
	$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -I$(INCDIR) -o $@ $<

clean:
	/bin/rm -f $(MDGX_OBJS) mdgx$(SFX)

