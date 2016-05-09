#   top-level $(MAKE)file for the AmberTools suite of programs
AMBERTOOLS=AmberTools13
CONFIG_FILE=config.h
CONFIG_COMMAND=./configure
PARALLEL_TOKEN=DMPI

include $(CONFIG_FILE)

.NOTPARALLEL:configured_serial THIRDPARTY $(MTKPP)
install: $(INSTALLTYPE)
.PHONY: INSTALL_MTKPP

serial: configured_serial THIRDPARTY $(MTKPP)
	@echo "Starting installation of ${AMBERTOOLS} serial at `date`".
#	utility routines and libraries:
	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
	(cd byacc && $(MAKE) install )
	(cd arpack && $(MAKE) install );
	(cd lapack && $(MAKE) $(LAPACK) )
	(cd blas && $(MAKE) $(BLAS) )
	(cd c9x-complex && $(MAKE) $(C9XCOMPLEX) )
	(cd etc && $(MAKE) install )
	(cd chamber && $(MAKE) install )
	(cd pbsa && $(MAKE) install )
	(cd cifparse && $(MAKE) install )

#	antechamber:
	(cd antechamber && $(MAKE) install )
	(cd sqm && $(MAKE) install )

#   ptraj:
	(cd ptraj && $(MAKE) install )

#   cpptraj:
	(cd cpptraj && $(MAKE) install)

#   miscellaneous:
	(cd reduce && $(MAKE) install )

#   leap
	(cd leap && $(MAKE) install )

#   nab:
	(cd sff && $(MAKE) install )
	(cd pbsa && $(MAKE) libinstall )
	(cd rism && $(MAKE) $(RISM) )
	(cd nab && $(MAKE) install )
	(cd nss && $(MAKE) install )

#   rism
	(if [ $(RISM) = "yes" ]; then \
		cd rism && $(MAKE) rism3d.snglpnt$(SFX); \
	fi;\
	)

#   mdgx:
	(cd mdgx && $(MAKE) $(MDGX) )

#   xtalutil and xray:
	(cd xtalutil/Clibs/ && $(MAKE) install )
	(cd xtalutil/CPrograms/ && $(MAKE) install )
	(cd xray && $(MAKE) install )

#	Parmed
	(cd parmed && $(MAKE) install)

#	mmpbsa_py
	(cd mmpbsa_py && $(MAKE) install)

#	amberlite
	(cd amberlite && $(MAKE) install)

#	paramfit
	(cd paramfit && $(MAKE) install)

	@echo ""
	@echo "Installation of ${AMBERTOOLS} serial is complete at `date`."
	@echo ""

nabonly: configured_serial $(NETCDF) $(XBLAS)  
#	utility routines and libraries:
	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
	(cd byacc && $(MAKE) install )
	(cd arpack && $(MAKE) install );
	(cd lapack && $(MAKE) $(LAPACK) )
	(cd blas && $(MAKE) $(BLAS) )
	(cd c9x-complex && $(MAKE) $(C9XCOMPLEX) )
	(cd pbsa && $(MAKE) libinstall )
	(cd rism && $(MAKE) $(RISM) )
	(cd cifparse && $(MAKE) install )
	(cd sff && $(MAKE) install )
	(cd nab && $(MAKE) install )
	(cd nss && $(MAKE) install )

	@echo
	@echo "Completed installation of NAB"
	@echo

openmp: nabonly
	(cd cpptraj && $(MAKE) install )
	(cd paramfit && $(MAKE) install )

mtkpponly: configured_serial install_mtkpp

install_mtkpp:
	@echo "Starting installation of MTK++ at `date`".
	(cd mtkpp && $(MAKE) && $(MAKE) -j 1 install )
	(cd mtkpp/tests && $(MAKE) create )

	@echo
	@echo "Completed installation of MTK++"
	@echo

.NOTPARALLEL:configured_parallel THIRDPARTY

parallel: configured_parallel THIRDPARTY
	@echo "Starting installation of ${AMBERTOOLS} MPI at `date`".
	(cd etc && $(MAKE) parallel)
	(cd sff && $(MAKE) libsff_mpi )
#	(cd pbsa && $(MAKE) clean && $(MAKE) install.parallel )
	(cd pbsa && $(MAKE) libinstall )
	(cd byacc && $(MAKE) install )
	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
	(cd mdgx && $(MAKE) $(MDGX) )
	(if [ "$(RISM)" = "yes" ]; then \
		cd rism && $(MAKE) install_mpi ; \
	fi; \
	)
	(cd nab && $(MAKE) install )
	(cd parmed && $(MAKE) parallel )
	(cd mmpbsa_py && $(MAKE) parallel )
	(if [ "$(RISM)" = "yes" ]; then \
		cd rism && $(MAKE) rism3d.snglpnt.MPI$(SFX) ; \
	fi; \
	)

	@echo ""
	@echo "Installation of ${AMBERTOOLS} MPI is complete at `date`."
	@echo ""

configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
		echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE) !" ;\
		exit 2 ;\ # $(CONFIG_COMMAND) ;\
	fi ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type parallel, not serial !" ;\
		echo "  Rerun $(CONFIG_COMMAND) without the -mpi flag." ;\
		exit 2 ;\
	fi ;\
	)

configured_parallel: configured
	@(grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null || \
		{ echo "Error: $(CONFIG_FILE) is of type serial, not parallel !" ;\
			echo "  Rerun $(CONFIG_COMMAND) with the -mpi flag." ;\
			exit 2 ;\
		} ;\
	)

clean::  netcdf_clean
	-(cd ucpp-1.3 && $(MAKE) clean )
	-(cd byacc && $(MAKE) clean )
	-(cd cifparse && $(MAKE) clean )
	-(cd nab && $(MAKE) clean )
	-(cd sff && $(MAKE) clean )
	-(cd nss && $(MAKE) clean )
	-(cd leap && $(MAKE) clean )
	-(cd arpack && $(MAKE) clean )
	-(cd blas && $(MAKE) clean )
	-(cd xblas && $(MAKE) clean )
	-(cd lapack && $(MAKE) clean )
	-(cd c9x-complex && $(MAKE) clean )
	-(cd etc && $(MAKE) clean )
	-(cd chamber && $(MAKE) clean )
	-(cd pbsa && $(MAKE) clean )
	-(cd antechamber && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )
	-(cd lib && $(MAKE) clean )
	-(cd ptraj && $(MAKE) clean )
	-(cd cpptraj && $(MAKE) clean)
	-(cd reduce && $(MAKE) clean )
	-(cd mtkpp && $(MAKE) clean )
	-(cd mdgx && $(MAKE) clean )
	-(cd xtalutil/CPrograms && $(MAKE) clean )
	-(cd xray && $(MAKE) clean )
	-(cd paramfit && $(MAKE) clean )
	-(cd rism && $(MAKE) clean )
	-(cd fftw-3.3 && $(MAKE) clean )
	-($(RM) -f fftw-3.3/mpi/fftw3-mpi.f03 )
	-(cd mmpbsa_py && $(MAKE) clean )
	-(cd parmed && $(MAKE) clean )

netcdf_clean:
	-(cd netcdf/src && $(MAKE) clean)
	-(cd netcdf/lib && rm -f *.a *.la)
	-(cd netcdf/include && rm -f *.mod netcdf.*)
	-(cd netcdf/bin && rm -f *)

uninstall:
	-(cd ucpp-1.3 && $(MAKE) uninstall )
	-(cd byacc && $(MAKE) uninstall )
	-(cd cifparse && $(MAKE) uninstall )
	-(cd nab && $(MAKE) uninstall )
	-(cd sff && $(MAKE) uninstall )
	-(cd nss && $(MAKE) uninstall )
	-(cd leap && $(MAKE) uninstall )
	-(cd arpack && $(MAKE) uninstall )
	-(cd blas && $(MAKE) uninstall )
	-($(RM) -f $(LIBDIR)/libxblas.a )
	-(cd lapack && $(MAKE) uninstall )
	-(cd c9x-complex && $(MAKE) uninstall )
	-(cd etc && $(MAKE) uninstall )
	-(cd chamber && $(MAKE) uninstall )
	-(cd pbsa && $(MAKE) uninstall )
	-(cd antechamber && $(MAKE) uninstall )
	-(cd sqm && $(MAKE) uninstall )
	-(cd lib && $(MAKE) uninstall )
	-(cd ptraj && $(MAKE) uninstall )
	-(cd reduce && $(MAKE) uninstall )
	-(cd mdgx && $(MAKE) uninstall )
	-(cd xtalutil/CPrograms && $(MAKE) uninstall )
	-(cd xray && $(MAKE) uninstall )
	-(cd mtkpp && $(MAKE) uninstall )
	-(cd paramfit && $(MAKE) uninstall )
	-(cd rism && $(MAKE) uninstall )
	-(cd netcdf/src && $(MAKE) uninstall)
	-(cd fftw-3.3 && $(MAKE) uninstall)
	-(cd cpptraj && $(MAKE) uninstall)
	-(cd mmpbsa_py && $(MAKE) uninstall )
	-(cd parmed && $(MAKE) uninstall )
	-(cd amberlite && $(MAKE) uninstall )
	-/bin/rm -f ../test/numprocs
	-/bin/rm -fr $(BINDIR)/ncmpidump $(BINDIR)/ncmpigen $(BINDIR)/ncvalid $(LIBDIR)/libnetcdf.a $(INCDIR)/mtkpp
	-(cd $(LIBDIR) && /bin/rm -f libxblas-amb.a libfftw3_mpi.a libfftw3_mpi.la)
	-(cd $(INCDIR) && /bin/rm -f fftw3-mpi.f03 fftw3-mpi.h)

$(INCDIR)/netcdf.mod:
	cd netcdf/src && $(MAKE) && $(MAKE) -j 1 install

$(LIBDIR)/libxblas-amb.a:
	(cd xblas && $(MAKE) -j 1 lib-amb && mv libxblas-amb.a $(LIBDIR) )

$(LIBDIR)/libfftw3.a $(LIBDIR)/libfftw3_mpi.a: 
	cd fftw-3.3 && $(MAKE) && $(MAKE) -j 1 install;

THIRDPARTY: $(NETCDF) $(XBLAS) $(FFTW3)

cuda:
	@echo "$(AMBERTOOLS) has no CUDA-enabled components"

	-(cd ../../src && $(MAKE) cuda )

cuda_parallel:
	@echo "$(AMBERTOOLS) has no CUDA-enabled components"

	-(cd ../../src && $(MAKE) cuda_parallel )
