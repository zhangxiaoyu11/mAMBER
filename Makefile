#  Simple top-level Makefile to point users to those hidden below:

include config.h
SHELL=/bin/sh

install:
	cd AmberTools/src && $(MAKE) install
	@(if [ -n "$(BUILDAMBER)" ] ; then \
	    cd src && $(MAKE) install; \
	  else \
		echo "==============================================================" ;\
		echo "$(AMBERHOME)/src/Makefile not found, or -noamber was set." ;\
		echo "This is expected if you have not installed Amber12." ;\
		echo "==============================================================" ;\
	fi ;\
	)

clean:
	-(cd AmberTools/src && $(MAKE) clean)
	-(cd src && $(MAKE) clean)

uninstall:
	-(cd AmberTools/src && $(MAKE) uninstall)
	-(cd src && $(MAKE) uninstall)

clean.test:
	-(cd AmberTools/test && $(MAKE) clean)
	-(cd test && $(MAKE) clean)

test::  test.$(INSTALLTYPE)

test.serial:
	-(cd AmberTools/test && $(MAKE) test)
	-@(if [ -n "$(BUILDAMBER)" ] ; then \
	    cd test && $(MAKE) test; \
	    echo "" ; \
	    echo "Summary of AmberTools serial tests:" ; \
	    echo "" ; cat ../logs/test_at_serial/at_summary; \
	  else \
		echo "==============================================================" ;\
		echo "$(AMBERHOME)/src/Makefile not found, or -noamber was set." ;\
		echo "This is expected if you have not installed Amber12." ;\
		echo "==============================================================" ;\
	fi ;\
	)

test.parallel:
	-(cd AmberTools/test && $(MAKE) test.parallel)
	-@(if [ -f test/Makefile ] ; then \
	    cd test && $(MAKE) test.parallel; \
	    echo "" ; \
	    echo "Summary of AmberTools parallel tests:" ; \
	    echo "" ; cat ../logs/test_at_parallel/at_summary; \
	  else \
		echo "==============================================================" ;\
		echo "$(AMBERHOME)/test/Makefile not found." ;\
		echo "This is expected if you have not installed Amber12." ;\
		echo "==============================================================" ;\
	fi ;\
	)

test.cuda:
	-(cd test && $(MAKE) test.cuda)

test.cuda_parallel:
	-(cd test && $(MAKE) test.cuda.parallel)

nabonly:
	-(cd AmberTools/src && $(MAKE) nabonly)

openmp:
	-(cd AmberTools/src && $(MAKE) openmp)
