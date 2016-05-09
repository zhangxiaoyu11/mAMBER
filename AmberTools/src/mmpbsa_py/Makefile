# Makefile for Python reimplementation of MM/PB(GB)SA written by Jason Swails
include ../config.h

SHELL = /bin/sh

install: 
	/bin/bash setup.sh $(BINDIR) $(PYTHON)
	$(BINDIR)/nab -o mmpbsa_py_nabnmode mmpbsa_entropy.nab
	/bin/mv mmpbsa_py_nabnmode $(BINDIR)
	$(BINDIR)/nab -o mmpbsa_py_energy mmpbsa_energy.nab
	/bin/mv mmpbsa_py_energy $(BINDIR)

serial: install

parallel:
	-/bin/bash setup.sh $(BINDIR) $(PYTHON) par

docs:
	doxygen mmpbsa.Doxyfile

clean:
	-(find . -name "*.pyc" | /usr/bin/xargs /bin/rm -f)
	/bin/rm -f mmpbsa_entropy.c mmpbsa_energy.c

uninstall:
	/bin/rm -fr $(BINDIR)/MMPBSA_mods $(BINDIR)/mmpbsa_py_nabnmode $(BINDIR)/MMPBSA $(BINDIR)/MMPBSA.MPI $(BINDIR)/mmpbsa_py_energy
	/bin/rm -fr $(BINDIR)/MMPBSA.py $(BINDIR)/MMPBSA.py.MPI $(BINDIR)/ante-MMPBSA.py $(BINDIR)/mpi4py
