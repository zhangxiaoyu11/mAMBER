# Makefile for PARMED written by Jason Swails
include ../config.h

SHELL = /bin/sh

install: 
	./setup.sh $(BINDIR) $(PYTHON)

serial: install

parallel:
	./setup.sh $(BINDIR) $(PYTHON) par

docs:
	doxygen parmed.Doxyfile

clean:
	-(find . -name "*.pyc" | /usr/bin/xargs /bin/rm)

uninstall:
	-/bin/rm -fr $(BINDIR)/chemistry $(BINDIR)/ParmedTools $(BINDIR)/parmed.py $(BINDIR)/xparmed.py
