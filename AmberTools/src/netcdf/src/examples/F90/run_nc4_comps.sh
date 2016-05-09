#!/bin/sh

# This shell script runs some test for netCDF-4. It checks that the
# HDF5 versions of the example data files match (according to ncdump)
# the netCDF-3 classic versions.

# $Id: run_nc4_comps.sh,v 10.0 2008/04/15 23:23:12 case Exp $

set -e
echo ""
echo "*** Comparing HDF5 example data files with netCDF classic format versions."
echo "*** checking simple_xy..."
../../ncdump/ncdump -n simple_xy nc4_simple_xy.nc > simple_xy.cdl
diff -w simple_xy.cdl $srcdir/../CDL/simple_xy.cdl
echo "*** checking sfc_pres_temp..."
../../ncdump/ncdump -n sfc_pres_temp nc4_sfc_pres_temp.nc > sfc_pres_temp.cdl
diff -w sfc_pres_temp.cdl $srcdir/../CDL/sfc_pres_temp.cdl
echo "*** checking pres_temp_4D..."
../../ncdump/ncdump -n pres_temp_4D nc4_pres_temp_4D.nc > pres_temp_4D.cdl
diff -w pres_temp_4D.cdl $srcdir/../CDL/pres_temp_4D.cdl

echo "*** All HDF5 example data files passed!"
exit 0
