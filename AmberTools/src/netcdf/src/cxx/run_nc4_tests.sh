#!/bin/sh

# This shell script runs some netCDF C++ API tests for netcdf-4 and
# netcdf-4 classic formats.

# $Id: run_nc4_tests.sh,v 10.0 2008/04/15 23:23:12 case Exp $

set -e
echo ""
echo "*** Testing C++ API test program output for netCDF-4."

echo "*** dumping netcdf-4 format file to nctst_netcdf4.cdl and comparing..."
../ncdump/ncdump -n ref_nctst nctst_netcdf4.nc > nctst_netcdf4.cdl
cmp nctst_netcdf4.cdl $srcdir/ref_nctst_netcdf4.cdl

echo "*** dumping netcdf-4 classic format file to nctst_netcdf4_classic.cdl and comparing..."
../ncdump/ncdump -n ref_nctst nctst_netcdf4_classic.nc > nctst_netcdf4_classic.cdl
cmp nctst_netcdf4_classic.cdl $srcdir/ref_nctst_netcdf4_classic.cdl

echo "*** All tests of C++ API test output passed!"
exit 0
