#!/bin/sh
# This shell script runs the cmp test on the example programs.
# $Id: do_comps.sh,v 10.0 2008/04/15 23:23:12 case Exp $

set -e
echo ""
echo "*** Testing that CXX examples produced same files as C examples."
echo "*** checking simple_xy.nc..."
cmp simple_xy.nc ../C/simple_xy.nc

echo "*** checking sfc_pres_temp.nc..."
cmp sfc_pres_temp.nc ../C/sfc_pres_temp.nc

echo "*** checking pres_temp_4D.nc..."
cmp pres_temp_4D.nc ../C/pres_temp_4D.nc

echo "*** All CXX example comparisons worked!"
exit 0
