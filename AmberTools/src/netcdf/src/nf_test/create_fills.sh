#!/bin/sh
# This shell script which creates the fill.nc file from fill.cdl.
# $Id: create_fills.sh,v 1.2 2007/11/15 21:44:49 jmongan Exp $

echo "*** Creating fills.nc."
set -e
../ncgen/ncgen -b $srcdir/fills.cdl
echo "*** SUCCESS!"
exit 0
