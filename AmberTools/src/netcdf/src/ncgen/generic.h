/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /home/case/cvsroot/amber11/AmberTools/src/netcdf/src/ncgen/generic.h,v 9.1 2007/11/15 21:44:48 jmongan Exp $
 *********************************************************************/

#ifndef UD_GENERIC_H
#define UD_GENERIC_H

union generic {			/* used to hold any kind of fill_value */
    float floatv;
    double doublev;
    int intv;
    short shortv;
    char charv;
};

#endif
