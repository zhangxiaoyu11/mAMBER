/*! 
   \file gaRanNumGen.h
   \brief Inline functions for random number generation
   \author Martin Peters

   $Date: 2007/09/14 10:28:10 $
   $Revision: 1.4 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   This file is part of MTK++.

   MTK++ is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   MTK++ is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lessser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   ----------------------------------------------------------------------------
*/

#ifndef GARANNUMGEN_H
#define GARANNUMGEN_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

namespace MTKpp
{

// ============================================================
// Function : setSeed()
// ------------------------------------------------------------
/*!
   \brief Set Random Number Generation Seed
*/
// ============================================================
inline void setSeed(int seed)
{
    if (seed < 0) {
      srand(time(NULL));
    }
    else {
      srand(seed);
    }
}

// ============================================================
// Function : ranNumBetweenZeroAndOne()
// ------------------------------------------------------------
/*!
   \brief Returns a random number (double) between 0 and 1
*/
// ============================================================
inline double ranNumBetweenZeroAndOne()
{
    return ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
}

// ============================================================
// Function : ranNumBetweenZeroAndX()
// ------------------------------------------------------------
/*!
   \brief Returns a random number (double) between 0 and X
*/
// ============================================================
inline double ranNumBetweenZeroAndX(int X)
{
    return ranNumBetweenZeroAndOne() * X;
}

#endif // GARANNUMGEN_H

} // MTKpp namespace
