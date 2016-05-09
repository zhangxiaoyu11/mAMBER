#ifndef mleRecipStructs
#define mleRecipStructs

#include "Matrix.h"

struct GridToGridMap {
  int ng;
  int ggordr;
  dmat s;
  imat m;
};
typedef struct GridToGridMap g2gmap;

#endif
