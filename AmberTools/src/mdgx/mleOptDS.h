#ifndef mleOptStructs
#define mleOptStructs

struct BasisFunctionSet {
  int nbss;
  int* bsslim;
  int* bsscrd;
};
typedef struct BasisFunctionSet bssf;

#endif
