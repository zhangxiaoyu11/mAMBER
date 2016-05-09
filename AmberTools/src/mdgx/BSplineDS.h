#ifndef BSplineStructs
#define BSplineStructs

struct BSplineCoeff {
  int m;
  double s;
  double d;
};
typedef struct BSplineCoeff bcof;

struct BSplineMap {
  int natom;
  int maxatom;
  int* atmid;
  bcof* xcof;
  bcof* ycof;
  bcof* zcof;
};
typedef struct BSplineMap bmap;

#endif
