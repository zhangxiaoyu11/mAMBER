#ifndef MatrixStructs
#define MatrixStructs

#include "fftw3.h"

struct IMatrix {
  int row;
  int col;
  int* data;
  int** map;
};
typedef struct IMatrix imat;

struct DMatrix {
  int row;
  int col;
  int pfft;
  double* data;
  double** map;
  fftw_complex* fdata;
  fftw_complex** fmap;
};
typedef struct DMatrix dmat;

struct CMatrix {
  int row;
  int col;
  char* data;
  char** map;
};
typedef struct CMatrix cmat;

#endif
