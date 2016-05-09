#ifndef GridStructs
#define GridStructs

#include "MatrixDS.h"
#include "fftw3.h"

struct IntegerBook {
  int pag;
  int row;
  int col;
  int* data;
  int*** map;
};
typedef struct IntegerBook ibook;

struct FloatBook {
  int pag;
  int row;
  int col;
  double orig[3];
  dmat lvec;
  float* data;
  float*** map;
};
typedef struct FloatBook fbook;

struct DoubleBook {
  int pag;
  int row;
  int col;
  int pfft;
  double* data;
  fftw_complex* fdata;
  double*** map;
  fftw_complex*** fmap;
};
typedef struct DoubleBook dbook;

struct CharacterBook {
  int pag;
  int row;
  int col;
  char* data;
  char*** map;
};
typedef struct CharacterBook cbook;

struct CompressedChargeBook {
  ibook qdata;
  int nexcp;
  int maxexcp;
  double tq;
  int* eidx;
  double* eval;
};
typedef struct CompressedChargeBook qbook;

#endif
