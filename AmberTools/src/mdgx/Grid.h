#ifndef GridHeadings
#define GridHeadings

#include "GridDS.h"
#include "MatrixDS.h"

ibook CreateIbook(int M, int N, int P);

dbook CreateDbook(int M, int N, int P, int prepFFT);

fbook CreateFbook(int M, int N, int P);

cbook CreateCbook(int M, int N, int P);

void DestroyIbook(ibook *A);

void DestroyFbook(fbook *A);

void DestroyDbook(dbook *A);

void DestroyCbook(cbook *A);

void CompressQ(qbook *cQ, dbook *Q);

void DecompressQ(dbook *Q, qbook *cQ);

void ExtractDpage(dbook *B, dmat *P, int n, int ll);

void ScribeDpage(dmat *P, dbook *B, int n);

dbook Dmat2Dbook(dmat *A);

double SumDbook(dbook *A);

void Interp2Grid(double* crd, double* val, int npts, fbook *A, dmat *U,
                 dmat *invU);
#endif
