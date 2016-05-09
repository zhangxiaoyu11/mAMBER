#ifndef TopologyHeadings
#define TopologyHeadings

#include "Constants.h"
#include "TopologyDS.h"
#include "CrdManipDS.h"
#include "pmeDirectDS.h"
#include "TrajectoryDS.h"

cmat Ascii2Mem(char* fname, int maxw, char* errmsg);

void OrderLJParameters(prmtop *tp);

#ifdef MPI
void GetPrmTop(prmtop *tp, trajcon *tj, int adjbnd);
#else
void GetPrmTop(prmtop *tp, int adjbnd);
#endif

prmtop CopyTopology(prmtop *tp);

prmtop InterpolateTopology(prmtop *tpA, prmtop *tpB, double lambda);

void FreeTopology(prmtop *tp);

void AdjustBondArray(int* A, int N, int P);

int FindAtom(prmtop *tp, int il, int ih, char* aname);

void PutPrmTop(prmtop *tp, char* fname, char* title);

int ProteinChiralityCheck(prmtop *tp, coord *tc, FILE *outp);

int FindDisulfides(prmtop *tp, coord *tc, FILE *outp);

int LocateResID(prmtop *tp, int ai, int llim, int hlim);

void LongRangeVDW(prmtop *tp, dircon *dcinp);

#endif
