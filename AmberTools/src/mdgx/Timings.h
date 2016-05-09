#ifndef TimingsHeaders
#define TimingsHeaders

#include "TimingsDS.h"
#include "pmeRecipDS.h"
#include "CellManipDS.h"

void InitExecon(execon *tm);

void mdgxStartTimer(execon *tm);

double mdgxStopTimer(execon *tm);

dmat GatherTimingData(execon *tm, cellgrid *CG, reccon *rcinp);

void PrintTimingData(execon *tm, reccon *rcinp, dmat *alltime, FILE *outp);

#endif
