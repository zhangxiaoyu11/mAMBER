#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "Constants.h"
#include "Parse.h"
#include "Macros.h"
#include "mdgxVector.h"
#include "Matrix.h"
#include "pmeDirect.h"
#include "Random.h"
#include "Command.h"

#include "pmeDirectDS.h"
#include "pmeRecipDS.h"
#include "TrajectoryDS.h"
#include "TopologyDS.h"
#include "CrdManipDS.h"
#include "ChargeFitDS.h"
#include "ParamFitDS.h"

/***=======================================================================***/
/*** PreExistingFileNames: for certain output trajectories, it may be      ***/
/***                       possible to specify multiple file base names    ***/
/***                       even while the first may be given by a default. ***/
/***                       This routine automates checks for the existence ***/
/***                       of all such files.                              ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   C:       matrix of pre-specified file names                         ***/
/***   def0:    the default for the first file name                        ***/
/***   fspec:   integer array filled with ones if file names have been     ***/
/***            specified or changed from the defaults, or zeros otherwise ***/
/***=======================================================================***/
static void PreExistingFileNames(cmat *C, char* def0, int* fspec)
{
  int i;

  /*** Test the existence of the first file ***/
  fspec[0] = (strcmp(C->map[0], def0) == 0) ? 0 : 1;
  for (i = 1; i < C->row; i++) {
    fspec[i] = (C->map[i][0] == '\0') ? 0 : 1;
  }
  for (i = C->row; i < MAXSYS; i++) {
    fspec[i] = 0;
  }
}

/***=======================================================================***/
/*** GetDataFileNames: this function reads in the names of all necessary   ***/
/***                   data (input and output) for a molecular dynamics    ***/
/***                   run.  This routine also seeks "base" names and the  ***/
/***                   corresponding suffixes, so mdgx can write numerous  ***/
/***                   output files with the same base and suffix but      ***/
/***                   different numbers (i.e. trj101.crd, where the base  ***/
/***                   is "trj", the suffix is ".crd", and the number is   ***/
/***                   101).  Note that the initialization of file names   ***/
/***                   is done in the CommandLineControl function.         ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tp:      the topology structure(s) (each topology stores its own    ***/
/***            source file name); tp is an array with two elements for up ***/
/***            to two systems                                             ***/
/***   tj:      the trajectory structure (stores the restart file name, as ***/
/***            well as many other output file names)                      ***/
/***   inp:     the input file                                             ***/
/***=======================================================================***/
static void GetDataFileNames(prmtop* tp, trajcon *tj, FILE *inp)
{
  int i, collect, EXoutbase, EXforcedump, EXparmfile, EXfmodfile; 
  int* EXinpcrd;
  int* EXfrcbase;
  int* EXrstbase;
  int* EXtrjbase;
  int* EXvelbase;
  int* EXfrcsuff;
  int* EXrstsuff;
  int* EXtrjsuff;
  int* EXvelsuff;
  int EXtp[2], EXeprules[2];
  char line[MAXLINE], searchtag[32], searchalias[32];
  cmat lwords;

  /*** Allocate memory for pre-existing file name/suffix flags ***/
  EXinpcrd = (int*)calloc(MAXSYS, sizeof(int));
  EXfrcbase = (int*)calloc(MAXSYS, sizeof(int));
  EXrstbase= (int*)calloc(MAXSYS, sizeof(int));
  EXtrjbase = (int*)calloc(MAXSYS, sizeof(int));
  EXvelbase = (int*)calloc(MAXSYS, sizeof(int));
  EXfrcsuff = (int*)calloc(MAXSYS, sizeof(int));
  EXrstsuff = (int*)calloc(MAXSYS, sizeof(int));
  EXtrjsuff = (int*)calloc(MAXSYS, sizeof(int));
  EXvelsuff = (int*)calloc(MAXSYS, sizeof(int));

  /*** First, verify that certain variables which could have    ***/
  /*** been given on the command line are not already specified ***/
  EXoutbase = (strcmp(tj->outbase, "mdout") != 0) ? 1 : 0;
  EXforcedump = (strcmp(tj->dumpname, "forcedump.dat") != 0) ? 1 : 0;
  EXparmfile = (tj->parmfile[0] == '\0') ? 0 : 1;
  EXfmodfile = (tj->fmodfile[0] == '\0') ? 0 : 1;
  EXtp[0] = (strcmp(tp[0].source, "prmtop") != 0) ? 1 : 0;
  EXtp[1] = (tp[1].source[0] == '\0') ? 0 : 1;
  for (i = 0; i < 2; i++) {
    EXeprules[i] = (tp[i].eprulesource[0] != '\0') ? 1 : 0;
  }
  PreExistingFileNames(&tj->ipcname, "inpcrd", EXinpcrd);
  PreExistingFileNames(&tj->trjbase, "mdcrd", EXtrjbase);
  PreExistingFileNames(&tj->velbase, "mdvel", EXvelbase);
  PreExistingFileNames(&tj->frcbase, "mdfrc", EXfrcbase);
  PreExistingFileNames(&tj->rstbase, "restrt", EXrstbase);

  /*** Search through the &files namelist ***/
  collect = AdvanceToSegment(inp, "files", 1);
  while (collect == 1) {

    /*** Read the next line ***/
    fgets(line, MAXLINE, inp);
    RemoveWhiteSpace(line, MAXLINE);
    RemoveComments(line);

    /*** Break if the line is "&end" ***/
    collect = DetectNamelistEnd(line, "GetDataFileNames");
    if (collect == 0) {
      continue;
    }

    /*** Eliminate and add spaces between special ***/
    /*** characters "=", "\n", and ","            ***/
    NixCommaCarriage(line);
    EqualSpace(line);
    lwords = ParseWords(line);

    /*** List of single string arguments to search for ***/
    if (EXtp[0] == 0) {
      SeekString(lwords, tp[0].source, "Topology", "-p");
    }
    if (EXeprules[0] == 0) {
      SeekString(lwords, tp[0].eprulesource, "EPRules", "-xpt");
    }
    for (i = 0; i < 2; i++) {
      sprintf(searchtag, "Topology%d", i+1);
      sprintf(searchalias, "-p%d", i+1);
      if (EXtp[i] == 0) {
	SeekString(lwords, tp[i].source, searchtag, searchalias);
      }
      sprintf(searchtag, "EPRules%d", i+1);
      sprintf(searchalias, "-xpt%d", i+1);
      if (EXeprules[i] == 0) {
	SeekString(lwords, tp[i].eprulesource, searchtag, searchalias);
      }
    }
    if (EXoutbase == 0) {
      SeekString(lwords, tj->outbase, "OutputBase", "-o");
    }
    if (EXforcedump == 0) {
      SeekString(lwords, tj->dumpname, "ForceDump", "-d");
    }
    if (EXparmfile == 0) {
      SeekString(lwords, tj->parmfile, "ParmFile", "-parm");
    }
    if (EXfmodfile == 0) {
      SeekString(lwords, tj->fmodfile, "FrcmodFile", "-fmod");
    }
    if (EXinpcrd[0] == 0) {
      SeekString(lwords, tj->ipcname.map[0], "StartCoord", "-c");
    }
    if (EXtrjbase[0] == 0) {
      SeekString(lwords, tj->trjbase.map[0], "CrdTrajBase", "-x");
    }
    if (EXvelbase[0] == 0) {
      SeekString(lwords, tj->velbase.map[0], "VelTrajBase", "-v");
    }
    if (EXfrcbase[0] == 0) {
      SeekString(lwords, tj->frcbase.map[0], "FrcTrajBase", "-f");
    }
    if (EXrstbase[0] == 0) {
      SeekString(lwords, tj->rstbase.map[0], "Restart", "-r");
    }
    SeekString(lwords, tj->outsuff, "OutputSuff", "-osf");
    SeekString(lwords, tj->rsrptname, "ResReport", "-rrp");
    SeekString(lwords, tj->trjsuff.map[0], "CrdTrajSuff", "-xsf");
    SeekString(lwords, tj->velsuff.map[0], "VelTrajSuff", "-vsf");
    SeekString(lwords, tj->frcsuff.map[0], "FrcTrajSuff", "-fsf");
    SeekString(lwords, tj->rstsuff.map[0], "RestartSuff", "-rsf");

    /*** List of group string arguments to search for ***/
    tj->ipcname = SeekNString(lwords, &tj->ipcname, EXinpcrd, "StartCoord",
			      "-c");
    tj->trjbase = SeekNString(lwords, &tj->trjbase, EXtrjbase, "CrdTrajBase",
			      "-x");
    tj->velbase = SeekNString(lwords, &tj->velbase, EXvelbase, "VelTrajBase",
			      "-v");
    tj->frcbase = SeekNString(lwords, &tj->frcbase, EXfrcbase, "FrcTrajBase",
			      "-f");
    tj->rstbase = SeekNString(lwords, &tj->rstbase, EXrstbase, "Restart",
			      "-r");
    tj->trjsuff = SeekNString(lwords, &tj->trjsuff, EXtrjsuff, "CrdTrajSuff",
			      "-xsf");
    tj->velsuff = SeekNString(lwords, &tj->velsuff, EXvelsuff, "VelTrajSuff",
			      "-vsf");
    tj->frcsuff = SeekNString(lwords, &tj->frcsuff, EXfrcsuff, "FrcTrajSuff",
			      "-fsf");
    tj->rstsuff = SeekNString(lwords, &tj->rstsuff, EXrstsuff, "Restart",
			      "-rsf");

    /*** Free allocated memory ***/
    DestroyCmat(&lwords);
  }

  /*** Free allocated memory ***/
  free(EXinpcrd);
  free(EXrstbase);
  free(EXtrjbase);
  free(EXvelbase);
  free(EXfrcbase);
  free(EXrstsuff);
  free(EXtrjsuff);
  free(EXvelsuff);
  free(EXfrcsuff);
}

/***=======================================================================***/
/*** CountOutputStreams: this function counts the number of output streams ***/
/***                     specified in a particular character matrix.  Also ***/
/***                     checks to make sure that there are no gaps in the ***/
/***                     list of file names.                               ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   T:       the character matrix containing file (base) names          ***/
/***=======================================================================***/
static int CountOutputStreams(cmat *T)
{
  int i, nout, nhigh;

  if (T->row == 0) {
    return 0;
  }

  nout = 0;
  for (i = 0; i < T->row; i++) {
    if (T->map[i][0] != '\0') {
      nout++;
      nhigh = i;
    }
  }
  if (nout != nhigh+1) {
    printf("CountOutputStreams >> Error.  List of output names invalid:\n");
    for (i = 0; i < nhigh+1; i++) {
      printf("CountOutputStreams >>  %s\n", T->map[i]);
    }
    exit(1);
  }

  return nout;
}

/***=======================================================================***/
/*** ExpandSuffixes: this function fills in blank suffix names when more   ***/
/***                 base file names have been specified.  If a particular ***/
/***                 suffix name is specified, all blank suffixes after it ***/
/***                 will be assigned the same string until a non-blank    ***/
/***                 suffix is encountered, at which point that becomes    ***/
/***                 the new replacement for blank suffixes.               ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   S:     the character matrix of suffixes                             ***/
/***   B:     the character matrix of base names                           ***/
/***=======================================================================***/
static void ExpandSuffixes(cmat *S, cmat *B)
{
  int i;
  char *ctmp;

  /*** Expand the suffixes to be as numerous as the bases ***/
  if (S->row < B->row) {
    *S = ReallocCmat(S, B->row, S->col);
  }

  /*** The first suffix must be specified or else ***/
  /*** we just wipe all the suffixes and bail out ***/
  if (S->map[0][0] == '\0') {
    for (i = 0; i < S->row*S->col; i++) {
      S->data[i] = '\0';
    }
    return;
  }
  for (i = 0; i < S->row; i++) {
    if (S->map[i][0] != '\0') {
      ctmp = S->map[i];
    }
    else {
      strcpy(S->map[i], ctmp);
    }
  }
}

/***=======================================================================***/
/*** SetOutputFrequency: this function ensures that the frequency of each  ***/
/***                     output stream is a factor of the total number of  ***/
/***                     steps in either the segment or full trajectory.   ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tj:     trajectory control information                              ***/
/***   ival:   the output frequency (also returned, possibly modified)     ***/
/***   vname:  the name of the variable being examined (printed in any     ***/
/***           warnings)                                                   ***/
/***=======================================================================***/
static int SetOutputFrequency(trajcon *tj, int ival, char* vname)
{
  int i;
  long long int stepcount, ficount;

  /*** Bail right out if the interval is just zero ***/
  if (ival == 0) {
    return ival;
  }

  /*** Make sure this interval is a clean one ***/
  stepcount = (tj->nfistep > 0) ? tj->nfistep : tj->nstep;
  if (stepcount % ival > 0) {
    printf("GetCntrlNamelist >> Warning.  When the run output is split into "
           "segments,\nGetCntrlNamelist >> all output frequencies must be "
	   "factors of nfistep.");
    if ((double)stepcount / (double)ival < 10.0) {
      for (i = ival+1; i <= stepcount; i++) {
	if (tj->nfistep % i == 0) {
	  ival = i;
	  break;
	}
      }
      printf("%s set to %d.\n", vname, ival);
    }
    else {
      if (tj->nfistep > 0) {
	ficount = tj->nstep / tj->nfistep;
      }
      stepcount = (ival+1) * (stepcount / ival);
      if (tj->nfistep > 0) {
	tj->nfistep = stepcount;
	tj->nstep = stepcount*ficount;
	printf("nfistep set to %lld,\nGetCntrlNamelist >> ", tj->nfistep);
	printf("nstep   set to %lld.\n", tj->nstep);
      }
      else {
	tj->nstep = stepcount;
	printf("nstep set to %lld.\n", stepcount);
      }
    }
  }

  return ival;
}

/***=======================================================================***/
/*** GetCntrlNamelist: this function reads in run parameters, as would be  ***/
/***                   found in the &cntrl namelist in PMEMD or SANDER     ***/
/***                   command input files.                                ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   dcinp:   direct space command information                           ***/
/***   rcinp:   reciprocal space command information                       ***/
/***   tp:      the topology structure array                               ***/
/***   tj:      the trajectory structure (stores the restart file name, as ***/
/***            well as many other output file names)                      ***/
/***   inp:     the input file                                             ***/
/***=======================================================================***/
static void GetCntrlNamelist(dircon *dcinp, prmtop* tp, trajcon *tj, FILE *inp)
{
  int i, collect, EXigseed, EXlambda, ntc;
  double EVcut, lampow, lampowm1, mfac, mxA, dmxA;
  double* pval;
  char line[MAXLINE];
  cmat lwords;

  /*** Default settings for dynamics / Hamiltonian ***/
  tp->lj14fac = 2.0;
  tp->elec14fac = 1.2;
  dcinp->Ecut = 8.0;
  dcinp->Vcut = 8.0;
  dcinp->MaxDens = 2.5;
  EVcut = -1.0;
  sprintf(tp->WaterName, "WAT ");
  tj->mode = 0;
  tj->starttime = 0.0;
  tj->dt = 0.001;
  tp->rattle = 0;
  tp->settle = 0;

  /*** Default settings for thermostat and barostat ***/
  ntc = 0;
  tj->Ttarget = -100.0;
  tj->Tinit = -100.0;
  tj->Ptarget = 1.0;
  tj->rattletol = 1.0e-6;
  tj->MaxRattleIter = 100;
  tj->RemoveMomentum = 0;
  tj->ioutfm = 0;
  tj->ntt = 0;
  tj->ntp = 0;
  tj->barostat = 1;
  tj->vrand = 1000;
  tj->BerendsenTCoupl = 0.4;
  tj->BerendsenPTime = 1.0;
  tj->BerendsenPCoupl = 44.6;
  tj->MCBarostatFac[0] = 2.0e-3;
  tj->MCBarostatFac[1] = -1.0;
  tj->MCBarostatFac[2] = -1.0;
  tj->MCBarostatFreq = 100;
  tj->npth.TauT = 1.0;
  tj->npth.TauP = 1.0;
  tj->lnth.gamma_ln = 0.0;

  /*** Default settings for trajectory control ***/
  tj->nstep = 1;
  tj->nfistep = 0;
  tj->currfi = 0;
  tj->ntwr = 0;
  tj->ntwx = 0;
  tj->ntwv = 0;
  tj->ntwf = 0;
  tj->ntpr = 0;
  tj->irest = 0;
  tj->topchk = 1;

  /*** Default settings for energy minimization ***/
  tj->EMinStep = 1.0;
  tj->EMinStep0 = 1.0;
  tj->EMinTol = 1.0e-4;

  /*** Default settings for thermodynamic integration ***/
  tj->TI = 0;
  tj->mxorder = 1;
  tj->nsynch = 1000;
  EXigseed = (tj->igseed == 72177) ? 0 : 1;
  EXlambda = (fabs(tj->lambda) < 1.0e-8) ? 0 : 1;

  /*** Allocate space for topology atom masks ***/
  for (i = 0; i < 2; i++) {
    tp[i].norattlemask = (char*)calloc(MAXNAME, sizeof(char));
    tp[i].rattlemask = (char*)calloc(MAXNAME, sizeof(char));
  }
  collect = AdvanceToSegment(inp, "cntrl", 1);
  while (collect == 1) {

    /*** Read the next line ***/
    fgets(line, MAXLINE, inp);
    RemoveWhiteSpace(line, MAXLINE);
    RemoveComments(line);

    /*** Break if the line is "&end" ***/
    collect = DetectNamelistEnd(line, "GetCntrlNamelist");
    if (collect == 0) {
      continue;
    }

    /*** Eliminate and add spaces between special characters ***/
    /*** "=", "\n", and ","                                  ***/
    NixCommaCarriage(line);
    EqualSpace(line);
    lwords = ParseWords(line);

    /*** Dynamics / Hamiltonian directives ***/
    SeekReal(lwords, &tp->lj14fac, "Vdw14Fac", "scnb");
    SeekReal(lwords, &tp->elec14fac, "Elec14Fac", "scee");
    SeekReal(lwords, &tj->rattletol, "RattleTol", "tol");
    SeekReal(lwords, &dcinp->Ecut, "ElecCut", "es_cutoff");
    SeekReal(lwords, &dcinp->Vcut, "VdwCut", "vdw_cutoff");
    SeekReal(lwords, &EVcut, "DirectCutoff", "cut");
    SeekReal(lwords, &tj->starttime, "StartTime", "t");
    SeekReal(lwords, &tj->dt, "TimeStep", "dt");
    SeekReal(lwords, &dcinp->MaxDens, "MaxDensity", "rho");
    SeekInt(lwords, &tj->irest, "RestartMD", "irest");
    SeekInt(lwords, &tp->rattle, "DoRATTLE", "rigidbond");
    SeekInt(lwords, &tp->settle, "DoSETTLE", "rigidwat");
    SeekInt(lwords, &tj->RemoveMomentum, "ZeroMomentum", "nscm");
    SeekInt(lwords, &tj->mode, "RunMode", "imin");
    SeekInt(lwords, &ntc, "ShakeLevel", "ntc");
    if (EXigseed == 0) {
      SeekInt(lwords, &tj->igseed, "RandomSeed", "ig");
    }

    /*** Thermodynamic integration directives ***/
    SeekInt(lwords, &tj->TI, "RunTI", "icfe");
    SeekInt(lwords, &tj->mxorder, "MixOrder", "klambda");
    SeekInt(lwords, &tj->nsynch, "SynchTI", "nsynch");
    if (EXlambda == 0) {
      SeekReal(lwords, &tj->lambda, "MixFactor", "clambda");
    }

    /*** Thermostat and barostat directives ***/
    SeekReal(lwords, &tj->Ttarget, "Temperature", "temp0");
    SeekReal(lwords, &tj->Ptarget, "Pressure", "pres0");
    SeekReal(lwords, &tj->Tinit, "StartTemp", "tempi");
    SeekReal(lwords, &tj->BerendsenTCoupl, "BerendsenTC", "taup");
    SeekReal(lwords, &tj->BerendsenPCoupl, "BerendsenPC", "comp");
    SeekReal(lwords, &tj->npth.TauT, "HooverTC", "tauthv");
    SeekReal(lwords, &tj->npth.TauP, "HooverPC", "tauphv");
    SeekReal(lwords, &tj->lnth.gamma_ln, "LangevinFreq", "gamma_ln");
    SeekReal(lwords, &tj->MCBarostatFac[0], "MCBarostatPC", "mccomp");
    SeekReal(lwords, &tj->MCBarostatFac[0], "MCBarostatPCX", "mccompx");
    SeekReal(lwords, &tj->MCBarostatFac[1], "MCBarostatPCY", "mccompy");
    SeekReal(lwords, &tj->MCBarostatFac[2], "MCBarostatPCZ", "mccompz");
    SeekInt(lwords, &tj->ntt, "Thermostat", "ntt");
    SeekInt(lwords, &tj->ntp, "CoordRscl", "ntp");
    SeekInt(lwords, &tj->barostat, "Barostat", "barostat");
    SeekInt(lwords, &tj->vrand, "RandomReset", "vrand");
    SeekInt(lwords, &tj->MCBarostatFreq, "MCBarostatFrq", "mcbfrq");

    /*** Output management directives ***/
    SeekLLInt(lwords, &tj->nstep, "StepCount", "nstlim");
    SeekLLInt(lwords, &tj->nfistep, "FileStepCount", "nfistep");
    SeekInt(lwords, &tj->ntpr, "WriteDiagnostics", "ntpr");
    SeekInt(lwords, &tj->ntwx, "WriteCoordinates", "ntwx");
    SeekInt(lwords, &tj->ntwv, "WriteVelocities", "ntwv");
    SeekInt(lwords, &tj->ntwf, "WriteForces", "ntwf");
    SeekInt(lwords, &tj->ntwr, "WriteRestart", "ntwr");
    SeekInt(lwords, &tj->topchk, "TopologyCheck", "tchk");
    SeekInt(lwords, &tj->ioutfm, "CoordFormat", "ioutfm");
    SeekString(lwords, tp[0].WaterName, "WaterName1", "watnam");
    SeekString(lwords, tp[1].WaterName, "WaterName2", "watnam");
    SeekString(lwords, tp[0].rattlemask, "RattleMask1", "rattlemask");
    SeekString(lwords, tp[1].rattlemask, "RattleMask2", "rattlemask");
    SeekString(lwords, tp[0].norattlemask, "NoRattleMask1", "norattlemask");
    SeekString(lwords, tp[1].norattlemask, "NoRattleMask2", "norattlemask");

    /*** Free allocated memory ***/
    DestroyCmat(&lwords);
  }

  /*** Restate control data in internal units/conventions ***/
  if (EVcut > 0.0) {
    dcinp->Ecut = EVcut;
    dcinp->Vcut = EVcut;
  }
  tp->lj14fac = 1.0 - 1.0/tp->lj14fac;
  tp->elec14fac = 1.0 - 1.0/tp->elec14fac;
  dcinp->Mcut = MAX(dcinp->Ecut, dcinp->Vcut);
  dcinp->invMcut = 1.0/dcinp->Mcut;
  dcinp->invEcut = 1.0/dcinp->Ecut;
  tj->BerendsenPCoupl *= 1.0e-6;
  if (tj->Ttarget < 0.0) {
    tj->Ttarget = 298.0;
    if (tj->Tinit < 0.0) {
      tj->Tinit = 0.0;
    }
  }
  else {
    if (tj->Tinit < 0.0) {
      tj->Tinit = tj->Ttarget;
    }
  }
  if (tj->TI == 1) {
    pval = PascalTriangle(tj->mxorder+1);
    lampow = 1.0;
    lampowm1 = 1.0/tj->lambda;
    mfac = 1.0;
    mxA = 0.0;
    dmxA = 0.0;
    for (i = 0; i <= tj->mxorder; i++) {
      mxA += mfac*lampow*pval[i];
      dmxA += i*mfac*lampowm1*pval[i];
      mfac *= -1.0;
      lampow *= tj->lambda;
      lampowm1 *= tj->lambda;
    }
    tj->mxA = mxA;
    tj->mxB = 1.0 - mxA;
    tj->dmxA = dmxA;
    tj->dmxB = -dmxA;
  }
  if (ntc > 0) {

    /*** ntc settings can determine the SHAKE and RATTLE ***/
    /*** preferences, and ntc will take precedence.      ***/
    if (ntc == 2) {
      if (tp->rattle == 0) {
	tp->rattle = 1;
      }
      if (tp->settle == 0) {
	tp->settle = 1;
      }
    }
    else if (ntc > 2) {
      printf("GetCntrlNamelist >> Error.  ntc > 2 is not supported.\n");
    }
  }

  /*** Check input: gamma_ln must for Langevin thermostat ***/
  if (tj->ntt == 3) {
    if (tj->lnth.gamma_ln < 1.0e-8) {
      printf("GetCntrlNamelist >> Error.  Langevin thermostat cannot be "
	     "called with a\nGetCntrlNamelist >> collision frequency of "
	     "%8.4lf / ps.\n", tj->lnth.gamma_ln);
      exit(1);
    }
  }

  /*** Check input: van-der Waals cutoff ***/
  /*** must be .GE. electrostatic cutoff ***/
  if (dcinp->Vcut < dcinp->Ecut) {
    printf("GetCntrlNamelist >> Error.  Vdw cutoff (%9.4lf) must be >= "
	   "electrostatic\nGetCntrlNamelist >> cutoff (%9.4lf)\n.",
	   dcinp->Vcut, dcinp->Ecut);
    exit(1);
  }

  /*** Check input: Monte-Carlo barostat is   ***/
  /*** the only permissible option with ntp=2 ***/
  if (tj->ntp == 2 && tj->barostat != 2) {
    printf("GetCntrlNamelist >> Error.  The Monte-Carlo barostat must be "
	   "used for\nGetCntrlNamelist >> anisotropic rescaling.\n");
    exit(1);
  }

  /*** Check input: certain RATTLE procedures are not allowed ***/
  if (tp->rattle == 2) {
    printf("GetCntrlNamelist >> Error.  Constraining all bonds in not "
	   "supported.\nGetCntrlNamelist >> Set DoRATTLE/ntc = 0 or 1.\n");
    exit(1);
  }

  /*** Check input: problems with file segmentation ***/
  if (tj->nfistep > 0) {
    if (tj->nstep % tj->nfistep != 0) {
      printf("GetCntrlNamelist >> Warning.  nstlim (mod) nfistep is nonzero.\n"
	     "GetCntrlNamelist >> nfistep / FileStepCount set to zero.\n");
      tj->nfistep = 0;
    }
    if (tj->ntwr > tj->nfistep || tj->ntwx > tj->nfistep ||
	tj->ntwv > tj->nfistep || tj->ntwf > tj->nfistep ||
	tj->ntpr > tj->nfistep) {
      printf("GetCntrlNamelist >> Warning.  When the run output is split into "
	     "segments,\nGetCntrlNamelist >> ntwr, ntpr, ntwx, ntwv, and "
	     "ntwf cannot exceed nfistep.\n");
      if (tj->ntwr > tj->nfistep) {
	tj->ntwr = tj->nfistep;
      }
      if (tj->ntwx > tj->nfistep) {
	tj->ntwx = tj->nfistep;
      }
      if (tj->ntwv > tj->nfistep) {
	tj->ntwv = tj->nfistep;
      }
      if (tj->ntwf > tj->nfistep) {
	tj->ntwf = tj->nfistep;
      }
      if (tj->ntpr > tj->nfistep) {
	tj->ntpr = tj->nfistep;
      }
    }
  }
  tj->ntpr = SetOutputFrequency(tj, tj->ntpr, "ntpr");
  tj->ntwx = SetOutputFrequency(tj, tj->ntwx, "ntwx");
  tj->ntwv = SetOutputFrequency(tj, tj->ntwv, "ntwv");
  tj->ntwf = SetOutputFrequency(tj, tj->ntwf, "ntwf");
  tj->ntwr = SetOutputFrequency(tj, tj->ntwr, "ntwr");

  /*** After reading the control namelist, we should know ***/
  /*** whether TI is active or not, and therefore whether ***/
  /*** we are dealing with two systems or just one.  We   ***/
  /*** also know exactly how many systems have been       ***/
  /*** specified, and how many topologies are available.  ***/
  tj->nsys = tj->ipcname.row;
  if (tj->TI == 1 && tj->nsys != 2) {
    printf("GetCntrlNamelist >> Error.  Thermodynamic integration has been "
	   "activated.\nGetCntrlNamelist >> Two sets of initial coordinates "
	   "are required.\n");
    exit(1);
  }
  if (tj->TI == 1 && tp[1].source[0] == '\0') {
    printf("GetCntrlNamelist >> Error.  Thermodynamic integration has been "
	   "activated.\nGetCntrlNamelist >> Two topologies are required.\n");
    exit(1);
  }
  tj->ntop = (tj->TI == 1) ? 2 : 1;
  if (tj->ntop == 1) {
    free(tp[1].rattlemask);
    free(tp[1].norattlemask);
  }

  /*** Check the output stream counts ***/
  if (tj->ntwr > 0 && CountOutputStreams(&tj->rstbase) != tj->nsys) {
    printf("GetCntrlNamelist >> Error. %d restart file base names "
	   "were specified\nGetCntrlNamelist >> for %d systems.\n",
	   CountOutputStreams(&tj->rstbase), tj->nsys);
    exit(1);
  }
  if (tj->ntwx > 0 && CountOutputStreams(&tj->trjbase) != tj->nsys) {
    printf("GetCntrlNamelist >> Error. %d trajectory file base names "
	   "were specified\nGetCntrlNamelist >> for %d systems.\n",
	   CountOutputStreams(&tj->trjbase), tj->nsys);
    exit(1);
  }
  if (tj->ntwv > 0 && CountOutputStreams(&tj->velbase) != tj->nsys) {
    printf("GetCntrlNamelist >> Error. %d velocity file base names "
	   "were specified\nGetCntrlNamelist >> for %d systems.\n",
	   CountOutputStreams(&tj->velbase), tj->nsys);
    exit(1);
  }
  if (tj->ntwf > 0 && CountOutputStreams(&tj->frcbase) != tj->nsys) {
    printf("GetCntrlNamelist >> Error. %d force file base names "
	   "were specified\nGetCntrlNamelist >> for %d systems.\n",
	   CountOutputStreams(&tj->frcbase), tj->nsys);
    exit(1);
  }

  /*** Check for output file suffixes ***/
  ExpandSuffixes(&tj->rstsuff, &tj->rstbase);
  ExpandSuffixes(&tj->trjsuff, &tj->trjbase);
  ExpandSuffixes(&tj->velsuff, &tj->velbase);
  ExpandSuffixes(&tj->frcsuff, &tj->frcbase);
}

/***=======================================================================***/
/*** GetTopolNamelist: this function reads in additional information about ***/
/***                   a topology that cannot be specified in the typical  ***/
/***                   AMBER topology format.  Essentially, information    ***/
/***                   provided in the &topol namelist alerts MDGX that    ***/
/***                   it's dealing with a non-standard topology, and to   ***/
/***                   look for certain things when it actually does read  ***/
/***                   in the topology file.                               ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tp:      the topology structure                                     ***/
/***=======================================================================***/
static void GetTopolNamelist(prmtop *tp, FILE *inp)
{
  int collect;
  char line[MAXLINE];
  cmat lwords;

  /*** Default settings ***/
  tp->ljbuck = 0;
  tp->qshape = 0;
  collect = AdvanceToSegment(inp, "cntrl", 1);
  while (collect == 1) {

    /*** Read the next line ***/
    fgets(line, MAXLINE, inp);
    RemoveWhiteSpace(line, MAXLINE);
    RemoveComments(line);

    /*** Break if the line is "&end" ***/
    collect = DetectNamelistEnd(line, "GetTopolNamelist");
    if (collect == 0) {
      continue;
    }

    /*** Eliminate and add spaces between special characters ***/
    /*** "=", "\n", and ","                                  ***/
    NixCommaCarriage(line);
    EqualSpace(line);
    lwords = ParseWords(line);

    /*** Information about the topology that cannot go in a ***/
    /*** standard AMBER prmtop file                         ***/
    SeekInt(lwords, &tp->ljbuck, "vdWstyle", "ljstyle");
    SeekInt(lwords, &tp->qshape, "ChargeStyle", "qshape");

    /*** Free allocated memory ***/
    DestroyCmat(&lwords);
  }
}

/***=======================================================================***/
/*** GetEwaldNamelist: this function reads in Ewald parameters, as would   ***/
/***                   be found in the &ewald namelist in PMEMD or SANDER  ***/
/***                   command input files.                                ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   dcinp:   direct space command information                           ***/
/***   rcinp:   reciprocal space command information                       ***/
/***   inp:     the input file                                             ***/
/***=======================================================================***/
static void GetEwaldNamelist(dircon *dcinp, reccon *rcinp, FILE *inp)
{
  int i, collect, genordr;
  char line[MAXLINE];
  cmat lwords;

  /*** Default settings for SPME ***/
  rcinp->ng = (int*)malloc(3*sizeof(int));
  rcinp->ng[0] = -1;
  rcinp->ng[1] = -1;
  rcinp->ng[2] = -1;
  rcinp->ordr[0] = 4;
  rcinp->ordr[1] = 4;
  rcinp->ordr[2] = 4;
  dcinp->Dtol = 6.0e-6;
  dcinp->lkpspc = 0.0625;
  dcinp->ewcoeff = -1.0;
  dcinp->LRvdw = 1;
  rcinp->S = -1.0;
  genordr = -1;

  /*** Default settings for MLE ***/
  rcinp->ggordr = 8;
  rcinp->nlev = 1;
  rcinp->nslab = 1;
  rcinp->nstrip = 1;
  for (i = 0; i < 4; i++) {
    rcinp->cfac[i] = -1.0;
    rcinp->PadYZ[i] = -1;
  }
  rcinp->cfac[0] = 1.0;
  rcinp->PadYZ[0] = 1;

  /*** Scan the input file ***/
  collect = AdvanceToSegment(inp, "ewald", 1);
  while (collect == 1) {

    /*** Read the next line ***/
    fgets(line, MAXLINE, inp);
    RemoveWhiteSpace(line, MAXLINE);
    RemoveComments(line);

    /*** Break if the line is "&end" ***/
    collect = DetectNamelistEnd(line, "GetEwaldNamelist");
    if (collect == 0) {
      continue;
    }

    /*** Eliminate and add spaces between special characters ***/
    /*** "=", "\n", and ","                                  ***/
    NixCommaCarriage(line);
    EqualSpace(line);
    lwords = ParseWords(line);

    /*** PME directives ***/
    SeekReal(lwords, &dcinp->Dtol, "DSumTol", "dsum_tol");
    SeekReal(lwords, &dcinp->ewcoeff, "EwaldCof", "ew_coeff");
    SeekReal(lwords, &rcinp->S, "Sigma", "sigma");
    SeekReal(lwords, &dcinp->lkpspc, "SplnSpc", "eedtbdns");
    SeekInt(lwords, &rcinp->ng[0], "MeshDimX", "nfft1");
    SeekInt(lwords, &rcinp->ng[1], "MeshDimY", "nfft2");
    SeekInt(lwords, &rcinp->ng[2], "MeshDimZ", "nfft3");
    SeekInt(lwords, &rcinp->ordr[0], "OrderX", "ordr1");
    SeekInt(lwords, &rcinp->ordr[1], "OrderY", "ordr2");
    SeekInt(lwords, &rcinp->ordr[2], "OrderZ", "ordr3");
    SeekInt(lwords, &genordr, "Order", "order");

    /*** Long-ranged vdW directives ***/
    SeekInt(lwords, &dcinp->LRvdw, "vdwmeth", "vdwMethod");

    /*** MLE directives ***/
    SeekInt(lwords, &rcinp->nlev, "EwaldLevels", "nlev");
    SeekInt(lwords, &rcinp->PadYZ[0], "Padding1", "lpad1");
    SeekInt(lwords, &rcinp->PadYZ[1], "Padding2", "lpad2");
    SeekInt(lwords, &rcinp->PadYZ[2], "Padding3", "lpad3");
    SeekReal(lwords, &rcinp->cfac[1], "Spread2", "cfac2");
    SeekReal(lwords, &rcinp->cfac[2], "Spread3", "cfac3");
    SeekReal(lwords, &rcinp->cfac[3], "Spread4", "cfac4");
    SeekInt(lwords, &rcinp->ggordr, "GridOrder", "ggordr");

    /*** Free allocated memory ***/
    DestroyCmat(&lwords);
  }

  /*** If a general order is specified, apply it ***/
  if (genordr > 0) {
    rcinp->ordr[0] = genordr;
    rcinp->ordr[1] = genordr;
    rcinp->ordr[2] = genordr;
  }

  /*** Compute the Ewald coefficient and then the Gaussian spread ***/
  if (rcinp->S <= 0.0 && dcinp->ewcoeff <= 0.0) {
    dcinp->ewcoeff = EwaldCoefficient(dcinp->Ecut, dcinp->Dtol);
    rcinp->S = 0.5/dcinp->ewcoeff;
  }
  else {
    if (rcinp->S > 0.0) {
      dcinp->ewcoeff = 0.5/rcinp->S;
    }
    else {
      rcinp->S = 0.5/dcinp->ewcoeff;
    }
    dcinp->Dtol = (1.0 - erf(dcinp->ewcoeff*dcinp->Ecut))/dcinp->Ecut;
  }

  /*** Check input ***/
  if (dcinp->Dtol < 0.0) {
    printf("GetEwaldNamelist >> Error.  Direct sum tolerance must be a "
	   "positive real value.GetEwaldNamelist >> Value of %9.6lf is "
	   "unacceptable.\n", dcinp->Dtol);
    exit(1);
  }
  for (i = 0; i < 3; i++) {
    if (rcinp->ordr[i] < 3) {
      printf("GetEwaldNamelist >> Error.  Order must be at least 3.\n"
	     "GetEwaldNamelist >> Order was specified as %d.\n",
	     rcinp->ordr[i]);
      exit(1);
    }
  }
}

/***=======================================================================***/
/*** GetForceNamelist: this function reads the &force namelist, which has  ***/
/***                   no analog other than perhaps &debug in SANDER.      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tj:      the trajectory control information (the force report file  ***/
/***              is considered a trajectory type of output)               ***/
/***   inp:     the input file                                             ***/
/***=======================================================================***/
static void GetForceNamelist(trajcon *tj, FILE *inp)
{
  int collect;
  char line[MAXLINE];
  cmat lwords;

  /*** By default, all output is written to the variable "S" ***/
  /*** in the force report                                   ***/
  sprintf(tj->DMPvar, "S");
  tj->DMPcrd = 1;
  tj->DMPbond = 1;
  tj->DMPangl = 1;
  tj->DMPdihe = 1;
  tj->DMPrelec = 1;
  tj->DMPdelec = 1;
  tj->DMPvdw = 1;
  tj->DMPall = 1;

  collect = AdvanceToSegment(inp, "force", 1);
  while (collect == 1) {

    /*** Read the next line ***/
    fgets(line, MAXLINE, inp);
    RemoveWhiteSpace(line, MAXLINE);
    RemoveComments(line);

    /*** Break if the line is "&end" ***/
    collect = DetectNamelistEnd(line, "GetForceNamelist");
    if (collect == 0) {
      continue;
    }

    /*** Eliminate and add spaces between special characters ***/
    /*** "=", "\n", and ","                                  ***/
    NixCommaCarriage(line);
    EqualSpace(line);
    lwords = ParseWords(line);

    /*** Directives for output organization ***/
    SeekString(lwords, tj->DMPvar, "VarName", "var");

    /*** Directives for output suppression ***/
    SeekInt(lwords, &tj->DMPcrd, "DumpCoord", "dumpcrd");
    SeekInt(lwords, &tj->DMPbond, "DumpBond", "dumpbond");
    SeekInt(lwords, &tj->DMPangl, "DumpAngl", "dumpangl");
    SeekInt(lwords, &tj->DMPdihe, "DumpDihe", "dumpdihe");
    SeekInt(lwords, &tj->DMPrelec, "DumpRElec", "dumprelec");
    SeekInt(lwords, &tj->DMPdelec, "DumpDElec", "dumpdelec");
    SeekInt(lwords, &tj->DMPvdw, "DumpVdw", "dumpvdw");
    SeekInt(lwords, &tj->DMPall, "DumpAll", "dumpall");

    /*** Free allocated memory ***/
    DestroyCmat(&lwords);
  }
}

/***=======================================================================***/
/*** GetParamNamelist: this function reads a &param namelist, which has no ***/
/***                   analog in sander.  It is distinct from the &fit     ***/
/***                   namelist reader, because although they both do      ***/
/***                   fitting the methods are distinct.                   ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tp:      the topology struct (stores its own source file)           ***/
/***   myfit:   the fitting struct                                         ***/
/***=======================================================================***/
static void GetParamNamelist(prmset *myparms, trajcon *tj, FILE *inp)
{
  int collect, currconf, nconf, maxconf;
  int maxbadj, maxaadj, maxhadj;
  char line[MAXLINE], etrgbuff[64];
  cmat lwords;

  /*** Default values ***/
  myparms->ep[0] = '\0';
  myparms->sr[0] = '\0';
  myparms->ao[0] = '\0';
  sprintf(myparms->WaterName, "WAT");
  sprintf(myparms->NrgUnits, "kcal/mol");
  myparms->ljbuck = 0;
  myparms->lj14fac = 2.0;
  myparms->elec14fac = 1.2;
  myparms->nbadj = 0;
  myparms->naadj = 0;
  myparms->nhadj = 0;
  myparms->cnstB = -1.0;
  myparms->cnstA = -1.0;
  myparms->cnstH = -1.0;
  myparms->FitAllBonds = 0;
  myparms->FitAllAngles = 0;
  myparms->FitAllTorsions = 0;
  myparms->mmtol = 30.0;
  myparms->esigtol = 5.0;
  myparms->fitscnb = 0;
  myparms->fitscee = 0;
  myparms->reportall = 1;
  myparms->verbose = 1;
  myparms->RemoveOutliers = 0;
  myparms->ititl = (char*)malloc(MAXLINE*sizeof(double));
  sprintf(myparms->ititl, "Generated by mdgx executing %s.\n", tj->inpname);
  myparms->icomm = (char*)malloc(MAXLINE*sizeof(double));
  sprintf(myparms->icomm, "FITTED BY MDGX");

  /*** If no parameter fitting input was received, de-allocate ***/
  /*** character matrices meant to hold such inputs.           ***/
  collect = AdvanceToSegment(inp, "param", 1);
  if (collect == 0) {
    return;
  }

  /*** Pre-allocate arrays for adjustable terms specifications ***/
  myparms->badj = (xbonddef*)malloc(32*sizeof(xbonddef));
  myparms->aadj = (xangldef*)malloc(32*sizeof(xangldef));
  myparms->hadj = (torterm*)malloc(32*sizeof(torterm));
  maxbadj = 32;
  maxaadj = 32;
  maxhadj = 32;

  /*** The presence of a &param namelist will ***/
  /*** override directives for dynamics.      ***/
  tj->mode = 4;
  nconf = 0;
  maxconf = 32;
  myparms->conf = (mmsys*)malloc(maxconf*sizeof(mmsys));
  while (collect == 1) {

    /*** Read the next line ***/
    fgets(line, MAXLINE, inp);
    RemoveWhiteSpace(line, MAXLINE);
    RemoveComments(line);

    /*** Break if the line is "&end" ***/
    collect = DetectNamelistEnd(line, "GetParamNamelist");
    if (collect == 0) {
      continue;
    }

    /*** Eliminate and add spaces between      ***/
    /*** special characters "=", "\n", and "," ***/
    NixCommaCarriage(line);
    EqualSpace(line);
    lwords = ParseWords(line);

    /*** Input directives ***/
    currconf = nconf;
    SeekStringTripletInc(lwords, myparms->conf[nconf].tp.source,
			 myparms->conf[nconf].crdsrc, etrgbuff, "System",
			 "sys", &nconf);
    SeekInt(lwords, &myparms->ljbuck, "vdWstyle", "ljstyle");
    SeekInt(lwords, &myparms->FitAllBonds, "FitBonds", "bonds");
    SeekInt(lwords, &myparms->FitAllAngles, "FitAngles", "angles");
    SeekInt(lwords, &myparms->FitAllTorsions, "FitTorsions", "torsions");
    SeekInt(lwords, &myparms->fitscnb, "FitLJ14", "fitscnb");
    SeekInt(lwords, &myparms->fitscee, "FitEE14", "fitscee");
    SeekInt(lwords, &myparms->reportall, "ReportAll", "repall");
    SeekInt(lwords, &myparms->verbose, "ShowProgress", "verbose");
    SeekInt(lwords, &myparms->RemoveOutliers, "ElimOutliers", "elimsig");
    SeekString(lwords, myparms->WaterName, "WaterName", "watnam");
    SeekString(lwords, myparms->NrgUnits, "EnergyUnits", "eunits");
    SeekString(lwords, myparms->sr, "SeriesReport", "srep");
    SeekString(lwords, myparms->ao, "AccReport", "accrep");
    SeekString(lwords, myparms->ititl, "ParmTitle", "title");
    SeekReal(lwords, &myparms->lj14fac, "Vdw14Fac", "scnb");
    SeekReal(lwords, &myparms->elec14fac, "Elec14Fac", "scee");
    SeekReal(lwords, &myparms->cnstB, "BondRest", "brst");
    SeekReal(lwords, &myparms->cnstA, "AngleRest", "arst");
    SeekReal(lwords, &myparms->cnstH, "TorsionRest", "hrst");
    SeekReal(lwords, &myparms->mmtol, "ConfTol", "ctol");
    SeekReal(lwords, &myparms->esigtol, "EOutlier", "esigtol");
    SeekTorsionID(lwords, myparms, "FitH", "fith", &maxhadj); 
    if (nconf > currconf) {

      /*** There was a new conformation read; commit ***/
      /*** the target energy to memory as a double   ***/
      myparms->conf[currconf].etrg = atof(etrgbuff);
    }
    if (nconf == maxconf) {
      maxconf += 32;
      myparms->conf = (mmsys*)realloc(myparms->conf, maxconf*sizeof(mmsys));
    }

    /*** Free allocated memory ***/
    DestroyCmat(&lwords);
  }

  /*** Remember the number of conformations ***/
  myparms->nconf = nconf;
}

/***=======================================================================***/
/*** GetFitNamelist: this function reads a &fit namelist, which has no     ***/
/***                 analog whatsoever in sander.                          ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tp:      the topology struct (stores its own source file)           ***/
/***   tj:      trajectory control information                             ***/
/***   myfit:   the fitting struct (stores lots of input data files)       ***/
/***=======================================================================***/
static void GetFitNamelist(fset *myfit, trajcon *tj, FILE *inp)
{
  int i, collect;
  int* exg;
  double* qsumtrg;
  char line[MAXLINE], maxmemstr[MAXNAME];
  cmat qeqstrings, qminstrings, qsumstrings, lwords;

  /*** If no charge fitting input was received, de-allocate ***/
  /*** character matrices meant to hold such inputs.        ***/
  collect = AdvanceToSegment(inp, "fit", 1);
  if (collect == 0) {
    return;
  }

  /*** Default input values ***/
  myfit->Rc = 3.0;
  myfit->Rmax = 6.0;
  myfit->psig = 3.16435;
  myfit->peps = 0.16275;
  myfit->prbarm = 0.9572;
  myfit->stericlim = 3.0;
  myfit->nfitpt = 1000;
  myfit->flimit = 0.4;
  myfit->totalq = (double*)calloc(MAXSYS, sizeof(double));
  myfit->wt = (double*)calloc(MAXSYS, sizeof(double));
  SetDVec(myfit->wt, MAXSYS, 1.0);
  myfit->qminwt = 1.0e-3;
  myfit->fhistbin = 0.1;
  myfit->DispAllDP = 0;
  myfit->epext[0] = '\0';
  myfit->confext[0] = '\0';
  myfit->histfile[0] = '\0';
  myfit->MaxMem = 1073741824;
  myfit->verbose = 1;
  maxmemstr[0] = '\0';

  /*** Allocate restraint and grid name matrices ***/
  myfit->gname = CreateCmat(1, MAXNAME);
  myfit->auxgname = CreateCmat(1, MAXNAME);
  myfit->tpname = CreateCmat(1, MAXNAME);
  myfit->eprule = CreateCmat(1, MAXNAME);
  qeqstrings = CreateCmat(1, MAXNAME);
  qminstrings = CreateCmat(1, MAXNAME);
  qsumstrings = CreateCmat(1, MAXNAME);
  qsumtrg = (double*)calloc(MAXSYS, sizeof(double));

  /*** The presence of a &fit namelist will ***/
  /*** override directives for dynamics.    ***/
  tj->mode = 3;
  while (collect == 1) {

    /*** Read the next line ***/
    fgets(line, MAXLINE, inp);
    RemoveWhiteSpace(line, MAXLINE);
    RemoveComments(line);

    /*** Break if the line is "&end" ***/
    collect = DetectNamelistEnd(line, "GetFitNamelist");
    if (collect == 0) {
      continue;
    }

    /*** Eliminate and add spaces between special characters ***/
    /*** "=", "\n", and ","                                  ***/
    NixCommaCarriage(line);
    EqualSpace(line);
    lwords = ParseWords(line);

    /*** General fitting parameters ***/
    SeekReal(lwords, &myfit->fitprob, "FitProb", "fprob");
    SeekReal(lwords, &myfit->testprob, "TestProb", "tprob");
    SeekReal(lwords, &myfit->flimit, "Proximity", "flim");
    SeekReal(lwords, &myfit->psig, "ProbeSig", "psig");
    SeekReal(lwords, &myfit->peps, "ProbeEps", "peps");
    SeekReal(lwords, &myfit->prbarm, "ProbeArm", "parm");
    SeekReal(lwords, &myfit->stericlim, "StericLimit", "pnrg");
    SeekReal(lwords, &myfit->qminwt, "MinQWeight", "minqwt");
    SeekReal(lwords, &myfit->fhistbin, "HistogramBin", "hbin");
    SeekReal(lwords, &myfit->totalq[0], "TotalQ", "qtot");
    SeekReal(lwords, &myfit->wt[0], "PhiWeight", "phiwt");
    SeekReal(lwords, &qsumtrg[0], "GroupTotalQ", "groupqtot");
    SeekReal(lwords, &myfit->Rc, "AcceptAll", "racc");
    SeekReal(lwords, &myfit->Rmax, "AcceptMax", "rmacc");
    SeekInt(lwords, &myfit->exclusive, "Exclusive", "excl");
    SeekInt(lwords, &myfit->nfitpt, "FitPoints", "nfpt");
    SeekInt(lwords, &myfit->DispAllDP, "ShowDipoles", "dpall");
    SeekInt(lwords, &myfit->verbose, "Verbose", "verbose");
    SeekNReal(lwords, myfit->wt, "PhiWeight", "phiwt", MAXSYS);
    SeekNReal(lwords, qsumtrg, "GroupTotalQ", "groupqtot", MAXSYS);
    SeekString(lwords, maxmemstr, "MaxMemory", "maxmem");

    /*** Output file parameters ***/
    SeekString(lwords, myfit->epext, "EPExtension", "epext");
    SeekString(lwords, myfit->confext, "CFExtension", "cfext");
    SeekString(lwords, myfit->histfile, "HistFile", "hist");

    /*** Lists of input data files and their respective weights ***/
    SeekString(lwords, myfit->gname.map[0], "QMPhi", "phi");
    exg = (int*)calloc(myfit->gname.row, sizeof(int));
    myfit->gname = SeekNString(lwords, &myfit->gname, exg, "QMPhi", "phi");
    free(exg);
    SeekString(lwords, myfit->auxgname.map[0], "AuxPhi", "auxphi");
    exg = (int*)calloc(myfit->auxgname.row, sizeof(int));
    myfit->auxgname = SeekNString(lwords, &myfit->auxgname, exg, "AuxPhi",
				  "auxphi");
    free(exg);

    /*** Lists of restraints ***/
    SeekString(lwords, qeqstrings.map[0], "EqualizeQ", "equalq");
    exg = (int*)calloc(qeqstrings.row, sizeof(int));
    qeqstrings = SeekNString(lwords, &qeqstrings, exg, "EqualizeQ", "equalq");
    free(exg);
    SeekString(lwords, qminstrings.map[0], "MinimizeQ", "minq");
    exg = (int*)calloc(qminstrings.row, sizeof(int));
    qminstrings = SeekNString(lwords, &qminstrings, exg, "MinimizeQ", "minq");
    free(exg);
    SeekString(lwords, qsumstrings.map[0], "GroupSumQ", "sumq");
    exg = (int*)calloc(qsumstrings.row, sizeof(int));
    qsumstrings = SeekNString(lwords, &qsumstrings, exg, "GroupSumQ", "sumq");
    free(exg);

    /*** Lists of System specifications (topologies) ***/
    SeekString(lwords, myfit->tpname.map[0], "Topology", "-p");
    exg = (int*)calloc(myfit->tpname.row, sizeof(int));
    myfit->tpname = SeekNString(lwords, &myfit->tpname, exg, "Topology", "-p");
    free(exg);
    SeekString(lwords, myfit->eprule.map[0], "GridEPRules", "-gxpt");
    exg = (int*)calloc(myfit->eprule.row, sizeof(int));
    myfit->eprule = SeekNString(lwords, &myfit->eprule, exg, "GridEPRules",
				"-gxpt");
    free(exg);

    /*** Free allocated memory ***/
    DestroyCmat(&lwords);
  }

  /*** Parse maximum memory string ***/
  if (maxmemstr[0] != '\0') {
    myfit->MaxMem = ReadNumericalShorthand(maxmemstr);
  }

  /*** Checks on some arrays of strings ***/
  if (myfit->gname.map[0][0] == '\0') {
    printf("GetFitNamelist >> Error.  First electrostatic potential grid "
	   "name is null.\n");
    exit(1);
  }
  if (myfit->auxgname.map[0][0] == '\0' && myfit->auxgname.row == 1) {
    myfit->auxgname.row = 0;
  }
  if (qeqstrings.map[0][0] == '\0' && qeqstrings.row == 1) {
    qeqstrings.row = 0;
  }
  if (qminstrings.map[0][0] == '\0' && qminstrings.row == 1) {
    qminstrings.row = 0;
  }
  if (qsumstrings.map[0][0] == '\0' && qsumstrings.row == 1) {
    qsumstrings.row = 0;
  }

  /*** Checks on other fitting data ***/
  if (myfit->Rc < 0.0) {
    printf("GetFitNamelist >> Error.  Invalid acceptance cutoff %9.6lf "
	   "specified.\n", myfit->Rc);
    exit(1);
  }
  if (myfit->Rmax < 0.0) {
    printf("GetFitNamelist >> Error.  Invalid acceptance maximum range %9.6lf "
	   ".\n", myfit->Rmax);
    exit(1);
  }

  /*** Allocate memory for the fitting data ***/
  myfit->ngrd = CountOutputStreams(&myfit->gname);
  myfit->fitpthist = (int*)calloc(10.0/myfit->fhistbin+1, sizeof(int));

  /*** If auxiliary grids have been specified for any ***/
  /*** of the primary grids, make sure that auxiliary ***/
  /*** grids are available for all primary grids.     ***/
  if (myfit->auxgname.map[0][0] != '\0' || myfit->auxgname.row > 0) {
    if (CountOutputStreams(&myfit->auxgname) != myfit->gname.row) {
      printf("GetFitNamelist >> Error.  If auxiliary grid files are "
	     "specified, they must\nGetFitNamelist >> be specified for all "
	     "primary grids. %d primary grids\nGetFitNamelist >> and %d "
	     "auxiliary grids were specified.\n", myfit->gname.row,
	     myfit->auxgname.row);
      exit(1);
    }
    for (i = 0; i < myfit->auxgname.row; i++) {
      if (myfit->auxgname.map[i] == '\0') {
	printf("GetFitNamelist >> Error.  Auxiliary grid %d not found.\n", i);
	exit(1);
      }
    }
  }

  /*** Process restraint requests ***/
  myfit->nqeq = CountOutputStreams(&qeqstrings);
  myfit->qeq = (nail*)malloc(myfit->nqeq*sizeof(nail));
  for (i = 0; i < qeqstrings.row; i++) {
    if (qeqstrings.map[i][0] == '\0') {
      continue;
    }
    myfit->qeq[i].maskstr = (char*)malloc(MAXNAME*sizeof(char));
    strcpy(myfit->qeq[i].maskstr, qeqstrings.map[i]);
  }
  myfit->nqmin = CountOutputStreams(&qminstrings);
  myfit->qmin = (nail*)malloc(myfit->nqmin*sizeof(nail));
  for (i = 0; i < qminstrings.row; i++) {
    if (qminstrings.map[i][0] == '\0') {
      continue;
    }
    myfit->qmin[i].maskstr = (char*)malloc(MAXNAME*sizeof(char));
    strcpy(myfit->qmin[i].maskstr, qminstrings.map[i]);
  }
  myfit->nqsum = CountOutputStreams(&qsumstrings);
  myfit->qsum = (nail*)malloc(myfit->nqsum*sizeof(nail));
  for (i = 0; i < qsumstrings.row; i++) {
    if (qsumstrings.map[i][0] == '\0') {
      continue;
    }
    myfit->qsum[i].maskstr = (char*)malloc(MAXNAME*sizeof(char));
    strcpy(myfit->qsum[i].maskstr, qsumstrings.map[i]);
    myfit->qsum[i].target = qsumtrg[i];
  }

  /*** Free allocated memory ***/
  free(qsumtrg);
  DestroyCmat(&qeqstrings);
  DestroyCmat(&qminstrings);
  DestroyCmat(&qsumstrings);
}

/***=======================================================================***/
/*** ReadCommandFile: function for reading an mdgx command file.           ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   dcinp:   direct space command information                           ***/
/***   rcinp:   reciprocal space command information                       ***/
/***   tp:      the topology struct (stores its own source file)           ***/
/***   tj:      the trajectory struct (stores the restart file name, as    ***/
/***              well as many other output file names)                    ***/
/***   source:  the command input file                                     ***/
/***=======================================================================***/
void ReadCommFile(dircon *dcinp, reccon *rcinp, prmtop* tp, trajcon *tj,
		  fset *myfit, prmset *myparms, char* source)
{
  int i, maxlinelen;
  char line[MAXLINE];
  FILE *inp;

  /*** Open the input file ***/
  if ((inp = fopen(source, "r")) == NULL) {
    printf("ReadCommFile >> Error.  Command file %s not specified.\n",
	   source);
    exit(1);
  }

  /*** First, scan for file information ***/
  GetDataFileNames(tp, tj, inp);

  /*** Next, scan for &cntrl namelist information ***/
  GetCntrlNamelist(dcinp, tp, tj, inp);

  /*** If there is a second sytem, it may be necessary to update ***/
  /*** certain topology information in the secondary topology.   ***/
  /*** Note that any such information (copied between topologies ***/
  /*** here) is not available to distinguish the initial and     ***/
  /*** final states in a TI calculation.                         ***/
  if (tj->TI == 1) {
    tp[1].lj14fac = tp[0].lj14fac;
    tp[1].elec14fac = tp[0].elec14fac;
    tp[1].rattle = tp[0].rattle;
    tp[1].settle = tp[0].settle;
    strcpy(tp[1].WaterName, tp[0].WaterName);
  }

  /*** Scan for any additional information about the topology. ***/
  /*** This is done here so that the topology reader can go on ***/
  /*** reading standard AMBER topology files, but can also be  ***/
  /*** queued to look for additional information that is not   ***/
  /*** found in a standard topology file.                      ***/
  for (i = 0; i < tj->ntop; i++) {
    GetTopolNamelist(&tp[i], inp);
  }

  /*** Now, scan for &ewald namelist information ***/
  GetEwaldNamelist(dcinp, rcinp, inp);

  /*** Scan for &force namelist information ***/
  GetForceNamelist(tj, inp);

  /*** Scan for &param namelist information ***/
  GetParamNamelist(myparms, tj, inp);

  /*** Scan for &fit namelist information ***/
  GetFitNamelist(myfit, tj, inp);

  /*** Finally, scan the file verbatim into memory ***/
  rewind(inp);
  i = 0;
  maxlinelen = 0;
  while (fgets(line, MAXLINE, inp) != NULL) {
    i++;
    if (strlen(line) > maxlinelen) {
      maxlinelen = strlen(line);
    }
  }
  rewind(inp);
  tj->inptext = CreateCmat(i, maxlinelen+1);
  rewind(inp);
  i = 0;
  while (fgets(line, MAXLINE, inp) != NULL) {
    strcpy(tj->inptext.map[i], line);
    i++;
  }
  fclose(inp);
}
