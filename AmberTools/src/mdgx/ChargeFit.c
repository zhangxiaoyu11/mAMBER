#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mdgxVector.h"
#include "Matrix.h"
#include "Grid.h"
#include "CellManip.h"
#include "VirtualSites.h"
#include "Macros.h"
#include "CrdManip.h"
#include "Parse.h"
#include "Manual.h"
#include "ChargeFit.h"
#include "Trajectory.h"
#include "Topology.h"
#include "Random.h"

/***=======================================================================***/
/*** AssingGridTopologies: when fitting electrostatic potential grids, it  ***/
/***                       may be that some grids pertain to different     ***/
/***                       models than the original topology.  In that     ***/
/***                       case, the topologies must be read into an       ***/
/***                       array and accessed with each grid in the fit.   ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   myfit:   the fitting control data (contains matrices of names for   ***/
/***            topologies and extra points rules files, number of grids)  ***/
/***   tp:      the default topology                                       ***/
/***   tj:      information that must be passed to GetPrmTop in MPI        ***/
/***            implementation                                             ***/
/***=======================================================================***/
#ifdef MPI
static void AssignGridTopologies(fset *myfit, prmtop *tp, trajcon *tj)
#else
static void AssignGridTopologies(fset *myfit, prmtop *tp)
#endif
{
  int i, j;
  int* epcorrchk;

  /*** Bail out if there are no additional topologies specified ***/
  if ((myfit->tpname.map[0][0] == '\0' && myfit->tpname.row == 1) &&
      (myfit->eprule.map[0][0] == '\0' && myfit->eprule.row == 1)) {
    myfit->tpname = ReallocCmat(&myfit->tpname, myfit->ngrd, MAXNAME);
    myfit->eprule = ReallocCmat(&myfit->eprule, myfit->ngrd, MAXNAME);
    for (i = 0; i < myfit->ngrd; i++) {
      strcpy(myfit->tpname.map[i], tp->source);
      strcpy(myfit->eprule.map[i], tp->eprulesource);
    }
    myfit->tpidx = (int*)calloc(myfit->ngrd, sizeof(int));
  }

  /*** Each grid corresponds to a topology.  If a ***/
  /*** topology is not specified for a particular ***/
  /*** grid, use the default topology from the    ***/
  /*** command line or &files namelist.           ***/
  if (myfit->tpname.row > myfit->ngrd) {
    printf("FillTopologyGridNames >> Error.  Topologies have been specified "
	   "for grids\nFillTopologyGridNames >> that do not exist.  Topology "
	   "specified for grid %d of %d.\n", myfit->tpname.row, myfit->ngrd);
    exit(1);
  }
  if (myfit->eprule.row > myfit->ngrd) {
    printf("FillTopologyGridNames >> Error.  Extra point rules have been "
	   "specified for\nFillTopologyGridNames >> grids that do not exist.  "
	   "Rules specified for grid %d of %d.\n", myfit->eprule.row,
	   myfit->ngrd);
    exit(1);
  }
  if (myfit->ngrd > myfit->tpname.row) {
    myfit->tpname = ReallocCmat(&myfit->tpname, myfit->ngrd, MAXNAME);
  }
  if (myfit->ngrd > myfit->eprule.row) {
    myfit->eprule = ReallocCmat(&myfit->eprule, myfit->ngrd, MAXNAME);
  }
  for (i = 0; i < myfit->ngrd; i++) {
    if (myfit->tpname.map[i][0] == '\0') {
      sprintf(myfit->tpname.map[i], tp->source);
      if (myfit->eprule.map[i][0] == '\0' && tp->eprulesource[0] != '\0') {
	sprintf(myfit->eprule.map[i], tp->eprulesource);
      }
    }
  }

  /*** Count the number of unique topologies. ***/
  myfit->tpidx = (int*)malloc(myfit->ngrd*sizeof(int));
  myfit->tpcount = 0;
  SetIVec(myfit->tpidx, myfit->ngrd, -1);
  for (i = 0; i < myfit->ngrd; i++) {
    if (myfit->tpidx[i] >= 0) {
      continue;
    }
    myfit->tpidx[i] = myfit->tpcount;
    for (j = 0; j < myfit->tpname.row; j++) {
      if (strcmp(myfit->tpname.map[i], myfit->tpname.map[j]) == 0) {
	myfit->tpidx[j] = myfit->tpidx[i];
      }
    }
    myfit->tpcount += 1;
  }

  /*** Scan topology names for correspondence with extra points ***/
  /*** rules files.  Do not allow multiple extra points rules   ***/
  /*** files to be specified with a single topology!            ***/
  epcorrchk = CpyIVec(myfit->tpidx, myfit->ngrd);
  for (i = 0; i < myfit->ngrd; i++) {
    if (epcorrchk[i] < 0) {
      continue;
    }
    for (j = 0; j < myfit->ngrd; j++) {
      if (myfit->tpidx[j] == myfit->tpidx[i]) {
	if (strcmp(myfit->eprule.map[i], myfit->eprule.map[j]) == 0) {
	  epcorrchk[j] = -1;
	}
	else {
	  printf("FillTopologyGridNames >> Error.  Extra points rule files\n"
		 "FillTopologyGridNames >> %s and\nFillTopologyGridNames >> "
		 "%s specified\nFillTopologyGridNames >> for topology %s.\n"
		 "FillTopologyGridNames >> Only one extra points rule file "
		 "may be specified\nFillTopologyGridNames >> for each "
		 "topology.\n", myfit->eprule.map[i], myfit->eprule.map[j],
		 myfit->tpname.map[i]);
	  exit(1);
	}
      }
    }
  }
  free(epcorrchk);

  /*** Read all topologies ***/
  myfit->TPbank = (prmtop*)malloc(myfit->tpcount*sizeof(prmtop));
  j = 0;
  for (i = 0; i < myfit->tpcount; i++) {

    /*** Data that gets added to the topology from input file ***/
    myfit->TPbank[i].lj14fac = tp->lj14fac;
    myfit->TPbank[i].elec14fac = tp->elec14fac;
    myfit->TPbank[i].rattle = tp->rattle;
    myfit->TPbank[i].settle = tp->settle;
    myfit->TPbank[i].ljbuck = tp->ljbuck;
    strcpy(myfit->TPbank[i].WaterName, tp->WaterName);
    myfit->TPbank[i].norattlemask = (char*)calloc(MAXNAME, sizeof(char));
    myfit->TPbank[i].rattlemask = (char*)calloc(MAXNAME, sizeof(char));
    strcpy(myfit->TPbank[i].norattlemask, tp->norattlemask);
    strcpy(myfit->TPbank[i].rattlemask, tp->rattlemask);
    while (myfit->tpidx[j] != i) {
      j++;
    }
    strcpy(myfit->TPbank[i].source, myfit->tpname.map[j]);
    strcpy(myfit->TPbank[i].eprulesource, myfit->eprule.map[i]);
    myfit->TPbank[i].lVDWc = CpyDVec(tp->lVDWc, 3*tp->ntypes);
#ifdef MPI
    GetPrmTop(&myfit->TPbank[i], tj, 1);
#else
    GetPrmTop(&myfit->TPbank[i], 1);
#endif
  }

  /*** If there are multiple topologies in this fit, take the     ***/
  /*** total charge constraints from the topologies themselves.   ***/
  /*** Charges will be refit, but their sum should stay the same. ***/
  if (myfit->tpcount > 1) {
    for (i = 0; i < myfit->tpcount; i++) {
      myfit->totalq[i] = DSum(myfit->TPbank[i].Charges,
			      myfit->TPbank[i].natom);
    }
  }
}

/***=======================================================================***/
/*** ColumnsForAtoms: there are a certain number of topologies involved in ***/
/***                  the fit, implying a number of unique atoms.  This    ***/
/***                  routine will determine which atoms are unique, with  ***/
/***                  consideration to user-specified charge equalization  ***/
/***                  constraints, and return a matrix of integers that    ***/
/***                  assigns each atom to a column of the fitting matrix. ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   myfit:      the fitting control data, contains a bank of topologies ***/
/***=======================================================================***/
static imat ColumnsForAtoms(fset *myfit)
{
  int i, j, k, maxatm, icol, natm, offset, maxcol;
  int* atmmask;
  int* atmresv;
  int* colscr;
  imat AC;
  coord crd;

  /*** Allocate the correspondence matrix ***/
  icol = 0;
  maxatm = 0;
  for (i = 0; i < myfit->tpcount; i++) {
    if (myfit->TPbank[i].natom > maxatm) {
      maxatm = myfit->TPbank[i].natom;
    }
    icol += myfit->TPbank[i].natom;
  }
  AC = CreateImat(myfit->tpcount, maxatm);
  maxcol = icol;

  /*** Determine the charge equalization restraints, which may ***/
  /*** span multiple systems, and count the number of columns  ***/
  atmresv = (int*)calloc(maxcol, sizeof(int));
  for (i = 0; i < myfit->nqeq; i++) {
    natm = 0;
    for (j = 0; j < myfit->tpcount; j++) {
      crd = CreateCoord(myfit->TPbank[j].natom);
      atmmask = ParseAmbMask(myfit->qeq[i].maskstr, &myfit->TPbank[j], &crd);
      natm += ISum(atmmask, myfit->TPbank[j].natom);
      free(atmmask);
      DestroyCoord(&crd);
    }
    myfit->qeq[i].atoms = (int*)malloc(natm*sizeof(int));
    myfit->qeq[i].natm = natm;
    natm = 0;
    offset = 0;
    for (j = 0; j < myfit->tpcount; j++) {
      crd = CreateCoord(myfit->TPbank[j].natom);
      atmmask = ParseAmbMask(myfit->qeq[i].maskstr, &myfit->TPbank[j], &crd);
      for (k = 0; k < myfit->TPbank[j].natom; k++) {
	if (atmmask[k] == 1) {
	  myfit->qeq[i].atoms[natm] = k + offset;
	  if (atmresv[k+offset] == 1) {
	    printf("ColumnsForAtoms >> Error.  Overlapping charge "
		   "equalization constraints\nConlumnsForAtoms >> detected.  "
		   "Terminating on restraint %s.\n", myfit->qeq[i].maskstr);
	    exit(1);
	  }
	  atmresv[k+offset] = 1;
	  natm++;
	}
      }
      offset += myfit->TPbank[j].natom;
      free(atmmask);
      DestroyCoord(&crd);
    }

    /*** The total number of columns is decreased by the ***/
    /*** number of atoms in this restraint minus one     ***/
    if (natm >= 2) {
      icol -= natm-1;
    }
  }
  free(atmresv);

  /*** Fill out the correspondence array ***/
  natm = 0;
  for (i = 0; i < myfit->tpcount; i++) {
    for (j = 0; j < myfit->TPbank[i].natom; j++) {
      AC.map[i][j] = natm;
      natm++;
    }
  }
  colscr = CountUp(maxcol);
  for (i = 0; i < myfit->nqeq; i++) {
    for (j = 1; j < myfit->qeq[i].natm; j++) {
      colscr[myfit->qeq[i].atoms[j]] = myfit->qeq[i].atoms[0];
    }
  }

  /*** Adjust the correspondence array to be sequential ***/
  k = 0;
  atmmask = (int*)malloc(maxcol*sizeof(int));
  SetIVec(atmmask, maxcol, 0);
  for (i = 0; i < maxcol; i++) {
    if (atmmask[i] == 1) {
      continue;
    }
    offset = colscr[i];
    for (j = i; j < maxcol; j++) {
      if (atmmask[j] == 0 && colscr[j] == offset) {
	colscr[j] = k;
	atmmask[j] = 1;
      }
    }
    k++;
  }
  free(atmmask);

  /*** Atom correspondence back into the matrix ***/
  k = 0;
  for (i = 0; i < myfit->tpcount; i++) {
    for (j = 0; j < myfit->TPbank[i].natom; j++) {
      AC.map[i][j] = colscr[k];
      k++;
    }
  }

  /*** Free allocated memory ***/
  free(colscr);

  myfit->q2fit = icol;
  return AC;
}

/***=======================================================================***/
/*** ReadEPotGrid: read an electrostatic potential grid as output by the   ***/
/***               Gaussian cubegen utility.  Only formatted cubegen       ***/
/***               outputs are read.                                       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   fname:     the name of the cubegen output file                      ***/
/***   tp:        the name of the topology used to interpret the molecular ***/
/***              coordiantes at the head of the cubegen output file       ***/
/***=======================================================================***/
static fbook ReadEPotGrid(char* fname, prmtop *tp, coord *crd)
{
  int i, j, k, headerfound;
  float *ftmp;
  double lvec[9], orig[3];
  double* nucchg;
  char line[MAXLINE];
  FILE *inp;
  fbook Ue;

  /*** Open the cubegen file ***/
  if ((inp = fopen(fname, "r")) == NULL) {
    printf("ReadEPotGrid >> Error.  File %s not found.\n", fname);
    exit(1);
  }

  /*** Process the cubegen output ***/
  headerfound = 1;
  while (headerfound == 1) {
    fgets(line, 128, inp);
    j = strlen(line);
    headerfound = 0;
    for (i = 0; i < strlen(line); i++) {
      if ((line[i] < '0' || line[i] > '9') && line[i] != '.' &&
	  line[i] != '-' && line[i] != ' ' && line[i] != '\n') {
	headerfound = 1;
	break;
      }
    }
  }
  sscanf(line, "%d%lf%lf%lf", &crd->natom, &orig[0], &orig[1], &orig[2]);
  if (crd->natom > tp->natom) {
    printf("ReadEPotGrid >> Error.  %d atoms present in %s,"
	   "\nReadEPotGrid >> %d in the associated topology %s.\n", crd->natom,
	   fname, tp->natom, tp->source);
    exit(1);
  }
  fscanf(inp, "%d%lf%lf%lf", &i, &lvec[0], &lvec[1], &lvec[2]);
  fscanf(inp, "%d%lf%lf%lf", &j, &lvec[3], &lvec[4], &lvec[5]);
  fscanf(inp, "%d%lf%lf%lf", &k, &lvec[6], &lvec[7], &lvec[8]);
  Ue = CreateFbook(i, j, k);
  ReflectDVec(Ue.lvec.data, lvec, 9);
  ReflectDVec(Ue.orig, orig, 3);
  nucchg = (double*)malloc(crd->natom*sizeof(double));
  for (i = 0; i < crd->natom; i++) {
    fscanf(inp, "%d%lf%lf%lf%lf", &j, &nucchg[i], &crd->loc[3*i],
	   &crd->loc[3*i+1], &crd->loc[3*i+2]);
  }
  if (tp->EPInserted == 1) {
    ExtendCoordinates(crd, tp);
  }
  if (crd->natom != tp->natom) {
    printf("ReadEPotGrid >> Error.  Atom count on grid %d, %d in topology "
	   "%s.\n", crd->natom, tp->natom, tp->source);
    exit(1);
  }
  for (i = 0; i < Ue.pag; i++) {
    for (j = 0; j < Ue.row; j++) {
      ftmp = Ue.map[i][j];
      for (k = 0; k < Ue.col; k++) {
	fscanf(inp, "%f", &ftmp[k]);
      }
    }
  }
  fclose(inp);

  /*** Free allocated memory ***/
  free(nucchg);

  return Ue;
}

/***=======================================================================***/
/*** FitPlaceXpt: place extra points around a molecule read from a cubegen ***/
/***              input file.  This routine places the molecule in a cell  ***/
/***              as part of a 1x1x1 cell grid for the purpose of feeding  ***/
/***              it into the standard extra point placement routines.     ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   crd:      coordinates of the molecule (extra point sites are        ***/
/***             included but all set to zero)                             ***/
/***   tp:       topology decribing the molecule                           ***/
/***=======================================================================***/
static void FitPlaceXpt(coord *crd, prmtop *tp)
{
  int i;
  double *loctmp, *aloc, *bloc, *cloc, *dloc, *eploc;
  expt *tmr;

  /*** Loop over all extra points and place them ***/
  for (i = 0; i < tp->nxtrapt; i++) {
    tmr = &tp->xtrapts[i];
    loctmp = crd->loc;
    aloc = &loctmp[3*tmr->fr1];
    bloc = &loctmp[3*tmr->fr2];
    if (tmr->frstyle > 1) {
      cloc = &loctmp[3*tmr->fr3];
    }
    if (tmr->frstyle == 6) {
      dloc = &loctmp[3*tmr->fr4];
    }
    eploc = &loctmp[3*tmr->atomid];
    XptLocator(aloc, aloc, bloc, cloc, dloc, eploc, eploc, tmr);
  }
}

/***=======================================================================***/
/*** PrepUPot: prepares an electrostatic potential grid read from a file   ***/
/***           for RESP fitting.  Units of the grid potential, scale, and  ***/
/***           molecular coordinates are fitted.                           ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   UPot:    the electrostatic potential grid                           ***/
/***   crd:     molecular coordinates associated with UPot                 ***/
/***   tp:      the topology                                               ***/
/***   UIdx:    grid of indices describing the accessibility and           ***/
/***            availability of points in UPot                             ***/
/***   myfit:   RESP fitting control data                                  ***/
/***   dcinp:   direct space interaction control data                      ***/
/***   rcinp:   reciprocal space interaction control data                  ***/
/***   isaux:   flag to make the function return after merely scaling grid ***/
/***            dimensions and coordinates                                 ***/
/***=======================================================================***/
static void PrepUPot(fbook *UPot, coord *crd, prmtop *tp, cbook *UIdx,
		     fset *myfit, int isaux)
{
  int h, i, j, k, ib, jb, kb, m, ljt;
  int mini, minj, mink, maxi, maxj, maxk, buffi, buffj, buffk;
  double r2, invr2, invr4, sig, eps, Afac, Bfac;
  double atmloc[3], ptloc[3], dr[3], cdepth[3];
  double *dtmp;
  dmat *Utbl;
  dbook Ulj;

  const int nelem = UPot->pag*UPot->row*UPot->col;
  for (j = 0; j < nelem; j++) {
    UPot->data[j] *= H2KCAL;
  }
  for (j = 0; j < 9; j++) {
    UPot->lvec.data[j] *= B2ANG;
  }
  for (j = 0; j < 3; j++) {
    UPot->orig[j] *= B2ANG;
  }
  for (j = 0; j < 3*crd->natom; j++) {
    crd->loc[j] *= B2ANG;
  }
  FitPlaceXpt(crd, tp);

  /*** Bail out if this is just an auxiliary grid ***/
  if (isaux == 1) {
    return;
  }
  Ulj = CreateDbook(UPot->pag, UPot->row, UPot->col, 0);

  /*** Compute the Lennard-Jones of each grid point ***/
  Utbl = &tp->LJutab;
  for (h = 0; h < tp->natom; h++) {
    for (m = 0; m < 3; m++) {
      atmloc[m] = crd->loc[3*h+m];
    }
    ljt = tp->LJIdx[h];
    if (ljt < 0) {
      sig = 0.0;
      eps = 0.0;
    }
    else {
      sig = pow(-1.0*Utbl->map[ljt][2*ljt]/Utbl->map[ljt][2*ljt+1], 1.0/6.0);
      eps = -0.25 * Utbl->map[ljt][2*ljt+1] / pow(sig, 6.0);
      sig = 0.5*(sig + myfit->psig);
      eps = sqrt(eps*myfit->peps);
      Afac = 4.0*eps*pow(sig, 12.0);
      Bfac = 4.0*eps*pow(sig, 6.0);
    }
    for (i = 0; i < Ulj.pag; i++) {
      for (j = 0; j < Ulj.row; j++) {
	dtmp = Ulj.map[i][j];
	for (k = 0; k < Ulj.col; k++) {
	  ptloc[0] = i;
	  ptloc[1] = j;
	  ptloc[2] = k;
	  RotateCrd(ptloc, 1, UPot->lvec);
	  for (m = 0; m < 3; m++) {
	    ptloc[m] += UPot->orig[m];
	    dr[m] = ptloc[m] - atmloc[m];
	  }
	  r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
	  invr2 = 1.0/(r2);
	  invr4 = invr2*invr2;
	  dtmp[k] += Afac*invr4*invr4*invr4 - Bfac*invr4*invr2;
	}
      }
    }
  }

  /*** Compute the grid buffer region for accessibility ***/
  HessianNorms(&UPot->lvec, cdepth);
  buffi = myfit->prbarm/cdepth[0] + 1;
  buffj = myfit->prbarm/cdepth[1] + 1;
  buffk = myfit->prbarm/cdepth[2] + 1;
  for (i = 0; i < Ulj.pag; i++) {
    for (j = 0; j < Ulj.row; j++) {
      for (k = 0; k < Ulj.col; k++) {
	
	/*** If this point is accessible to the probe, ***/
	/*** declare it so and mark all surrounding    ***/
	/*** points accessible if they are within the  ***/
	/*** probe arm's reach.                        ***/
	if (Ulj.map[i][j][k] < myfit->stericlim) {
	  UIdx->map[i][j][k] = 1;
	  continue;
	}

	/*** If we're still here, this point is not  ***/
	/*** accessible to probe but may nonetheless ***/
	/*** be accessible to the probe arm.         ***/
	mini = MAX(i-buffi, 0);
	minj = MAX(j-buffj, 0);
	mink = MAX(k-buffk, 0);
	maxi = MIN(i+buffi, Ulj.pag-1);
	maxj = MIN(j+buffj, Ulj.row-1);
	maxk = MIN(k+buffk, Ulj.col-1);
	for (ib = mini; ib <= maxi; ib++) {
	  for (jb = minj; jb <= maxj; jb++) {
	    for (kb = mink; kb <= maxk; kb++) {

	      /*** In order for this point to be worth investigating, ***/
	      /*** it has to be possible for a probe to be situated   ***/
	      /*** on this point with an acceptable steric energy and ***/
	      /*** then reach its arm out to the point that we've     ***/
	      /*** thus far declared inaccessible.                    ***/
	      if (Ulj.map[ib][jb][kb] > myfit->stericlim) {
		continue;
	      }
	      ptloc[0] = ib - i;
	      ptloc[1] = jb - j;
	      ptloc[2] = kb - k;
	      RotateCrd(ptloc, 1, UPot->lvec);
	      if (ptloc[0]*ptloc[0] + ptloc[1]*ptloc[1] + ptloc[2]*ptloc[2] <
		  myfit->prbarm) {

		/*** Mark the point accessible and  ***/
		/*** bail out of these nested loops ***/
		UIdx->map[i][j][k] = 1;
		ib = maxi+1;
		jb = maxj+1;
		kb = maxk+1;
	      }
	    }
	  }
	}
      }
    }
  }

  /*** Free allocated memory ***/
  DestroyDbook(&Ulj);
}

/***=======================================================================***/
/*** MaskAllAtoms: restraints may span multiple systems.  This function    ***/
/***               examines the restraint mask and applies it to each      ***/
/***               system individually to determine how many unique atoms  ***/
/***               it involves.                                            ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   myfit:       the fitting control data                               ***/
/***   qmask:       the mask of all charges (length the number of columns  ***/
/***                in the fitting matrix)                                 ***/
/***   maskstr:     the mask string (ambmask format)                       ***/
/***   Atm2Col:     matrix describing the mapping of atoms in each system  ***/
/***                to columns of the fitting constraint matrix            ***/
/***   stack:       flag to record counts of how many atoms share each     ***/
/***                unique fitted charge state; set to zero for boolean    ***/
/***                indications as to whether each fitted charge is        ***/
/***                present in the mask, 1 to record counts                ***/
/***   spec:        do the mask only for a specific topology; set to -1 to ***/
/***                make the mask for all topologies                       ***/
/***=======================================================================***/
static int MaskAllAtoms(fset *myfit, int* qmask, char* maskstr,
			imat *Atm2Col, int stack, int spec)
{
  int i, j, natm;
  int *itmp;
  int* atmmask;
  prmtop *tp;
  coord crd;

  SetIVec(qmask, myfit->q2fit, 0);
  for (i = 0; i < myfit->tpcount; i++) {
    if (spec >= 0 && i != spec) {
      continue;
    }
    tp = &myfit->TPbank[i];
    crd = CreateCoord(tp->natom);
    atmmask = ParseAmbMask(maskstr, tp, &crd);
    itmp = Atm2Col->map[i];
    for (j = 0; j < tp->natom; j++) {
      if (atmmask[j] == 1) {
	qmask[itmp[j]] = (stack == 1) ? qmask[itmp[j]] + 1 : 1;
      }
    }
    free(atmmask);
    DestroyCoord(&crd);
  }
  natm = ISum(qmask, myfit->q2fit);

  return natm;
}

/***=======================================================================***/
/*** MakeCnstMatrix: make a constraint matrix for charge fitting based on  ***/
/***                 a series of ambmask strings and other data specified  ***/
/***                 by the user.  This constraint matrix has n+1 columns  ***/
/***                 for fitting n charges, the final column being the     ***/
/***                 corresponding targets on the right hand side of the   ***/
/***                 linear least squares problem (the "b" vector in       ***/
/***                 Ax=b).                                                ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   myfit:    the fitting control data                                  ***/
/***   Atm2Col:  matrix describing the mapping of atoms in each system to  ***/
/***             columns of the fitting constraint matrix                  ***/
/***=======================================================================***/
static dmat MakeCnstMatrix(fset *myfit, imat *Atm2Col, dmat *fitmat)
{
  int i, j, k, nrst, rstcon, natm, tpkey, maxatom;
  int* allqmask;
  int* tmpqmask;
  int* nsamp;
  double *dtmp;
  dmat cnstmat;

  /*** The obligatory total charge restraints, one per system ***/
  cnstmat = CreateDmat(myfit->tpcount, myfit->q2fit+1, 0);
  maxatom = 0;
  for (i = 0; i < myfit->tpcount; i++) {
    maxatom += myfit->TPbank[i].natom;
    for (j = 0; j < myfit->TPbank[i].natom; j++) {
      cnstmat.map[i][Atm2Col->map[i][j]] += 1.0e7;
    }
    cnstmat.map[i][myfit->q2fit] = (1.0e7)*myfit->totalq[i];
  }
  nrst = myfit->tpcount;
  rstcon = myfit->tpcount;

  /*** Count the number of times each atom is sampled in the ***/
  /*** fitting matrix, the number of equations involved, in  ***/
  /*** order to properly scale the restraint weights.        ***/
  nsamp = (int*)calloc(myfit->q2fit, sizeof(int));
  for (i = 0; i < fitmat->row; i++) {
    dtmp = fitmat->map[i];
    for (j = 0; j < myfit->q2fit; j++) {
      if (fabs(dtmp[j]) > 1.0e-8) {
	nsamp[j] += 1;
      }
    }
  }

  /*** Charge minimization restraints ***/
  allqmask = (int*)malloc(maxatom*sizeof(int));
  for (i = 0; i < myfit->nqmin; i++) {

    /*** This restraint may span multiple systems.  We must    ***/
    /*** determine how many unique atoms are being restrained. ***/ 
    natm = MaskAllAtoms(myfit, allqmask, myfit->qmin[i].maskstr,
			Atm2Col, 0, -1);
    if (natm == 0) {
      continue;
    }

    /*** Add new restraints for each specified charge ***/
    nrst += natm;
    cnstmat = ReallocDmat(&cnstmat, nrst, myfit->q2fit+1);
    for (j = 0; j < myfit->q2fit; j++) {
      if (allqmask[j] == 1) {
	dtmp = cnstmat.map[rstcon];
	dtmp[j] = myfit->qminwt*nsamp[j];
	rstcon++;
      }
    }
  }

  /*** Charge group sum restraints ***/
  tmpqmask = (int*)malloc(maxatom*sizeof(int));
  for (i = 0; i < myfit->nqsum; i++) {

    /*** Check to see that this charge group ***/
    /*** does exist in some system.          ***/                  
    natm = 0;
    for (j = 0; j < myfit->tpcount; j++) {
      natm = MaskAllAtoms(myfit, allqmask, myfit->qsum[i].maskstr,
			  Atm2Col, 1, j);
      if (natm > 0) {
	tpkey = j;
	break;
      }
    }
    if (natm == 0) {
      continue;
    }

    /*** It is now known that this charge group exists in ***/
    /*** one or more of the systems, in whole or in part. ***/
    /*** Now we must determine whether the charge group   ***/
    /*** exists entirely or not at all in every system.   ***/
    /*** If this condition is not met, it would result in ***/
    /*** very strange and contradictory restraints, so    ***/
    /*** stop the program if a bad case is detected.      ***/
    for (j = 0; j < myfit->tpcount; j++) {
      natm = MaskAllAtoms(myfit, tmpqmask, myfit->qsum[i].maskstr,
			  Atm2Col, 1, j);
      if (ISum(tmpqmask, natm) == 0) {
	continue;
      }
      for (k = 0; k < myfit->q2fit; k++) {
	if (tmpqmask[k] != allqmask[k]) {
	  printf("MakeCnstMatrix >> Error.  Mismatch in masks generated for "
		 "topologies\nMakeCnstMatrix >> %s and %s.\n",
		 myfit->TPbank[tpkey].source, myfit->TPbank[j].source);
	  exit(1);
	}
      }
    }

    /*** Add the new restraint for this group ***/
    nrst += 1;
    cnstmat = ReallocDmat(&cnstmat, nrst, myfit->q2fit+1);
    dtmp = cnstmat.map[rstcon];
    for (j = 0; j < myfit->q2fit; j++) {
      if (allqmask[j] > 0) {
	dtmp[j] = allqmask[j] * 1.0e7;
      }
    }
    dtmp[myfit->q2fit] = 1.0e7 * myfit->qsum[i].target;
    rstcon++;
  }

  /*** Free allocated memory ***/
  free(allqmask);
  free(tmpqmask);
  free(nsamp);

  return cnstmat;
}

/***=======================================================================***/
/*** SortPtDistance: function called by quicksort for comparing minimum    ***/
/***                 distances between a grid point and the molecule.      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   pt[A,B]:    the atomc structs                                       ***/
/***=======================================================================***/
static int SortPtDistance(const void *ptA, const void *ptB)
{
  double disA = ((fitpt*)ptA)[0].minr;
  double disB = ((fitpt*)ptB)[0].minr;

  if (disA < disB) {
    return -1;
  }
  else if (disA > disB) {
    return 1;
  }
  else {
    return 0;
  }
}

/***=======================================================================***/
/*** ContributeFitPt: contribute this fitting point to the matrix of       ***/
/***                  candidates.                                          ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   grdpt:    pointer to grid point catalog                             ***/
/***   crd:      coordinates                                               ***/
/***   tp:       topology                                                  ***/
/***   fitline:  the next line of the growing fitting matrix to be written ***/
/***=======================================================================***/
static void ContributeFitPt(fitpt *grdpt, coord *crd, double* fitline,
			    fbook *UPot, double wt, int* colmap, int ncol)
{
  int i;
  double dx, dy, dz, r;
  double dijk[3];

  dijk[0] = grdpt->ix;
  dijk[1] = grdpt->iy;
  dijk[2] = grdpt->iz;
  RotateCrd(dijk, 1, UPot->lvec);
  for (i = 0; i < 3; i++) {
    dijk[i] += UPot->orig[i];
  }
  for (i = 0; i < crd->natom; i++) {
    dx = crd->loc[3*i] - dijk[0];
    dy = crd->loc[3*i+1] - dijk[1];
    dz = crd->loc[3*i+2] - dijk[2];
    r = sqrt(dx*dx + dy*dy + dz*dz);

    if (i < 0 || colmap[i] >= ncol) {
      printf("Whoa... i = %4d, maps to %4d out of %4d.\n", i, colmap[i], ncol);
    }

    fitline[colmap[i]] += wt*BIOQ/r;
  }
  fitline[ncol] = wt*UPot->map[grdpt->ix][grdpt->iy][grdpt->iz];
}

/***=======================================================================***/
/*** FlagFitPt: flag all points within a short cutoff around a fitting     ***/
/***            point that was just committed to the matrix of candidate   ***/
/***            fitting points.  These flags, which prevent points from    ***/
/***            being used in the fit, make it so that extremely closely   ***/
/***            positioned points do not both contribute to the fitting    ***/
/***            matrix.                                                    ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   grdpts:   array of candidate fitting points                         ***/
/***   ptid:     ID number of the fitting point within the grdpts catalog  ***/
/***   npt:      the total number of points catalogged in grdpts           ***/
/***   UPot:     the electrostatic potential grid                          ***/
/***   UIdx:     indexes for whether the grid points are flagged or not    ***/
/***   catIdx:   maps the points in UPot/UIdx to their entries in grdpts   ***/
/***=======================================================================***/
static void FlagFitPt(fset *myfit, fitpt* grdpts, int ptid, fbook *UPot,
		      cbook *Uflag, ibook *UIdx, int buffi, int buffj,
		      int buffk)
{
  int i, j, k, imin, imax, jmin, jmax, kmin, kmax, homei, homej, homek;
  double r;
  double dijk[3];
  int *itmp;
  char *ctmp;

  homei = grdpts[ptid].ix;
  homej = grdpts[ptid].iy;
  homek = grdpts[ptid].iz;
  imin = MAX(0, homei - buffi);
  imax = MIN(UPot->pag, homei + buffi);
  jmin = MAX(0, homej - buffj);
  jmax = MIN(UPot->pag, homej + buffj);
  kmin = MAX(0, homek - buffk);
  kmax = MIN(UPot->pag, homek + buffk);
  for (i = imin; i < imax; i++) {
    for (j = jmin; j < jmax; j++) {
      ctmp = Uflag->map[i][j];
      itmp = UIdx->map[i][j];
      for (k = kmin; k < kmax; k++) {
	if (ctmp[k] == 0) {
	  continue;
	}
	dijk[0] = i - homei;
	dijk[1] = j - homej;
	dijk[2] = k - homek;
	RotateCrd(dijk, 1, UPot->lvec);	
	r = sqrt(dijk[0]*dijk[0] + dijk[1]*dijk[1] + dijk[2]*dijk[2]);
	if (r < myfit->flimit) {
	  ctmp[k] = 0;
	  grdpts[itmp[k]].flagged = 1;
	}
      }
    }
  }
}

/***=======================================================================***/
/*** HistogramFitPt: place a fitting point in a histogram based on its     ***/
/***                 distance from the molecule.                           ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   myfit:    the fitting control structure (contains the pre-allocated ***/
/***             histogram)                                                ***/
/***   gpt:      the fitting point                                         ***/
/***=======================================================================***/
static void HistogramFitPt(fset *myfit, fitpt *gpt)
{
  int ir;

  ir = gpt->minr/myfit->fhistbin;
  myfit->fitpthist[ir] += 1;
}

/***=======================================================================***/
/*** WriteConformation: write a conformation of the molecule, for purposes ***/
/***                    of visualizing the charge distribution, in PDB     ***/
/***                    format.                                            ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   h:       the number of the conformation                             ***/
/***   crd:     the coordinates of the conformation                        ***/
/***   tp:      the topology                                               ***/
/***   mfit:    the fitting set                                            ***/
/***   tj:      trajectory control data (directive to overwrite outputs)   ***/
/***=======================================================================***/
static void WriteConformation(int h, double* crd, prmtop *tp, fset *myfit,
			      trajcon *tj)
{
  int i, ires;
  char outname[MAXNAME];
  FILE *outp;

  /*** Only do this for the first instance of each system ***/
  for (i = 0; i < h; i++) {
    if (myfit->tpidx[i] == myfit->tpidx[h]) {
      return;
    }
  }

  /*** Print the conformation iin PDB format ***/
  sprintf(outname, "%s.%s", tp->source, myfit->confext);
  outp = FOpenSafe(outname, tj->OverwriteOutput);
  for (i = 0; i < tp->natom; i++) {
    ires = LocateResID(tp, i, 0, tp->nres);
    fprintf(outp, "ATOM %6d %.4s %.4s%c%4d    %8.3lf%8.3lf%8.3lf\n",
	    i, &tp->AtomNames[4*i], &tp->ResNames[4*ires], 'A', 1,
	    crd[3*i], crd[3*i+1], crd[3*i+2]);
  }
  fprintf(outp, "END\n");
  fclose(outp);
}

/***=======================================================================***/
/*** SelectFitPoint: introduce a criterion to limit the number of points   ***/
/***                 far from the molecular surface which enter the fit.   ***/
/***                 The criterion is that, after a certain cutoff Rc, the ***/
/***                 probability of accepting a point at a distance r from ***/
/***                 the molecular surface drops off such that the number  ***/
/***                 of points accepted at r would be equal to the number  ***/
/***                 at Rc if the molecule were perfectly spherical.       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   gpt:      the grid point                                            ***/
/***   myfit:    the fitting control structure                             ***/
/***   tj:       trajectory control data (contains random number counter)  ***/
/***=======================================================================***/
static int SelectFitPoint(fitpt *gpt, fset *myfit, trajcon *tj)
{
  double r;

  /*** Accept if the point is within Rc of the molecule ***/
  if (gpt->minr < myfit->Rc) {
    return 1;
  }

  /*** Accept a roughly constant number of points for ***/
  /*** any given distance from the molecule beyond Rc ***/
  r = myfit->Rc / gpt->minr;
  if (gpt->minr < myfit->Rmax && ran2(&tj->rndcon) < r*r) {
    return 1;
  }

  /*** Reject the point ***/
  return 0;
}

/***=======================================================================***/
/*** MakeFittingMatrix: compute the fitting matrix based on points taken   ***/
/***                    from all grids, seived by various user-specified   ***/
/***                    cutoffs.                                           ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   myfit:    the fitting control structure                             ***/
/***   Atm2Col:  correspondence between atoms and matrix columns           ***/
/***   allcrd:   matrix to accumulate coordinates of all molecular         ***/
/***             conformations                                             ***/
/***   tj:       trajectory control information (contains random number    ***/
/***             generator counter)                                        ***/
/***=======================================================================***/
static dmat MakeFittingMatrix(fset *myfit, imat *Atm2Col, dmat *allcrd,
			      trajcon *tj)
{
  int h, i, j, k, m, buffi, buffj, buffk, gcon, ngpt, namelen;
  int congruent, setfitpt, totalfitpt;
  long long int totalmem;
  double dx, dy, dz, r, minr;
  double cdepth[3], dijk[3];
  char *ctmp;
  dmat fitmat;
  fitpt* grdpts;
  ibook UIdx;
  cbook Uflag;
  fbook UPot, auxUPot;
  coord crd, auxcrd;
  prmtop *tp;

  /*** Check to ensure that there is sufficient ***/
  /*** room in memory for the fitting matrix    ***/
  totalmem = myfit->ngrd;
  totalmem *= myfit->nfitpt;
  totalmem *= (myfit->q2fit+1) * 24;
  if (totalmem > myfit->MaxMem) {
    printf("MakeFittingMatrix >> Fitting matrix would require %lld bytes of "
	   "memory.\nMakeFittingMatrix >> This exceeds the allowed %lld byte "
	   "limit.  Increase\nMakeFittingMatrix >> the maximum allowed memory "
	   "with the maxmem keyword in the\nMakeFittingMatrix >> &fit "
	   "namelist or reduce the number of fitting points.\n", totalmem,
	   myfit->MaxMem);
    exit(1);
  }
  if (myfit->verbose == 1) {
    printf("mdgx >> Fitting matrices will occupy %lld bytes of memory.\n",
	   totalmem);
    namelen = 0;
    for (i = 0; i < myfit->ngrd; i++) {
      namelen = MAX(namelen, strlen(myfit->gname.map[i]));
    }
  }

  /*** Allocate space for the fitting matrix ***/
  totalfitpt = 0;
  fitmat = CreateDmat(myfit->ngrd*myfit->nfitpt, myfit->q2fit+1, 0);

  /*** Loop over all grids ***/
  for (h = 0; h < myfit->ngrd; h++) {

    /*** Titillate the user ***/
    if (myfit->verbose == 1) {
      fprintf(stderr, "\rmdgx >> Composing fit for %s", myfit->gname.map[h]);
      for (i = 0; i < namelen - strlen(myfit->gname.map[h]); i++) {
	fprintf(stderr, " ");
      }
      fflush(stderr);
    }

    /*** Allocate coordinates for this grid's system ***/
    tp = &myfit->TPbank[myfit->tpidx[h]];
    crd = CreateCoord(tp->natom);
    auxcrd = CreateCoord(tp->natom);

    /*** Order the list of points as a function of distance from   ***/
    /*** the solute.  The list will then be searched in increasing ***/
    /*** order of distance for candidate fitting points.           ***/
    UPot = ReadEPotGrid(myfit->gname.map[h], tp, &crd);
    Uflag = CreateCbook(UPot.pag, UPot.row, UPot.col);
    UIdx = CreateIbook(UPot.pag, UPot.row, UPot.col);
    PrepUPot(&UPot, &crd, tp, &Uflag, myfit, 0);
    if (myfit->auxgname.row == myfit->gname.row) {
      auxUPot = ReadEPotGrid(myfit->auxgname.map[h], tp, &auxcrd);
      PrepUPot(&auxUPot, &auxcrd, tp, &Uflag, myfit, 1);
      congruent = 1;
      for (i = 0; i < 9; i++) {
	if (fabs(auxUPot.lvec.data[i] - UPot.lvec.data[i]) > 1.0e-8) {
	  congruent = 0;
	}
      }
      for (i = 0; i < 3; i++) {
	if (fabs(auxUPot.orig[i] - UPot.orig[i]) > 1.0e-8) {
	  congruent = 0;
	}
      }
      if (congruent == 0) {
	printf("MakeFittingMatrix >> Error.  Non-correspondence in grid "
	       "dimensions for\nMakeFittingMatrix >> %s and %s.\n",
	       myfit->auxgname.map[h], myfit->gname.map[h]);
	exit(1);
      }
      for (i = 0; i < 3*crd.natom; i++) {
	if (fabs(auxcrd.loc[i] - crd.loc[i]) > 1.0e-5) {
	  congruent = 0;
	}
      }
      if (congruent == 0) {
	printf("MakeFittingMatrix >> Error.  Non-correspondence in atom "
	       "positions for\nMakeFittingMatrix >> %s and %s.\n",
	       myfit->auxgname.map[h], myfit->gname.map[h]);
	exit(1);
      }

      /*** If we're still here, this auxiliary grid can be folded       ***/
      /*** into the primary grid; the result can be processed normally. ***/
      for (i = 0; i < UPot.pag*UPot.row*UPot.col; i++) {
	UPot.data[i] = 0.5*(UPot.data[i] + auxUPot.data[i]);
      }

      /*** Free allocated memory ***/
      DestroyFbook(&auxUPot);
    }
    grdpts = (fitpt*)malloc(UPot.pag*UPot.row*UPot.col*sizeof(fitpt));
    gcon = 0;
    for (i = 0; i < UPot.pag; i++) {
      for (j = 0; j < UPot.row; j++) {
	ctmp = Uflag.map[i][j];
	for (k = 0; k < UPot.col; k++) {

	  /*** Skip this point if it is inaccessible ***/
	  if (ctmp[k] == 0) {
	    continue;
	  }
	  dijk[0] = i;
	  dijk[1] = j;
	  dijk[2] = k;
	  RotateCrd(dijk, 1, UPot.lvec);
	  for (m = 0; m < 3; m++) {
	    dijk[m] += UPot.orig[m];
	  }
	  for (m = 0; m < tp->natom; m++) {
	    dx = crd.loc[3*m] - dijk[0];
	    dy = crd.loc[3*m+1] - dijk[1];
	    dz = crd.loc[3*m+2] - dijk[2];
	    if (m == 0) {
	      minr = sqrt(dx*dx + dy*dy + dz*dz);
	    }
	    else {
	      r = sqrt(dx*dx + dy*dy + dz*dz);
	      minr = MIN(r, minr);
	    }
	  }

	  /*** Catalog this grid point ***/
	  grdpts[gcon].ix = i;
	  grdpts[gcon].iy = j;
	  grdpts[gcon].iz = k;
	  grdpts[gcon].flagged = 0;
	  grdpts[gcon].minr = minr;
	  gcon++;
	}
      }
    }
    ngpt = gcon;
    qsort(grdpts, ngpt, sizeof(fitpt), SortPtDistance);

    /*** Label points in the index map based on the new catalog order ***/
    for (i = 0; i < ngpt; i++) {
      UIdx.map[grdpts[i].ix][grdpts[i].iy][grdpts[i].iz] = i;
    }

    /*** Determine the buffer region for testing exclusions ***/
    HessianNorms(&UPot.lvec, cdepth);
    buffi = myfit->prbarm/cdepth[0] + 1;
    buffj = myfit->prbarm/cdepth[1] + 1;
    buffk = myfit->prbarm/cdepth[2] + 1;

    /*** Loop over each grid point, selecting points for fitting  ***/
    /*** when they are accessible.  Selected points will mask out ***/
    /*** others nearby acccording to a minimum distance (flimit)  ***/
    /*** specified in the fitting struct.                         ***/
    setfitpt = 0;
    for (i = 0; i < ngpt; i++) {

      /*** Bail out if the quota for this set is reached ***/
      if (setfitpt == myfit->nfitpt) {
	break;
      }

      /*** Continue if this point is already flagged ***/
      if (grdpts[i].flagged == 1) {
	continue;
      }

      /*** Accept the point with some probability based ***/
      /*** on its distance from the molecular surface   ***/
      if (SelectFitPoint(&grdpts[i], myfit, tj) == 1) {
	ContributeFitPt(&grdpts[i], &crd, fitmat.map[totalfitpt], &UPot,
			myfit->wt[h], Atm2Col->map[myfit->tpidx[h]],
			myfit->q2fit);
	FlagFitPt(myfit, grdpts, i, &UPot, &Uflag, &UIdx, buffi, buffj, buffk);
	HistogramFitPt(myfit, &grdpts[i]);
	totalfitpt++;
	setfitpt++;
      }
      else {

	/*** Even if the point is rejected, the region around ***/
	/*** it must be flagged in order to achieve the point ***/
	/*** spread that is desired.                          ***/
        FlagFitPt(myfit, grdpts, i, &UPot, &Uflag, &UIdx, buffi, buffj, buffk);
      }
    }

    /*** Commit coordinates to buffer for printing if necessary ***/
    for (i = 0; i < 3*tp->natom; i++) {
      allcrd->map[h][i] = crd.loc[i];
    }

    /*** Print coordinates to conformation file ***/
    WriteConformation(h, crd.loc, tp, myfit, tj);

    /*** Free allocated memory ***/
    DestroyFbook(&UPot);
    DestroyIbook(&UIdx);
    DestroyCbook(&Uflag);
    free(grdpts);
    DestroyCoord(&crd);
    DestroyCoord(&auxcrd);
  }

  /*** Titillate the user ***/
  if (myfit->verbose == 1) {
    printf("\nmdgx >> Fitting matrix composed.\n");
  }

  return fitmat;
}

/***=======================================================================***/
/*** SnapFittedQ: snap the fitted set of charges to the nearest 1.0e-5     ***/
/***              proton units, ensure that all charges that were intended ***/
/***              to be equalized are perfectly equalized, and ensure that ***/
/***              the sum of all charges adds to the requested value.      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   q:         the fitted charges                                       ***/
/***   myfit:     fitting control data                                     ***/
/***=======================================================================***/
static void SnapFittedQ(double* q, fset *myfit, imat *Atm2Col)
{
  int i, j, k, chksum, colnum;
  prmtop *tp;

  /*** Round all charges ***/
  for (i = 0; i < myfit->q2fit; i++) {
    if (q[i] >= 0.0) {
      j = q[i]*1.0e5;
      q[i] = j*1.0e-5;
    }
    else {
      j = -q[i]*1.0e5;
      q[i] = -j*1.0e-5;
    }
  }

  /*** Map the fitted charges back to their topologies.   ***/
  /*** In order to make the total charges on each         ***/
  /*** individual topology come out to exactly the stated ***/
  /*** values it may be necessary to slightly break the   ***/
  /*** equality between related charges on different      ***/
  /*** systems, but only by a few hundred-thousandths of  ***/
  /*** a proton charge.                                   ***/
  for (i = 0; i < myfit->tpcount; i++) {
    tp = &myfit->TPbank[i];
    for (j = 0; j < tp->natom; j++) {
      tp->Charges[j] = q[Atm2Col->map[i][j]];
    }
  }

  /*** Adjust charges on all molecules to exactly the specified values ***/
  for (i = 0; i < myfit->tpcount; i++) {
    tp = &myfit->TPbank[i];
    chksum = DSum(tp->Charges, tp->natom) - myfit->totalq[i];
    for (j = 0; j < tp->natom; j++) {
      colnum = Atm2Col->map[i][j];
      for (k = 0; k < tp->natom; k++) {
	if (j != k && Atm2Col->map[i][k] == colnum) {
	  colnum = -1;
	  break;
	}
      }

      /*** If this atom has no equivalents, remove ***/
      /*** any excess charge by adding it here.    ***/
      if (colnum >= 0) {
	tp->Charges[j] -= chksum;
	break;
      }
    }
  }
}

/***=======================================================================***/
/*** SolvRESP: solve a linear least squares problem and report statistics  ***/
/***           on the results.  The fitting matrix is stored in the first  ***/
/***           n-1 columns of the fitting and constraint matrices; the     ***/
/***           target vector is stored in the final column, for simplified ***/
/***           data passing.                                               ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   fitmat:     the fitting matrix and target vector derived from       ***/
/***               electrostatic potential data                            ***/
/***   cnstmat:    the constraints matrix and target vector                ***/
/***   myfit:      the fitting command data, used here to determine the    ***/
/***               number of jackknife fits and their style                ***/
/***=======================================================================***/
static dmat SolveRESP(dmat *fitmat, dmat *cnstmat, fset *myfit, imat *Atm2Col)
{
  int i, j, ip;
  double* bfit;
  double* btst;
  double* bpred;
  dmat Afit, Atst, qres;

  /*** Titillate the user ***/
  if (myfit->verbose == 1) {
    printf("mdgx >> Solving linear least-squares problem for %d independent "
	   "charges.\n", myfit->q2fit);
  }

  /*** Allocate memory for fitting and testing data ***/
  qres = CreateDmat(1, myfit->q2fit+2, 0);
  Afit = CreateDmat(fitmat->row+cnstmat->row, myfit->q2fit, 0);
  bfit = (double*)malloc((fitmat->row+cnstmat->row)*sizeof(double));
  Atst = CreateDmat(fitmat->row, myfit->q2fit, 0);
  btst = (double*)malloc(fitmat->row*sizeof(double));
  bpred = (double*)malloc(fitmat->row*sizeof(double));
  for (i = 0; i < fitmat->row; i++) {
    for (j = 0; j < myfit->q2fit; j++) {
      Afit.map[i][j] = fitmat->map[i][j];
      Atst.map[i][j] = fitmat->map[i][j];
    }
    bfit[i] = fitmat->map[i][myfit->q2fit];
    btst[i] = fitmat->map[i][myfit->q2fit];
  }
  ip = i;
  for (i = 0; i < cnstmat->row; i++) {
    for (j = 0; j < myfit->q2fit; j++) {
      Afit.map[ip][j] = cnstmat->map[i][j];
    }
    bfit[ip] = cnstmat->map[i][myfit->q2fit];
    ip++;
  }
  AxbQRRxc(Afit, bfit, myfit->verbose);
  BackSub(Afit, bfit);
  SnapFittedQ(bfit, myfit, Atm2Col);
  DMatVecMult(&Atst, bfit, bpred);

  /*** Store this set of charges along with its statistics ***/
  ReflectDVec(qres.map[0], bfit, myfit->q2fit);
  qres.map[0][myfit->q2fit] = Pearson(btst, bpred, Atst.row);
  qres.map[0][myfit->q2fit+1] = VecRMSD(btst, bpred, Atst.row);

  /*** Titillate the user ***/
  if (myfit->verbose == 1) {
    printf("mdgx >> Fit complete.\n");
  }

  /*** Free allocated memory ***/
  DestroyDmat(&Afit);
  DestroyDmat(&Atst);
  free(bfit);
  free(btst);
  free(bpred);

  return qres;
}

/***=======================================================================***/
/*** CalculateDipoles: this function takes the best fitted charge set and  ***/
/***                   puts it on each of the fitting conformations to     ***/
/***                   compute dipole moments.  Generates additional text  ***/
/***                   in the output file.                                 ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***                                                                       ***/
/***=======================================================================***/
static void CalculateDipoles(fset *myfit, double *qres, dmat *allcrd,
			     imat *Atm2Col, FILE *outp)
{
  int i, j, k, i3, gnamelen, tpnamelen;
  int *qidx;
  double DPdev;
  double* dp;
  double* sqdp;
  double* sqscr;
  prmtop *tp;

  /*** Allocate memory and compute dipoles ***/
  dp = (double*)calloc(3*myfit->ngrd, sizeof(double));
  sqdp = (double*)malloc(myfit->ngrd*sizeof(double));
  sqscr = (double*)malloc(myfit->ngrd*sizeof(double));
  for (i = 0; i < myfit->ngrd; i++) {

    /*** The dipole is not relevant if the molecule is charged ***/
    if (fabs(myfit->totalq[myfit->tpidx[i]]) > 1.0e-8) {
      continue;
    }
    tp = &myfit->TPbank[myfit->tpidx[i]];
    qidx = Atm2Col->map[myfit->tpidx[i]];
    for (j = 0; j < tp->natom; j++) {
      for (k = 0; k < 3; k++) {
	dp[3*i+k] += qres[qidx[j]]*allcrd->map[i][3*j+k];
      }
    }
  }

  /*** Print outputs ***/
  HorizontalRule(outp, 0);
  fprintf(outp, "(3.) Dipole moments on all conformations, (units of Debye)"
	  "\n\n");
  tpnamelen = 0;
  for (i = 0; i < myfit->tpcount; i++) {
    j = strlen(myfit->TPbank[i].source);
    if (tpnamelen < j) {
      tpnamelen = j;
    }
  }
  for (i = 0; i < myfit->ngrd; i++) {
    i3 = 3*i;
    sqdp[i] = sqrt(dp[i3]*dp[i3] + dp[i3+1]*dp[i3+1] + dp[i3+2]*dp[i3+2]) *
      EA2DEBYE;
  }
  if (myfit->DispAllDP == 1) {
    gnamelen = 0;
    for (i = 0; i < myfit->ngrd; i++) {
      j = strlen(myfit->gname.map[i]);
      if (gnamelen < j) {
	gnamelen = j;
      }
    }
    for (i = 0; i < myfit->ngrd; i++) {
      fprintf(outp, "  %-*s  %-*s  %9.5lf\n", gnamelen, myfit->gname.map[i],
	      tpnamelen, myfit->TPbank[myfit->tpidx[i]].source, sqdp[i]);
    }
    fprintf(outp, "\n");
  }
  fprintf(outp, " System");
  for (i = 0; i < tpnamelen-6; i++) {
    fprintf(outp, " ");
  }
  fprintf(outp, "   Mean Dipole    Std. Dev.\n ");
  for (i = 0; i < tpnamelen; i++) {
    fprintf(outp, "-");
  }
  fprintf(outp, "   -----------   -----------\n");
  for (i = 0; i < myfit->tpcount; i++) {
    k = 0;
    for (j = 0; j < myfit->ngrd; j++) {
      if (myfit->tpidx[j] == i) {
	sqscr[k] = sqdp[j];
	k++;
      }
    }
    DPdev = (k >= 2) ? DStDev(sqscr, k) : 0.0;
    fprintf(outp, " %-*s   %9.5lf   %9.5lf\n", tpnamelen,
	    myfit->TPbank[i].source, DAverage(sqscr, k), DPdev);
  }
  HorizontalRule(outp, 1);

  /*** Free allocated memory ***/
  free(dp);
  free(sqdp);
  free(sqscr);
}

/***=======================================================================***/
/*** PrintEPRules: print a file of extra point rules that will modify each ***/
/***               topology and add any needed extra points.               ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   myfit:      fitting control data                                    ***/
/***   tpnum:      the number of the topology for which to print rules     ***/
/***               (the fitted charges for each topology are now stored in ***/
/***               the topology structs themselves)                        ***/
/***   tj:         trajectory control information                          ***/
/***=======================================================================***/
static void PrintEPRules(fset *myfit, int tpnum, trajcon *tj)
{
  int i, j, natm, nres;
  eprule *tmr;
  char outname[MAXNAME];
  FILE *outp;
  prmtop *tp;

  tp = &myfit->TPbank[tpnum];
  sprintf(outname, "%s.%s", tp->source, myfit->epext);
  outp = FOpenSafe(outname, tj->OverwriteOutput);
  fprintf(outp, "%% Extra points rules for charge model fitted as described "
	  "in\n%% %s.\n%%\n\n", tj->outbase);
  fprintf(outp, "%% The following rules will modify the charges of\n%% atoms "
	  "given in the topology %s.\n", tp->source);
  for (i = 0; i < tp->natom; i++) {
    nres = LocateResID(tp, i, 0, tp->nres);
    if (tp->EPInserted == 0 || tp->OldAtomNum[i] >= 0) {
      fprintf(outp, "&rule\n  ResidueName  %.4s\n  ExtraPoint   %.4s\n  "
	      "FrameStyle   0\n  Charge       %9.5lf\n&end\n\n",
	      &tp->ResNames[4*nres], &tp->AtomNames[4*i],
	      tp->Charges[i]);
    }
  }
  if (tp->EPInserted == 0) {
    return;
  }
  fprintf(outp, "%% The following rules will add charged extra points\n%% "
	  "to the original topology.\n");
  for (i = 0; i < tp->neprule; i++) {
    tmr = &tp->eprules[i];
    for (j = 0; j < tp->natom; j++) {
      if (strncmp(&tp->AtomNames[4*j], tmr->epname, 4) == 0) {
	natm = j;
      }
    }
    fprintf(outp, "&rule\n  ResidueName  %.4s\n  ExtraPoint   %.4s\n  "
	     "FrameStyle   %d\n", tmr->resname, tmr->epname, tmr->frstyle);
    fprintf(outp, "  FrameAtom1   %.4s\n  FrameAtom2   %.4s\n", tmr->fr1,
	    tmr->fr2);
    if (tmr->frstyle > 1) {
      fprintf(outp, "  FrameAtom3   %.4s\n", tmr->fr3);
    }
    if (tmr->frstyle == 6) {
      fprintf(outp, "  FrameAtom4   %.4s\n", tmr->fr4);
    }
    fprintf(outp, "  Vector12     %16.10lf\n", tmr->d1);
    if (tmr->frstyle == 2 || tmr->frstyle == 5 || tmr->frstyle == 6) {
      fprintf(outp, "  Vector13     %16.10lf\n", tmr->d2);
    }
    else if (tmr->frstyle == 4) {
      fprintf(outp, "  Theta        %16.10lf\n", tmr->d2);
    }
    if (tmr->frstyle == 3) {
      fprintf(outp, "  Vector23     %16.10lf\n", tmr->d3);
    }
    else if (tmr->frstyle == 5) {
      fprintf(outp, "  Vector12x13  %16.10lf\n", tmr->d3);
    }
    fprintf(outp, "  Charge       %9.5lf\n", tp->Charges[natm]);
    if (tmr->sig >= 0.0) {
      fprintf(outp, "  Sigma        %16.10lf\n", tmr->sig);
    }
    if (tmr->eps >= 0.0) {
      fprintf(outp, "  Epsilon      %16.10lf\n", tmr->eps);
    }
    fprintf(outp, "&end\n\n");
  }
  fclose(outp);
}

/***=======================================================================***/
/*** PrintRespHistogram: print a histogram of the RESP fitting points,     ***/
/***                     indicating their minimum distance to the solute.  ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   myfit:      the fitting command data, used here to determine the    ***/
/***               number of jackknife fits and their style                ***/
/***=======================================================================***/
static void PrintRespHistogram(fset *myfit, trajcon *tj)
{
  int i;
  FILE *outp;

  outp = FOpenSafe(myfit->histfile, tj->OverwriteOutput);
  fprintf(outp, "%% Molecule:fit point distance histogram for charge model "
	  "fitted\n%% as described in %s.\n%%\n%% Counts are normalized by "
	  "the total number of fitting points.\n\n%% Bin Center     Count\n"
	  "%% ----------   ---------\n", tj->outbase);
  for (i = 0; i < 10.0/myfit->fhistbin; i++) {
    fprintf(outp, "  %10.6lf   %9.6lf\n", (i+0.5)*myfit->fhistbin,
	    (double)myfit->fitpthist[i]/(myfit->ngrd*myfit->nfitpt));
  }
  fclose(outp);
}

/***=======================================================================***/
/*** OutputRESP: produce formatted output to present the results of the    ***/
/***             RESP calculation.  If an extra points file was specified, ***/
/***             then a new extra points file with the appropriate names   ***/
/***             and charges assigned to each point will be written.  A    ***/
/***             set of scaled charges will be written in 5 x %16.8e       ***/
/***             format for inclusion in the topology file.                ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   myfit:      the fitting command data, used here to determine the    ***/
/***               number of jackknife fits and their style                ***/
/***   tj:         trajectory control data, used here for directives on    ***/
/***               output overwriting and input file reprinting            ***/
/***   qres:       matrix filled with putative solutions to the RESP fit   ***/
/***   allcrd:     coordinates of all fitting conformations                ***/
/***   Atm2Col:    correspondence of atoms in each system and columns in   ***/
/***               the fitting matrix                                      ***/
/***=======================================================================***/
static void OutputRESP(fset *myfit, trajcon *tj, dmat *qres, dmat *allcrd,
		       imat *Atm2Col)
{
  int i, j, k, minloc;
  double maxcorr, minerr, qvar;
  double* qvalues;
  FILE *outp;
  time_t ct;
  prmtop *tp;

  /*** Open the output file ***/
  ct = time(&ct);
  outp = FOpenSafe(tj->outbase, tj->OverwriteOutput);
  fprintf(outp, "Run on %s", asctime(localtime(&ct)));

  /*** Reprint the input file ***/
  HorizontalRule(outp, 0);
  fprintf(outp, "\nINPUT LINE TEXT:\n\n");
  PrintParagraph(tj->inpline, 79, outp);
  fprintf(outp, "\nINPUT FILE TEXT:\n\n");
  for (i = 0; i < tj->inptext.row; i++) {
    fprintf(outp, "%s", tj->inptext.map[i]);
  }
  HorizontalRule(outp, 1);

  /*** Print the best resulting charge set(s) ***/
  HorizontalRule(outp, 0);
  fprintf(outp, "(1.) Charges on all atoms, proton units\n\n");
  for (i = 0; i < qres->row; i++) {
    if (i == 0 || qres->map[i][myfit->q2fit+1] < minerr) {
      minerr = qres->map[i][myfit->q2fit+1];
      maxcorr = qres->map[i][myfit->q2fit];
      minloc = i;
    }
  }
  fprintf(outp, " Correlation of fitted to original electrostatic "
	  "potential: %8.5lf\n", maxcorr);
  fprintf(outp, " Error fitted versus original electrostatic "
	  "potential (kcal/mol): %8.5lf\n", minerr);
  for (i = 0; i < myfit->tpcount; i++) {
    fprintf(outp, "\n System %d: %s\n Atom    Charge   Variance\n ----  "
	    "----------  ----------\n", i+1, myfit->TPbank[i].source);
    tp = &myfit->TPbank[i];
    if (qres->row == 1) {
      qvar = 0.0;
    }
    else {
      qvalues = (double*)malloc(qres->row*sizeof(double));
      for (j = 0; j < qres->row; j++) {
	qvalues[j] = qres->map[j][i];
      }
      qvar = DStDev(qvalues, qres->row);
      free(qvalues);
    }
    for (j = 0; j < tp->natom; j++) {
      fprintf(outp, " %.4s  %10.5lf  %10.5lf\n", &tp->AtomNames[4*j],
	      qres->map[minloc][Atm2Col->map[i][j]], qvar);
    }
  }
  HorizontalRule(outp, 1);

  /*** Print the best resulting charge set in prmtop format ***/
  HorizontalRule(outp, 0);
  fprintf(outp, "(2.) Charges on all atoms, prmtop format (proton units "
	  "scaled by 18.2223)\n");
  for (i = 0; i < myfit->tpcount; i++) {
    tp = &myfit->TPbank[i];
    fprintf(outp, "\n System %d: %s\n", i+1, myfit->TPbank[i].source);
    k = 0;
    for (j = 0; j < tp->natom; j++) {
      fprintf(outp, "%16.8e",
	      qres->map[minloc][Atm2Col->map[i][j]]*sqrt(BIOQ));
      k++;
      if (k == 5) {
	k = 0;
	fprintf(outp, "\n");
      }
    }
    if (k > 0) {
      fprintf(outp, "\n");
    }
  }
  HorizontalRule(outp, 1);

  /*** Compute and print dipole moments ***/
  CalculateDipoles(myfit, qres->map[minloc], allcrd, Atm2Col, outp);
  fclose(outp);

  /*** Print an extra points file ***/
  if (myfit->epext[0] != '\0') {
    for (i = 0; i < myfit->tpcount; i++) {
      PrintEPRules(myfit, i, tj);
    }
  }

  /*** Print the histogram of fitting points ***/
  if (myfit->histfile[0] != '\0') {
    PrintRespHistogram(myfit, tj);
  }
}

/***=======================================================================***/
/*** DestroyFitData: cleanup for all data related to charge fitting.       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   myfit:  the fitting control data                                    ***/
/***=======================================================================***/
static void DestroyFitData(fset *myfit)
{
  int i;

  /*** Simple frees ***/
  free(myfit->fitpthist);
  free(myfit->tpidx);
  free(myfit->totalq);
  free(myfit->wt);

  /*** Constraint structs ***/
  for (i = 0; i < myfit->nqeq; i++) {
    free(myfit->qeq[i].atoms);
    free(myfit->qeq[i].maskstr);
  }
  free(myfit->qeq);
  for (i = 0; i < myfit->nqmin; i++) {
    free(myfit->qmin[i].maskstr);
  }
  free(myfit->qmin);
  for (i = 0; i < myfit->nqsum; i++) {
    free(myfit->qsum[i].maskstr);
  }
  free(myfit->qsum);

  /*** File names ***/
  DestroyCmat(&myfit->gname);
  DestroyCmat(&myfit->auxgname);
  DestroyCmat(&myfit->tpname);
  DestroyCmat(&myfit->eprule);

  /*** Topologies ***/
  for (i = 0; i < myfit->tpcount; i++) {
    FreeTopology(&myfit->TPbank[i]);
  }
  free(myfit->TPbank);
}

/***=======================================================================***/
/*** FitCharges: the main function for parameter optimization.             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tp:     the system topology (the initial guess in the case of full  ***/
/***           nonlinear parameter optimization)                           ***/
/***   tj:     trajectory control information                              ***/
/***   myfit:  the fitting control data                                    ***/
/***=======================================================================***/
void FitCharges(prmtop *tp, trajcon *tj, fset *myfit)
{
  int i, maxatm;
  imat Atm2Col;
  dmat cnstmat, fitmat, qres;
  dmat allcrd;

  /*** Match up topology names with grids ***/
#ifdef MPI
  AssignGridTopologies(myfit, tp, tj);
#else
  AssignGridTopologies(myfit, tp);
#endif

  /*** Match up atoms in all topologies with matrix columns ***/
  Atm2Col = ColumnsForAtoms(myfit);

  /*** Allocate space for coordinate buffer ***/
  maxatm = 0;
  for (i = 0; i < myfit->tpcount; i++) {
    if (myfit->TPbank[i].natom > maxatm) {
      maxatm = myfit->TPbank[i].natom;
    }
  }
  allcrd = CreateDmat(myfit->ngrd, 3*maxatm, 0);

  /*** Create the matrix of all fitting points ***/
  fitmat = MakeFittingMatrix(myfit, &Atm2Col, &allcrd, tj);

  /*** Create the matrix of constraints ***/
  cnstmat = MakeCnstMatrix(myfit, &Atm2Col, &fitmat);

  /*** Solve the linear least squares problem ***/
  qres = SolveRESP(&fitmat, &cnstmat, myfit, &Atm2Col);

  /*** Output results ***/
  OutputRESP(myfit, tj, &qres, &allcrd, &Atm2Col);

  /*** Free allocated memory ***/
  DestroyDmat(&fitmat);
  DestroyDmat(&cnstmat);
  DestroyDmat(&allcrd);
  DestroyImat(&Atm2Col);
  DestroyFitData(myfit);

  /*** Exit!  We are done. ***/
  exit(1);
}
