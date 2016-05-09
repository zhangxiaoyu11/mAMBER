#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mdgxVector.h"
#include "Matrix.h"
#include "ParamFit.h"
#include "Topology.h"
#include "CrdManip.h"
#include "Parse.h"
#include "ChargeFit.h"
#include "VirtualSites.h"
#include "Manual.h"
#include "Nonbonded.h"
#include "Trajectory.h"
#include "Macros.h"
#include "ParamOut.h"
#include "ParamRead.h"

#include "CrdManipDS.h"

/***=======================================================================***/
/*** TypeCompare: compare the types of two pairs of atoms, returning 2 if  ***/
/***              the types can be connected directly, 1 if they can be    ***/
/***              connected by a wildcard, and 0 if they cannot connect.   ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   T[1,2][a-d]:   the first and second sets' various types a-d         ***/
/***   order:         the order of the comparison (2 = bond, 3 = angle,    ***/
/***                  4 = torsion)                                         ***/
/***   strict:        if set to 1, this will require that the atom type    ***/
/***                  names strictly match in either forward or reverse    ***/
/***                  order; no wildcards!                                 ***/
/***=======================================================================***/
int TypeCompare(char* T1a, char* T1b, char* T1c, char* T1d, char* T2a,
		char* T2b, char* T2c, char* T2d, int order, int strict)
{
  int aworks, bworks, cworks, dworks;

  /*** Possible problem with this routine: cannot interpret   ***/
  /*** atom type sequences with wildcards in middle positions ***/
  /*** like CT-X -O  or CT-X -CT-H1, if such things exist     ***/
  if (order == 2 || order == 3) {

    /*** Angles are compared by supplying the atom types ***/
    /*** to this function in the order A, C, B           ***/
    if (order == 3 && strncmp(T1c, T2c, 4) != 0) {
      return 0;
    }
    if ((strncmp(T1a, T2a, 4) == 0 && strncmp(T1b, T2b, 4) == 0) ||
	(strncmp(T1a, T2b, 4) == 0 && strncmp(T1b, T2a, 4) == 0)) {
      return 2;
    }
    if (strict == 1) {
      return 0;
    }
    if ((strncmp(T1a, T2a, 4) == 0 || strncmp(T2a, "X   ", 4) == 0) && 
	(strncmp(T1b, T2b, 4) == 0 || strncmp(T2b, "X   ", 4) == 0)) {
      return 1;
    }
    if ((strncmp(T1a, T2b, 4) == 0 || strncmp(T2b, "X   ", 4) == 0) && 
	(strncmp(T1b, T2a, 4) == 0 || strncmp(T2a, "X   ", 4) == 0)) {
      return 1;
    }
    if (strncmp(T2a, "X   ", 4) == 0 && strncmp(T2b, "X   ", 4) == 0) {
      return 1;
    }
  }
  if (order == 4) {

    /*** Check for explicit correspondence in the types ***/
    if ((strncmp(T1a, T2a, 4) == 0 && strncmp(T1b, T2b, 4) == 0 &&
	 strncmp(T1c, T2c, 4) == 0 && strncmp(T1d, T2d, 4) == 0) ||
	(strncmp(T1a, T2d, 4) == 0 && strncmp(T1b, T2c, 4) == 0 &&
         strncmp(T1c, T2b, 4) == 0 && strncmp(T1d, T2a, 4) == 0)) {
      return 2;
    }
    if (strict == 1) {
      return 0;
    }

    /*** Check for parallel correspondence, with wildcards ***/
    aworks = (strncmp(T2a, "X   ", 4) == 0 ||
	      strncmp(T1a, T2a, 4) == 0) ? 1 : 0;
    bworks = (strncmp(T2b, "X   ", 4) == 0 ||
              strncmp(T1b, T2b, 4) == 0) ? 1 : 0;
    cworks = (strncmp(T2c, "X   ", 4) == 0 ||
              strncmp(T1c, T2c, 4) == 0) ? 1 : 0;
    dworks = (strncmp(T2d, "X   ", 4) == 0 ||
              strncmp(T1d, T2d, 4) == 0) ? 1 : 0;
    if (aworks ==  1 && bworks == 1 && cworks == 1 && dworks == 1) {
      return 1;
    }

    /*** Check for anti-parallel correspondence, with wildcards ***/
    dworks = (strncmp(T2a, "X   ", 4) == 0 ||
	      strncmp(T1d, T2a, 4) == 0) ? 1 : 0;
    cworks = (strncmp(T2b, "X   ", 4) == 0 ||
              strncmp(T1c, T2b, 4) == 0) ? 1 : 0;
    bworks = (strncmp(T2c, "X   ", 4) == 0 ||
              strncmp(T1b, T2c, 4) == 0) ? 1 : 0;
    aworks = (strncmp(T2d, "X   ", 4) == 0 ||
              strncmp(T1a, T2d, 4) == 0) ? 1 : 0;
    if (aworks ==  1 && bworks == 1 && cworks == 1 && dworks == 1) {
      return 1;
    }
  }

  return 0;
}

/***=======================================================================***/
/*** PruneDuplicateAdj: prune duplicate adjustable terms.                  ***/
/***=======================================================================***/
static void PruneDuplicateAdj(prmset *mp, int order)
{
  int i, j, nadj, nuniadj;
  int* duplicate;

  nadj = (order == 2) ? mp->nbadj : (order == 3) ? mp->naadj : mp->nhadj;
  if (nadj <= 1) {
    return;
  }

  /*** Detect any duplicate specifications ***/
  duplicate = (int*)calloc(nadj, sizeof(int));
  for (i = 0; i < nadj-1; i++) {
    if (duplicate[i] > 0) {
      continue;
    }
    for (j = i+1; j < nadj; j++) {
      if (duplicate[j] > 0) {
	continue;
      }
      if (order == 2) {
	duplicate[j] = TypeCompare(mp->badj[i].atype, mp->badj[i].btype,
				   "    ", "    ", mp->badj[j].atype,
				   mp->badj[j].btype, "    ", "    ", 2, 1);
      }
      else if (order == 3) {
	duplicate[j] = TypeCompare(mp->aadj[i].atype, mp->aadj[i].btype,
				   mp->aadj[i].ctype, "    ",
				   mp->aadj[j].atype, mp->aadj[j].btype,
				   mp->aadj[j].ctype, "    ", 3, 1);
      }
      else {
	duplicate[j] = TypeCompare(mp->hadj[i].atype, mp->hadj[i].btype,
				   mp->hadj[i].ctype, mp->hadj[i].dtype,
				   mp->hadj[j].atype, mp->hadj[j].btype,
				   mp->hadj[j].ctype, mp->hadj[j].dtype,
				   4, 1);
      }
    }
  }

  /*** Preserve only the unique parameter specifications ***/
  nuniadj = 0;
  for (i = 0; i < nadj; i++) {
    if (duplicate[i] == 0) {
      if (order == 2) {
	mp->badj[nuniadj] = mp->badj[i];
      }
      else if (order == 3) {
	mp->aadj[nuniadj] = mp->aadj[i];
      }
      else {
	mp->hadj[nuniadj] = mp->hadj[i];
      }
      nuniadj++;
    }
  }
  if (order == 2) {
    mp->nbadj = nuniadj;
  }
  else if (order == 3) {
    mp->naadj = nuniadj;
  }
  else {
    mp->nhadj = nuniadj;
  }

  /*** Free allocated memory ***/
  free(duplicate);
}

/***=======================================================================***/
/*** CleanCommandInput: clean up command input information, to avoid some  ***/
/***                    problems later in the fit.  Removes duplicate      ***/
/***                    specifications of adjustable bond, angle, and      ***/
/***                    dihedral terms.                                    ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:      the parameter fitting data                                 ***/
/***=======================================================================***/
static void CleanCommandInput(prmset *mp)
{

  /*** Check for duplicate adjustable term specifications ***/
  PruneDuplicateAdj(mp, 2);
  PruneDuplicateAdj(mp, 3);
  PruneDuplicateAdj(mp, 4);
}

/***=======================================================================***/
/*** ReadSystemConf: read a system conformation and store it in memory,    ***/
/***                 placing extra points if needed.                       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   conf:    the conformation being read                                ***/
/***   n:       the number of the conformation in the master parameter set ***/
/***   mp:      the master parameter set                                   ***/
/***=======================================================================***/
#ifdef MPI
static void ReadSystemConf(mmsys *conf, int n, prmset *mp, trajcon *tj)
#else
static void ReadSystemConf(mmsys *conf, int n, prmset *mp)
#endif
{
  /*** Augment the topology with universal settings ***/
  sprintf(conf->tp.WaterName, mp->WaterName);
  conf->tp.ljbuck = mp->ljbuck;
  if (mp->ep[0] != '\0') {
    sprintf(conf->tp.eprulesource, mp->ep);
  }
  else {
    conf->tp.eprulesource[0] = '\0';
  }

  /*** Settle and Rattle are disabled in      ***/
  /*** parameter fitting.  Nothing is moving. ***/ 
  conf->tp.settle = 0;
  conf->tp.rattle = 0;
  conf->tp.lj14fac = mp->lj14fac;
  conf->tp.elec14fac = mp->elec14fac;

  /*** Read the topology ***/
#ifdef MPI
  if (n > 0 && strcmp(mp->conf[n-1].tp.source, conf->tp.source) == 0) {
    conf->tp = CopyTopology(&mp->conf[n-1].tp);
  }
  else {
    GetPrmTop(&conf->tp, tj, 1);
  }
#else
  if (n > 0 && strcmp(mp->conf[n-1].tp.source, conf->tp.source) == 0) {
    conf->tp = CopyTopology(&mp->conf[n-1].tp);
  }
  GetPrmTop(&conf->tp, 1);
#endif

  /*** Placeholders for certain topology fields;    ***/
  /*** they are not used in this module but do need ***/
  /*** to be allocated so they can later be freed.  ***/
  conf->tp.lVDWc = (double*)calloc(conf->tp.ntypes, sizeof(double));
  conf->tp.rattlemask = (char*)calloc(MAXNAME, sizeof(char));
  conf->tp.norattlemask = (char*)calloc(MAXNAME, sizeof(char));

  /*** Read the coordinates ***/
  conf->crd = ReadRst(&conf->tp, conf->crdsrc);
}

/***=======================================================================***/
/*** FreeConformation: free all memory associated with a conformation.     ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   conf:   the conformation to free                                    ***/
/***=======================================================================***/
static void FreeConformation(mmsys *conf)
{
  FreeTopology(&conf->tp);
  DestroyCoord(&conf->crd);
  free(conf->bmap.id);
  free(conf->bmap.val);
  free(conf->bmap.kernel);
  free(conf->bmap.contrib);
  free(conf->amap.id);
  free(conf->amap.val);
  free(conf->amap.kernel);
  free(conf->amap.contrib);
  free(conf->hmap.id);
  free(conf->hmap.val);
  free(conf->hmap.kernel);
  free(conf->hmap.contrib);
  DestroyDmat(&conf->excl);
  DestroyDmat(&conf->nbnrg);
}

/***=======================================================================***/
/*** GroupSystems: group the systems according to similar topology files,  ***/
/***               and compute enorm for each system.  This function also  ***/
/***               converts the energies to internal units of kcal/mol, if ***/
/***               they are not already in that format.                    ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:      the fitting data (contains a list of all systems)          ***/
/***=======================================================================***/
static void GroupSystems(prmset *mp)
{
  int i, j, ngrp, grpsize;

  /*** Determine the energy units ***/
  j = strlen(mp->NrgUnits);
  for (i = 0; i < j; i++) {
    mp->NrgUnits[i] = ToUpper(mp->NrgUnits[i]);
  }
  if (strcmp(mp->NrgUnits, "HARTREE") == 0 ||
      strcmp(mp->NrgUnits, "ATOMIC") == 0) {
    for (i = 0; i < mp->nconf; i++) {
      mp->conf[i].etrg *= 627.509469;
    }
  }
  else if (strcmp(mp->NrgUnits, "KJ") == 0 ||
	   strcmp(mp->NrgUnits, "KILOJOULES") == 0) {
    for (i = 0; i < mp->nconf; i++) {
      mp->conf[i].etrg /= 4.184;
    }
  }
  else if (strcmp(mp->NrgUnits, "J") == 0 ||
	   strcmp(mp->NrgUnits, "JOULES") == 0) {
    for (i = 0; i < mp->nconf; i++) {
      mp->conf[i].etrg /= 4184.0;
    }
  }

  /*** Enumerate all topology groups ***/
  for (i = 0; i < mp->nconf; i++) {
    mp->conf[i].GroupNum = -1;
  }
  ngrp = 0;
  for (i = 0; i < mp->nconf; i++) {
    if (mp->conf[i].GroupNum < 0) {
      mp->conf[i].GroupNum = ngrp;
      grpsize = 1;
      for (j = i+1; j < mp->nconf; j++) {
	if (strcmp(mp->conf[i].tp.source, mp->conf[j].tp.source) == 0) {
	  mp->conf[j].GroupNum = ngrp;
	  grpsize++;
	}
      }
      ngrp++;
    }
  }
  mp->nunisys = ngrp;
  mp->GroupCount = (int*)calloc(ngrp, sizeof(int));
  for (i = 0; i < mp->nconf; i++) {
    mp->GroupCount[mp->conf[i].GroupNum] += 1;
  }
}

/***=======================================================================***/
/*** CompEnorm: compute enorm for all systems, a normalized target energy  ***/
/***            adjusted to put the average of all target energies for     ***/
/***            each system at ther average molecular mechanics energy     ***/
/***            according to the original model.                           ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:     the fitting data set                                        ***/
/***=======================================================================***/
static void CompEnorm(prmset *mp)
{
  int i, j, k, ngrp, maxsize;
  int* grpsize;
  double gspread, gmean;
  double* eave;
  double* eadj;
  imat allid;
  dmat allval;

  eave = (double*)calloc(mp->nunisys, sizeof(double));
  eadj = (double*)calloc(mp->nunisys, sizeof(double));
  grpsize = (int*)calloc(mp->nunisys, sizeof(int));
  for (i = 0; i < mp->nconf; i++) {
    ngrp = mp->conf[i].GroupNum;
    eave[ngrp] += mp->conf[i].etrg;
    eadj[ngrp] += mp->conf[i].eorig;
    grpsize[ngrp] += 1;
  }
  maxsize = 0;
  for (i = 0; i < mp->nunisys; i++) {
    eave[i] /= grpsize[i];
    eadj[i] /= grpsize[i];
    maxsize = (grpsize[i] > maxsize) ? grpsize[i] : maxsize;
    grpsize[i] = 0;
  }

  /*** Compute enorm for all groups ***/
  for (i = 0; i < mp->nconf; i++) {
    ngrp = mp->conf[i].GroupNum;
    mp->conf[i].enorm = mp->conf[i].etrg - eave[ngrp] + eadj[ngrp];
  }

  /*** Check for outliers and warn if they exist ***/
  allval = CreateDmat(mp->nunisys, maxsize, 0);
  allid = CreateImat(mp->nunisys, maxsize);
  for (i = 0; i < mp->nconf; i++) {
    ngrp = mp->conf[i].GroupNum;
    allval.map[ngrp][grpsize[ngrp]] = mp->conf[i].enorm;
    grpsize[ngrp] += 1;
    allid.map[ngrp][grpsize[ngrp]] = i;
  }
  for (i = 0; i < mp->nunisys; i++) {
    gspread = DStDev(allval.map[i], grpsize[i]);
    gmean = DAverage(allval.map[i], grpsize[i]);
    for (j = 0; j < grpsize[i]; j++) {
      if (fabs(allval.map[i][j] - gmean) / gspread > mp->esigtol) {
	if (mp->verbose == 1) {
	  printf("CompENorm >> Warning.  Conformation %6d (system %s)\n"
		 "CompENorm >> energy is %7.4lf sigma from the mean.\n"
		 "CompEnorm >> Sigma = %9.4lf, target energy = %9.4lf (mean "
		 "%9.4lf)\n\n", allid.map[i][j],
		 mp->conf[allid.map[i][j]].tp.source,
		 fabs(allval.map[i][j]-gmean) / gspread, gspread,
		 allval.map[i][j], gmean);
	}

	/*** Destroy this offending conformation ***/
	if (mp->RemoveOutliers == 1) {
	  FreeConformation(&mp->conf[allid.map[i][j]]);
	  for (k = allid.map[i][j]; k < mp->nconf-1; k++) {
	    mp->conf[k] = mp->conf[k+1];
	  }
	  for (k = 0; k < mp->nunisys*maxsize; k++) {
	    if (allid.data[k] > allid.map[i][j]) {
	      allid.data[k] -= 1;
	    }
	  }
	  mp->nconf -= 1;
	}
      }
    }
  }
}

/***=======================================================================***/
/*** MakeTermKey: make a table of unique bonded interactions across all    ***/
/***              systems, and in so doing create a table of bonded terms  ***/
/***              for each system with indexing into the master key.       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:     the set of systems and parameter guidelines                 ***/
/***   order:  the order of the terms, 2 for bonds, 3 for angles, and 4    ***/
/***           for dihedral terms                                          ***/
/***=======================================================================***/
static void MakeTermKey(prmset *mp, int order)
{
  int i, j, k, m, n, klim, nmpterm, nterm, found;
  int atmA, atmB, atmC, atmD, resA, resB, resC, resD;
  int compval, maxcomp;
  double diheK, dihePhi, diheN, singlet;
  char typeA[8], typeB[8], typeC[8], typeD[8];
  prmtop *tp;
  coord *crd;

  /*** Allocate space for bonds in the fitting set's master key ***/
  if (order == 2) {
    nmpterm = mp->nbond;
    for (i = 0; i < nmpterm; i++) {
      mp->bonds[i].ninst = 0;
      mp->bonds[i].fitcol = -1;
    }
  }
  else if (order == 3) {
    nmpterm = mp->nangl;
    for (i = 0; i < nmpterm; i++) {
      mp->angls[i].ninst = 0;
      mp->angls[i].fitcol = -1;
    }
  }
  else if (order == 4) {
    nmpterm = mp->ntor;
    for (i = 0; i < nmpterm; i++) {
      mp->torsions[i].ninst = 0;
      mp->torsions[i].fitcol = -1;
    }
  }

  /*** Loop over all conformations ***/
  for (i = 0; i < mp->nconf; i++) {

    /*** Unpack this system ***/
    tp = &mp->conf[i].tp;
    crd = &mp->conf[i].crd;

    /*** Allocate the bonded terms table for this system ***/
    nterm = 0;
    if (order == 2) {
      for (j = 0; j < tp->natom; j++) {
	nterm += tp->BLC[j].nbond;
      }
      mp->conf[i].bmap.id = (bidx*)malloc(nterm*sizeof(bidx));
      mp->conf[i].bmap.nbond = nterm;
    }
    else if (order == 3) {
      for (j = 0; j < tp->natom; j++) {
	nterm += tp->ALC[j].nangl;
      }
      mp->conf[i].amap.id = (aidx*)malloc(nterm*sizeof(aidx));
      mp->conf[i].amap.nangl = nterm;
    }
    if (order == 4) {
      for (j = 0; j < tp->natom; j++) {
	for (k = 0; k < tp->HLC[j].ndihe; k++) {
	  nterm += tp->HLC[j].HC[k].nt;
	}
      }
      mp->conf[i].hmap.id = (hidx*)malloc(nterm*sizeof(hidx));
      mp->conf[i].hmap.ntterm = nterm;
    }

    /*** Examine each bonded term individually ***/
    nterm = 0;
    for (j = 0; j < tp->natom; j++) {
      klim = (order == 2) ? tp->BLC[j].nbond :
	(order == 3) ? tp->ALC[j].nangl : tp->HLC[j].ndihe;
      for (k = 0; k < klim; k++) {
	if (order == 2) {
	  atmA = tp->BLC[j].BC[k].a;
	  atmB = tp->BLC[j].BC[k].b;
	}
	else if (order == 3) {
	  atmA = tp->ALC[j].AC[k].a;
	  atmB = j;
	  atmC = tp->ALC[j].AC[k].c;
	}
	else if (order == 4) {
	  atmA = tp->HLC[j].HC[k].a;
	  atmB = j;
	  atmC = tp->HLC[j].HC[k].c;
	  atmD = tp->HLC[j].HC[k].d;
	}
	strncpy(typeA, &tp->AtomTypes[4*atmA], 4);
	strncpy(typeB, &tp->AtomTypes[4*atmB], 4);
	if (order > 2) {
	  strncpy(typeC, &tp->AtomTypes[4*atmC], 4);
	}
	if (order == 4) {
	  strncpy(typeD, &tp->AtomTypes[4*atmD], 4);
	}

	/*** Reference this bond against the master key ***/
	found = 0;
	maxcomp = 0;
	if (order == 2) {
	  for (m = 0; m < nmpterm; m++) {
	    compval = TypeCompare(typeA, typeB, "    ", "    ",
				  mp->bonds[m].atype, mp->bonds[m].btype,
				  "    ", "    ", 2, 0);
	    if (compval > maxcomp) {

	      /*** Undo the previous assignment if one exists ***/
	      if (maxcomp == 1) {
		mp->bonds[mp->conf[i].bmap.id[nterm].key].ninst -= 1;
	      }
	      maxcomp = compval;

	      /*** Assign or reassign the bond to its place in the key ***/
	      mp->conf[i].bmap.id[nterm].key = m;
	      mp->bonds[m].ninst += 1;
	      found = 1;
	    }
	  }
	  if (found == 0) {

	    /*** This bond is not present in the parameter set ***/
	    printf("MakeTermKey >> Error.  A bond between %.4s and %.4s is "
		   "not found in the\nMakeTermKey >> parameter files.\n",
		   typeA, typeB);
	    exit(1);
	  }
	  mp->conf[i].bmap.id[nterm].a = atmA;
	  mp->conf[i].bmap.id[nterm].b = atmB;
	  nterm++;
	}
	else if (order == 3) {
	  for (m = 0; m < nmpterm; m++) {
	    compval = TypeCompare(typeA, typeC, typeB, "    ",
			 	  mp->angls[m].atype, mp->angls[m].ctype,
				  mp->angls[m].btype, "    ", 3, 0);
	    if (compval > maxcomp) {

              /*** Undo the previous assignment if one exists ***/
              if (maxcomp == 1) {
                mp->angls[mp->conf[i].amap.id[nterm].key].ninst -= 1;
              }
              maxcomp = compval;

              /*** Assign or reassign the angle to its place in the key ***/
	      mp->conf[i].amap.id[nterm].key = m;
	      mp->angls[m].ninst += 1;
	      found = 1;
	    }
	  }
	  if (found == 0) {

	    /*** This angle is not present in the parameter set ***/
	    printf("MakeTermKey >> Error.  An angle between %.4s, %.4s, and "
		   "%.4s is not found in the\nMakeTermKey >> parameter files."
		   "\n", typeA, typeB, typeC);
	    exit(1);
	  }
	  mp->conf[i].amap.id[nterm].a = atmA;
	  mp->conf[i].amap.id[nterm].b = atmB;
	  mp->conf[i].amap.id[nterm].c = atmC;
	  nterm++;
	}
	else if (order == 4) {
	  for (m = 0; m < tp->HLC[j].HC[k].nt; m++) {

	    /*** For dihedrals, there is another layer of indexing. ***/
	    /*** We can identify the stiffness K, phase angle Phi,  ***/
	    /*** and periodicity N of the mth torsion term of this  ***/
	    /*** kth dihedral fourier series controlled by atom j   ***/
	    /*** of system i.                                       ***/ 
	    diheK = tp->HParam[tp->HLC[j].HC[k].t[m]].K;
	    dihePhi = tp->HParam[tp->HLC[j].HC[k].t[m]].Phi;
	    diheN = tp->HParam[tp->HLC[j].HC[k].t[m]].N;
	    singlet = (m == tp->HLC[j].HC[k].nt-1) ? -1.0 : 1.0;

	    /*** Now scan through the list of unique ***/
	    /*** torsion terms identified thus far.  ***/
	    found = 0;
	    maxcomp = 0;
	    for (n = 0; n < nmpterm; n++) {
	      compval = TypeCompare(typeA, typeB, typeC, typeD,
				    mp->torsions[n].atype,
				    mp->torsions[n].btype,
				    mp->torsions[n].ctype,
				    mp->torsions[n].dtype, 4, 0);
	      if (compval > maxcomp) {

		/*** This conditional is nested for clarity ***/
		if (fabs(mp->torsions[n].phase - dihePhi) < 1.0e-4 &&
		    fabs(mp->torsions[n].pn - diheN) < 1.0e-4) {

		  /*** Undo the previous assignment if one exists ***/
		  if (maxcomp == 1) {
		    mp->torsions[mp->conf[i].hmap.id[nterm].key].ninst -= 1;
		  }
		  maxcomp = compval;

		  /*** Assign the torsion to its place in the key ***/
		  mp->conf[i].hmap.id[nterm].key = n;
		  mp->torsions[n].ninst += 1;
		  found = 1;
		}
	      }
	    }
	    if (found == 0) {

	      /*** This torsion is not present in the parameter set ***/
	      printf("MakeTermKey >> Error.  A torsion term for %.4s %.4s "
		     "%.4s %.4s is not found\nMakeTermKey >> in the parameter "
		     "files.  Amplitude %9.4lf, periodicity\nMakeTermKey >> "
		     "%9.4lf, and phase angle %9.4lf were required.  Data "
		     "point\nMakeTermKey >> %s is not completely represented "
		     "by the parameter files.\n", typeA, typeB, typeC, typeD,
		     diheK, diheN, dihePhi, mp->conf[i].crdsrc);
	      resA = LocateResID(tp, atmA, 0, tp->nres);
	      resB = LocateResID(tp, atmB, 0, tp->nres);
	      resC = LocateResID(tp, atmC, 0, tp->nres);
	      resD = LocateResID(tp, atmD, 0, tp->nres);
	      printf("MakeTermKey >> Atoms in this dihedral:\n"
		     "MakeTermKey >>    %.4s %.4s\n"
		     "MakeTermKey >>    %.4s %.4s\n"
		     "MakeTermKey >>    %.4s %.4s\n"
		     "MakeTermKey >>    %.4s %.4s\n", &tp->AtomNames[4*resA],
		     &tp->AtomNames[4*atmA], &tp->AtomNames[4*resB],
		     &tp->AtomNames[4*atmB], &tp->AtomNames[4*resC],
		     &tp->AtomNames[4*atmC], &tp->AtomNames[4*resD],
		     &tp->AtomNames[4*atmD]);
	      exit(1);
	    }

	    /*** The number of torsion terms in this ***/
	    /*** system must be incremented here     ***/
	    mp->conf[i].hmap.id[nterm].a = atmA;
	    mp->conf[i].hmap.id[nterm].b = atmB;
	    mp->conf[i].hmap.id[nterm].c = atmC;
	    mp->conf[i].hmap.id[nterm].d = atmD;
	    nterm++;
	  }
	}
      }
    }

    /*** Store the number of bond, angle, or dihedral terms  ***/
    /*** in the conformation data.  This is NOT the number   ***/
    /*** of unique bond, angle, or torsion types--it is the  ***/
    /*** number of bonds, angles, or torsions IN a molecule. ***/
    if (order == 2) {
      mp->conf[i].bmap.nbond = nterm;
    }
    else if (order == 3) {
      mp->conf[i].amap.nangl = nterm;
    }
    else if (order == 4) {
      mp->conf[i].hmap.ntterm = nterm;
    }
  }

  /*** Store the number of unique terms in the parameter set ***/
  if (order == 2) {
    mp->nbond = nmpterm;
  }
  else if (order == 3) {
    mp->nangl = nmpterm;
  }
  else if (order == 4) {
    mp->ntor = nmpterm;
  }
}

/***=======================================================================***/
/*** FindAdjustableTerms: this routine parses the lists of bond, angle,    ***/
/***                      torsion, and scee / scnb terms to obtain a list  ***/
/***                      of adjustable parameters.  Adjustable parameters ***/
/***                      are suggested by the user in the form of atom    ***/
/***                      types or names of atoms making up each bonded    ***/
/***                      term or subsystem containing a scaling factor.   ***/
/***                      By matching these types or names to bonds in     ***/
/***                      the list of structures, the program identifies   ***/
/***                      the contributors to a linear least-squares fit   ***/
/***                      and, if necessary, iteratively refines that fit  ***/
/***                      by adjustment of nonlinear contributions such as ***/
/***                      bond length or dihedral phase angle.             ***/
/***                                                                       ***/
/*** Argument:                                                             ***/
/***   mp:     the set of systems and parameter guidelines                 ***/
/***=======================================================================***/
static void FindAdjustableTerms(prmset *mp)
{
  int i, j, found, ncol;

  /*** The fitting matrix column counter ***/
  ncol = 0;

  /*** All bonds, angles, and torsions default ***/
  /*** to no place in the fitting matrix       ***/
  mp->nbvar = 0;
  mp->navar = 0;
  mp->nhvar = 0;
  for (i = 0; i < mp->nbond; i++) {
    mp->bonds[i].fitcol = -1;
  }
  for (i = 0; i < mp->nangl; i++) {
    mp->angls[i].fitcol = -1;
  }
  for (i = 0; i < mp->ntor; i++) {
    mp->torsions[i].fitcol = -1;
  }

  /*** If all bonds are selected, enumerate them ***/
  if (mp->FitAllBonds == 1) {
    for (i = 0; i < mp->nbond; i++) {
      if (mp->bonds[i].ninst > 0) {
	mp->bonds[i].fitcol = ncol;
	mp->nbvar += 1;
	ncol++;
      }
    }
  }

  /*** Check adjustable bonds ***/
  else {
    for (i = 0; i < mp->nbadj; i++) {
      found = 0;
      for (j = 0; j < mp->nbond; j++) {
	if (mp->bonds[j].ninst == 0) {
	  continue;
	}
        if ((strncmp(mp->badj[i].atype, mp->bonds[j].atype, 4) == 0 &&
	     strncmp(mp->badj[i].btype, mp->bonds[j].btype, 4) == 0)||
	    (strncmp(mp->badj[i].btype, mp->bonds[j].atype, 4) == 0 &&
	     strncmp(mp->badj[i].atype, mp->bonds[j].btype, 4) == 0)) {
	  mp->bonds[j].fitcol = ncol;
	  mp->nbvar += 1;
	  found = 1;
	  ncol++;
        }
      }
      if (found == 0) {
        printf("mdgx >> Bond type %.4s %.4s marked for optimization but not "
	       "found in systems.\n", mp->badj[i].atype, mp->badj[i].btype);
      }
    }
  }

  /*** If all angles are selected, enumerate them ***/
  if (mp->FitAllAngles == 1) {
    for (i = 0; i < mp->nangl; i++) {
      if (mp->angls[i].ninst > 0) {
	mp->angls[i].fitcol = ncol;
	ncol++;
	mp->navar += 1;
      }
    }
  }

  /*** Check adjustable angles ***/
  else {
    for (i = 0; i < mp->naadj; i++) {
      found = 0;
      for (j = 0; j < mp->nangl; j++) {
	if (mp->angls[j].ninst == 0) {
	  continue;
	}
        if (strncmp(mp->aadj[i].btype, mp->angls[j].btype, 4) == 0 && 
	    ((strncmp(mp->aadj[i].atype, mp->angls[j].atype, 4) == 0 &&
	    strncmp(mp->aadj[i].ctype, mp->angls[j].ctype, 4) == 0)||
  	     (strncmp(mp->aadj[i].ctype, mp->angls[j].atype, 4) == 0 &&
	      strncmp(mp->aadj[i].atype, mp->angls[j].ctype, 4) == 0))) {
	  mp->angls[j].fitcol = ncol;
	  mp->navar += 1;
	  found = 1;
	  ncol++;
	}
      }
      if (found == 0) {
	printf("mdgx >> Angle type %.4s %.4s %.4s marked for optimization but "
	       "not found\nmdgx >> in systems.\n", mp->aadj[i].atype,
	       mp->aadj[i].btype, mp->aadj[i].ctype);
      }
    }
  }

  /*** If all torsions are selected, enumerate them ***/
  if (mp->FitAllTorsions == 1) {
    for (i = 0; i < mp->ntor; i++) {
      if (mp->torsions[i].ninst > 0 && mp->torsions[i].impr == 0) {
	mp->torsions[i].fitcol = ncol;
	mp->nhvar += 1;
	ncol++;
      }
    }
  }

  /*** Check adjustable torsions ***/
  else {
    for (i = 0; i < mp->nhadj; i++) {
      found = 0;
      for (j = 0; j < mp->ntor; j++) {
        if (TypeCompare(mp->torsions[j].atype, mp->torsions[j].btype,
			mp->torsions[j].ctype, mp->torsions[j].dtype,
			mp->hadj[i].atype, mp->hadj[i].btype,
			mp->hadj[i].ctype, mp->hadj[i].dtype,
			4, 1) == 2 && mp->torsions[j].ninst > 0) {
	  mp->torsions[j].fitcol = ncol;
	  mp->nhvar += 1;
	  found = 1;
	  ncol++;
        }
      }
      if (found == 0) {
        printf("mdgx >> Torsion type %.4s %.4s %.4s %.4s marked for "
	       "optimization but not found\nmdgx >> in systems.\n",
	       mp->hadj[i].atype, mp->hadj[i].btype, mp->hadj[i].ctype,
	       mp->hadj[i].dtype);
      }
    }
  }

  /*** Store the number of adjustable parameters ***/
  mp->nparm = ncol;

  /*** Check whether 1:4 scaling factors are in play ***/
  mp->nparm += mp->fitscnb + mp->fitscee;
}

/***=======================================================================***/
/*** BondTermParse:  compute the contributions due to a bonded term,       ***/
/***                 parsed into the kernel and coefficient-scaled values. ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   conf:   the molecular mechanics system, complete with topology,     ***/
/***           coordinates, and energy tables                              ***/
/***   idx:    the index of the bond in the conformation                   ***/
/***   order:  the order of the bonded term (2 = bond, 3 = angle...)       ***/
/***   mp:     the parameter set, containing master tables of terms        ***/
/***   krnl:   the kernel of the bond contribution (returned)              ***/
/***   cntrb:  the contribution, krnl*stiffness (returned)                 ***/
/***=======================================================================***/
static int BondTermParse(mmsys *conf, int idx, int order, prmset *mp,
			 double *val, double *krnl, double *cntrb)
{
  int i;
  double dx, dy, dz, r, dl, sangle, theta, costheta, mgba, mgbc, invbabc;
  double dtheta;
  double ab[3], ba[3], bc[3], ac[3], cd[3], crabbc[3], crbccd[3], scr[3];
  double *aptr, *bptr, *cptr, *dptr;
  xbonddef *bndtmp;
  xangldef *angtmp;
  torterm *tortmp;
  coord *crd;


  crd = &conf->crd;

  /*** Bond ***/
  if (order == 2) {

    /*** Compute displacement ***/
    aptr = &crd->loc[3*conf->bmap.id[idx].a];
    bptr = &crd->loc[3*conf->bmap.id[idx].b];
    dx = bptr[0] - aptr[0];
    dy = bptr[1] - aptr[1];
    dz = bptr[2] - aptr[2];
    r = sqrt(dx*dx + dy*dy + dz*dz);

    /*** Index into the bond table to get the ***/
    /*** stiffness and equilibrium length     ***/
    bndtmp = &mp->bonds[conf->bmap.id[idx].key];
    dl = bndtmp->l0 - r;
    *val = dl;
    *krnl = dl*dl;
    *cntrb = bndtmp->K*dl*dl;
  }

  /*** Angle ***/
  else if (order == 3) {

    /*** Compute displacements ***/
    aptr = &crd->loc[3*conf->amap.id[idx].a];
    bptr = &crd->loc[3*conf->amap.id[idx].b];
    cptr = &crd->loc[3*conf->amap.id[idx].c];
    for (i = 0; i < 3; i++) {
      ba[i] = aptr[i] - bptr[i];
      bc[i] = cptr[i] - bptr[i];
      ac[i] = cptr[i] - aptr[i];
    }
    mgba = ba[0]*ba[0] + ba[1]*ba[1] + ba[2]*ba[2];
    mgbc = bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2];
    invbabc = 1.0/sqrt(mgba*mgbc);
    costheta = (ba[0]*bc[0] + ba[1]*bc[1] + ba[2]*bc[2]) * invbabc;
    costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
    theta = acos(costheta);
    angtmp = &mp->angls[conf->amap.id[idx].key];
    dtheta = theta - angtmp->th0;
    *val = dtheta;
    *krnl = dtheta*dtheta;
    *cntrb = angtmp->K*dtheta*dtheta;
  }

  /*** Torsion ***/
  else if (order == 4) {

    /*** Compute displacements ***/
    aptr = &crd->loc[3*conf->hmap.id[idx].a];
    bptr = &crd->loc[3*conf->hmap.id[idx].b];
    cptr = &crd->loc[3*conf->hmap.id[idx].c];
    dptr = &crd->loc[3*conf->hmap.id[idx].d];
    for (i = 0; i < 3; i++) {
      ab[i] = bptr[i] - aptr[i];
      bc[i] = cptr[i] - bptr[i];
      cd[i] = dptr[i] - cptr[i];
    }
    CrossP(ab, bc, crabbc);
    CrossP(bc, cd, crbccd);
    costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
    costheta /=
      sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
           (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
    CrossP(crabbc, crbccd, scr);
    costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
    if (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) {
      theta = acos(costheta);
    }
    else {
      theta = -acos(costheta);
    }
    tortmp = &mp->torsions[conf->hmap.id[idx].key];
    sangle = tortmp->pn*theta - tortmp->phase;
    *val = theta;
    *krnl = 1.0 + cos(sangle);
    *cntrb = tortmp->K * (1.0 + cos(sangle));
  }

  /*** Alert if the contribution exceeds a specified tolerance ***/
  return (*cntrb > mp->mmtol) ? 1 : 0;
}

/***=======================================================================***/
/*** AllBondedTerms: compute the molecular mechanics contributions due to  ***/
/***                 bonded terms.                                         ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   conf:   the molecular mechanics system, complete with topology,     ***/
/***           coordinates, and energy tables                              ***/
/***   mp:     the parameter set, containing master tables of bonded terms ***/
/***   strain: indicator that there may be a bad interaction in this       ***/
/***           configuration, grounds to explore alternative geometries    ***/
/***=======================================================================***/
static double AllBondedTerms(mmsys *conf, prmset *mp, int *strain)
{
  int i, scon;
  double mmnrg;
  coord *crd;

  /*** Unpack the system; only the coordinates are needed, as all ***/
  /*** topology information is now stored in the prmset struct    ***/
  crd = &conf->crd;
  mmnrg = 0.0;

  /*** Flag for bad interactions, which may      ***/
 /*** indicate a misplaced atom in the molecule ***/
  scon = 0;

  /*** Bonds ***/
  for (i = 0; i < conf->bmap.nbond; i++) {
    scon += BondTermParse(conf, i, 2, mp, &conf->bmap.val[i],
			  &conf->bmap.kernel[i], &conf->bmap.contrib[i]);
    mmnrg += conf->bmap.contrib[i];
  }

  /*** Angles ***/
  for (i = 0; i < conf->amap.nangl; i++) {
    scon += BondTermParse(conf, i, 3, mp, &conf->amap.val[i],
			  &conf->amap.kernel[i], &conf->amap.contrib[i]);
    mmnrg += conf->amap.contrib[i];
  }

  /*** Dihedral fourier terms ***/
  for (i = 0; i < conf->hmap.ntterm; i++) {
    scon += BondTermParse(conf, i, 4, mp, &conf->hmap.val[i],
			  &conf->hmap.kernel[i], &conf->hmap.contrib[i]);
    mmnrg += conf->hmap.contrib[i];
  }

  /*** Are there any strained interactions? ***/
  *strain = scon;

  return mmnrg;
}

/***=======================================================================***/
/*** PlaceAtoms: this function references a list of atoms to rearrange a   ***/
/***             molecule.  A partial molecule is built, unless the list   ***/
/***             contains aliases for all atoms.                           ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   currvec:   the current state vector of the system... currvec[k]     ***/
/***              contains the identity of atmorder[k] in the real system  ***/
/***   atmorder:  the atom ordering vector... it is not simply 0, 1, 2,    ***/
/***              ... n because the molecule must be filled out one atom   ***/
/***              at a time, always with a means of checking whether the   ***/
/***              latest addition has introduced some impossible situation ***/
/***   asmts:     contains a more intelligible list of which atom is in    ***/
/***              what slot; asmts[k] is the number of the atom from the   ***/
/***              original coordinate which is now being used as atom k    ***/
/***              according to the topology.  This array is initialized to ***/
/***              -1 at the beginning of this routine, to carry additional ***/
/***              information about which atoms have NOT been placed.      ***/
/***   conf:      the molecular system, contains topology and coordinates  ***/
/***   ref:       the reference coordinates                                ***/
/***=======================================================================***/
static void PlaceAtoms(int* currvec, int* atmorder, int* asmts, int nassigned,
		       mmsys *conf, coord *ref)
{
  int i, j;

  SetIVec(asmts, conf->tp.natom, -1);
  for (i = 0; i < nassigned; i++) {
    for (j = 0; j < 3; j++) {
      conf->crd.loc[3*atmorder[i]+j] = ref->loc[3*currvec[i]+j];
    }
    asmts[atmorder[i]] = currvec[i];
  }
}

/***=======================================================================***/
/*** IncrementState: increment the state of the guess held in currvec.     ***/
/***                 Returns 1 if the highest possible state has been      ***/
/***                 reached.                                              ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   currvec:   the current state vector of the system... currvec[k]     ***/
/***              contains the identity of atmorder[k] in the real system  ***/
/***=======================================================================***/
static int IncrementState(int* currvec, int *nassigned, int natom)
{
  int i, pivot, inci;

  /*** Increment the ith counter of the state vector ***/
  /*** and then check backwards to verify that no    ***/
  /*** duplicate assignments have been made.  Lower- ***/
  /*** numbered counters take precedence if any      ***/
  /*** duplicates exist.                             ***/
  pivot = *nassigned-1;
  inci = 1;
  while (inci == 1 && pivot >= 0) {
    inci = 0;
    currvec[pivot] += 1;

    /*** If we have hit the roof, erase the most recently ***/
    /*** placed atoms and back-track until counters can   ***/
    /*** be incremented again.                            ***/
    if (currvec[pivot] == natom) {
      currvec[pivot] = -1;
      inci = 1;
      pivot--;
      continue;
    }

    /*** If we have successfully incremented, check to be ***/
    /*** sure that there are no duplicate atoms.          ***/
    for (i = 0; i < pivot; i++) {
      if (currvec[i] == currvec[pivot]) {
	inci = 1;
      }
    }
  }

  /*** The pivot may have shifted backwards ***/
  *nassigned = pivot+1;

  /*** Incrementing was successful and the new    ***/
  /*** state is a candidate for energy evaluation ***/      
  if (inci == 0) {
    return 0;
  }

  /*** The pivot had to be pushed all the way back to  ***/
  /*** -1, indicating that all states have been tested ***/
  else {
    return 1;
  }
}

/***=======================================================================***/
/*** WriteInpcrd: write an inpcrd file in AMBER format.  This routine is   ***/
/***              here because it is not currently needed anywhere else in ***/
/***              the program.                                             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   crd:      the coordinates                                           ***/
/***   bdim:     the box dimensions                                        ***/
/***   finame:   the name of the inpcrd file to write                      ***/
/***   tj:       the trajectory control data (for overwrite flag)          ***/
/***=======================================================================***/
static void WriteInpcrd(coord *crd, double* bdim, char* finame, trajcon *tj)
{
  int h, i;
  FILE *outp;

  outp = FOpenSafe(finame, tj->OverwriteOutput);
  fprintf(outp, "Coordinates generated by mdgx\n%6d\n", crd->natom);
  h = 0;
  for (i = 0; i < crd->natom; i++) {
    fprintf(outp, "%12.7lf%12.7lf%12.7lf", crd->loc[3*i], crd->loc[3*i+1],
	    crd->loc[3*i+2]);
    h++;
    if (h == 2) {
      fprintf(outp, "\n");
      h = 0;
    }
  }
  if (h != 0) {
    fprintf(outp, "\n");
  }
  for (i = 0; i < 6; i++) {
    fprintf(outp, "%12.7lf", bdim[i]);
  }
  fprintf(outp, "\n");
  fclose(outp);
}

/***=======================================================================***/
/*** RearrangeMol: rearrange the atoms of a molecule to minimize molecular ***/
/***               mechanics bonded term contributions.                    ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   conf:   the molecular mechanics system, complete with topology,     ***/
/***           coordinates, and energy tables                              ***/
/***   mp:     the parameter set, containing master tables of bonded terms ***/
/***   tj:     trajectory control data (for overwrite flag)                ***/
/***=======================================================================***/
static void RearrangeMol(mmsys *conf, prmset *mp, trajcon *tj)
{
  int i, j, atma, atmb, natom, nfound, nadd, nassigned, finished, strain;
  int iter, relaxed;
  int* atmorder;
  int* atmfound;
  int* currvec;
  int* bestvec;
  int* sampled;
  int* assignments;
  double dx, dy, dz, currnrg;
  double bdim[6];
  double* bestnrg;
  dmat bndmat, dstmat, stiffmat;
  coord locbuff;
  coord *crd;
  xbonddef *bndtmp;

  /*** Make a matrix of bonded terms ***/
  natom = conf->tp.natom;
  bndmat = CreateDmat(natom, natom, 0);
  stiffmat = CreateDmat(natom, natom, 0);
  for (i = 0; i < conf->bmap.nbond; i++) {
    atma = conf->bmap.id[i].a;
    atmb = conf->bmap.id[i].b;
    bndtmp = &mp->bonds[conf->bmap.id[i].key];
    bndmat.map[atma][atmb] = bndtmp->l0;
    bndmat.map[atmb][atma] = bndtmp->l0;
    stiffmat.map[atma][atmb] = bndtmp->K;
    stiffmat.map[atmb][atma] = bndtmp->K;
  }

  /*** The best energy after adding each atom ***/
  bestnrg = (double*)malloc(natom*sizeof(double));
  SetDVec(bestnrg, natom, 1.0e30);

  /*** Compute the inter-atomic distances in the coordinate set ***/
  crd = &conf->crd;
  dstmat = CreateDmat(natom, natom, 0);
  for (i = 0; i < natom; i++) {
    for (j = 0; j < natom; j++) {
      dx = crd->loc[3*j] - crd->loc[3*i];
      dy = crd->loc[3*j+1] - crd->loc[3*i+1];
      dz = crd->loc[3*j+2] - crd->loc[3*i+2];
      dstmat.map[i][j] = sqrt(dx*dx + dy*dy + dz*dz);
    }
  }

  /*** Plot a course through the molecule: in       ***/
  /*** what order shall we add atoms to the system? ***/
  atmorder = (int*)malloc(natom*sizeof(int));
  atmfound = (int*)calloc(natom, sizeof(int));
  atmorder[0] = 0;
  atmfound[0] = 1;
  nfound = 1;
  while (nfound < natom) {
    nadd = 0;
    for (i = 0; i < nfound; i++) {
      for (j = 0; j < conf->bmap.nbond; j++) {
	if (conf->bmap.id[j].a == atmorder[i] &&
	    atmfound[conf->bmap.id[j].b] == 0) {
	  atmfound[conf->bmap.id[j].b] = 1;
	  atmorder[nfound+nadd] = conf->bmap.id[j].b;
	  nadd++;
	}
	if (conf->bmap.id[j].b == atmorder[i] &&
	    atmfound[conf->bmap.id[j].a] == 0) {
	  atmfound[conf->bmap.id[j].a] = 1;
	  atmorder[nfound+nadd] = conf->bmap.id[j].a;
	  nadd++;
	}
      }
    }
    nfound += nadd;
  }

  /*** Copy the coordinates for safe keeping ***/
  locbuff = CopyCoord(&conf->crd);

  /*** Find putative identities for each atom based on bonding. ***/
  /*** The vectors currvec and bestvec record the current and   ***/
  /*** best configurations found thus far.  These vectors do    ***/
  /*** NOT hold the identity of kth atom in position k, rather  ***/
  /*** they hold the identity of the atmorder[k]th atom at      ***/
  /*** position k.  The molecule is filled out according to     ***/
  /*** the array atmorder.                                      ***/
  currvec = (int*)malloc(natom*sizeof(int));
  bestvec = (int*)malloc(natom*sizeof(int));
  sampled = (int*)calloc(natom, sizeof(int));
  assignments = (int*)calloc(natom, sizeof(int));
  SetIVec(assignments, natom, -1);
  SetIVec(currvec, natom, -1);
  nassigned = 2;
  currvec[0] = 0;
  currvec[1] = 1;
  finished = 0;
  iter = 0;
  relaxed = 0;
  while (finished == 0) {

    /*** Test this configuration to the extent  ***/
    /*** that atom positions have been assigned ***/
    PlaceAtoms(currvec, atmorder, assignments, nassigned, conf, &locbuff);
    currnrg = 0.0;
    strain = 0;
    for (i = 0; i < conf->bmap.nbond; i++) {
      if (assignments[conf->bmap.id[i].a] >= 0 &&
	  assignments[conf->bmap.id[i].b] >= 0) {
	strain += BondTermParse(conf, i, 2, mp, &dz, &dx, &dy);
	currnrg += dy;
      }
    }
    for (i = 0; i < conf->amap.nangl; i++) {
      if (assignments[conf->amap.id[i].a] >= 0 &&
	  assignments[conf->amap.id[i].b] >= 0 &&
	  assignments[conf->amap.id[i].c] >= 0) {
	strain += BondTermParse(conf, i, 3, mp, &dz, &dx, &dy);
	currnrg += dy;
      }
    }
    for (i = 0; i < conf->hmap.ntterm; i++) {
      if (assignments[conf->hmap.id[i].a] >= 0 &&
	  assignments[conf->hmap.id[i].b] >= 0 &&
	  assignments[conf->hmap.id[i].c] >= 0 &&
	  assignments[conf->hmap.id[i].d] >= 0) {
	strain += BondTermParse(conf, i, 4, mp, &dz, &dx, &dy);
	currnrg += dy;
      }
    }

    /*** Now, the energy has been summed; is it ***/
    /*** acceptable, and if so is it the best?  ***/
    if (strain == 0) {
      if (currnrg < bestnrg[nassigned-1]) {
	bestnrg[nassigned-1] = currnrg;
	if (nassigned == natom) {
	  relaxed = 1;
	  ReflectIVec(bestvec, currvec, natom);
	}
      }
      if (nassigned < natom) {
	currvec[nassigned] = -1;
	nassigned++;
      }
    }
    finished = IncrementState(currvec, &nassigned, natom);
    iter++;
  }
  if (relaxed == 0) {
    printf("mdgx >>   Optimization unsuccessful after %d iterations.\n", iter);
    conf->PassedEtol = 0;
  }
  else {
    printf("mdgx >>   Optimization complete after %d iterations.\n", iter);
    printf("mdgx >>   Best energy %12.6lf\n", bestnrg[natom-1]);
    conf->PassedEtol = 2;
    for (i = 0; i < natom; i++) {
      for (j = 0; j < 3; j++) {
	conf->crd.loc[3*atmorder[i]+j] = locbuff.loc[3*bestvec[i]+j];
      }
    }
    AllBondedTerms(conf, mp, &strain);
    for (i = 0; i < 3; i++) {
      bdim[i] = 100.0;
      bdim[i+3] = 90.0;
    }
    if (tj->OverwriteOutput == 1) {
      WriteInpcrd(&conf->crd, bdim, conf->crdsrc, tj);
    }
  }

  /*** Free allocated memory ***/
  free(bestnrg);
  free(currvec);
  free(bestvec);
  free(sampled);
  DestroyCoord(&locbuff);
}

/***=======================================================================***/
/*** CompBondedMME: compute molecular mechanics energy due to bonded       ***/
/***                interactions.                                          ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   conf:   the molecular mechanics system, complete with topology,     ***/
/***           coordinates, and energy tables                              ***/
/***   mp:     the parameter set, containing master tables of bonded terms ***/
/***   tj:     trajectory control data (for passing overwrite flag to      ***/
/***           RearrangeMol if needed)                                     ***/
/***=======================================================================***/
static void CompBondedMME(mmsys *conf, prmset *mp, trajcon *tj)
{
  int strain;
  double mmnrg;

  /*** Allocate tables ***/
  conf->bmap.val = (double*)malloc(conf->bmap.nbond*sizeof(double));
  conf->bmap.kernel = (double*)malloc(conf->bmap.nbond*sizeof(double));
  conf->bmap.contrib = (double*)malloc(conf->bmap.nbond*sizeof(double));
  conf->amap.val = (double*)malloc(conf->amap.nangl*sizeof(double));
  conf->amap.kernel = (double*)malloc(conf->amap.nangl*sizeof(double));
  conf->amap.contrib = (double*)malloc(conf->amap.nangl*sizeof(double));
  conf->hmap.val = (double*)malloc(conf->hmap.ntterm*sizeof(double));
  conf->hmap.kernel = (double*)malloc(conf->hmap.ntterm*sizeof(double));
  conf->hmap.contrib = (double*)malloc(conf->hmap.ntterm*sizeof(double));

  /*** Check for high molecular mechanics energy ***/
  mmnrg = AllBondedTerms(conf, mp, &strain);
  if (strain > 0) {
    printf("mdgx >> Rearrangement required for conformation:\nmdgx >> "
	   "  Topology %s / Coordinates %s.\n"
	   "mdgx >>   %d strained interactions, total bonded energy "
	   "%12.6lf\n", conf->tp.source, conf->crdsrc, strain, mmnrg);
    RearrangeMol(conf, mp, tj);
  }
  else {
    conf->PassedEtol = 1;
  }
}

/***=======================================================================***/
/*** CompNBMatrix: compute the matrix of non-bonded interactions for this  ***/
/***               molecule.  After this computation, the matrix will be   ***/
/***               parsed for exclusions and 1:4 interactions.             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:     the parameter set, containing fitting directives            ***/
/***   conf:   the molecular mechanics system, complete with topology,     ***/
/***           coordinates, and energy tables                              ***/
/***=======================================================================***/
static void CompNBMatrix(prmset *mp, mmsys *conf)
{
  int i, j, ni, nj, natom;
  double atmx, atmy, atmz, atmq, dx, dy, dz, invr, invr2, invr6, invr12;
  double Aval, Bval;
  double* ljAtmp;
  imat partnb;
  dmat excl, nrg;
  prmtop *tp;
  coord *crd;

  /*** Unpack the system ***/
  tp = &conf->tp;
  crd = &conf->crd;
  natom = tp->natom;

  /*** Allocate tables ***/
  nrg = CreateDmat(natom, natom, 0);
  excl = CreateDmat(natom, natom, 0);
  partnb = CreateImat(natom, natom);

  /*** Allocate the matrix.  Electrostatic interactions ***/
  /*** go above the diagonal, van-der Waals below it.   ***/
  for (i = 0; i < tp->natom-1; i++) {
    atmx = crd->loc[3*i];
    atmy = crd->loc[3*i+1];
    atmz = crd->loc[3*i+2];
    atmq = 332.0636*tp->Charges[i];
    if (tp->LJIdx[i] >= 0) {
      ljAtmp = tp->LJutab.map[tp->LJIdx[i]];
    }
    for (j = i+1; j < tp->natom; j++) {
      dx = crd->loc[3*j] - atmx;
      dy = crd->loc[3*j+1] - atmy;
      dz = crd->loc[3*j+2] - atmz;
      invr2 = 1.0/(dx*dx + dy*dy + dz*dz);
      invr = sqrt(invr2);
      invr6 = invr2*invr2*invr2;
      invr12 = invr6*invr6;

      /*** Compute electrostatic interaction ***/
      nrg.map[i][j] = atmq*tp->Charges[j]*invr;

      /*** Compute Lennard-Jones interaction ***/
      nrg.map[j][i] = 0.0;
      if (tp->LJIdx[i] >= 0 && tp->LJIdx[j] >= 0) {
	Aval = ljAtmp[2*tp->LJIdx[j]];
	Bval = ljAtmp[2*tp->LJIdx[j]+1];
	nrg.map[j][i] = Aval*invr12 + Bval*invr6;
      }
    }
  }

  /*** Exclusions as necessary ***/
  for (i = 0; i < tp->natom-1; i++) {
    for (j = i+1; j < tp->natom; j++) {
      excl.map[i][j] = 1.0 - TestBondedExclusion(i, j, tp);
      excl.map[j][i] = excl.map[i][j];
    }
  }

  /*** Mark 1:4 exclusions ***/
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->HLC[i].ndihe; j++) {
      if (tp->HLC[i].HC[j].eval14 == 1) {
	ni = tp->HLC[i].HC[j].a;
	nj = tp->HLC[i].HC[j].d;
	excl.map[ni][nj] = 1.0 - tp->HLC[i].HC[j].scee;
	excl.map[nj][ni] = 1.0 - tp->HLC[i].HC[j].scnb;
	partnb.map[ni][nj] = 1;
	partnb.map[nj][ni] = 1;
      }
    }
  }

  /*** Mark any other 1:1, 1:2, 1:3, and 1:4 eliminations ***/
  /*** (due to virtual sites, mainly)                     ***/
  if (tp->nxtrapt > 0) {
    for (i = 0; i < tp->natom; i++) {
      for (j = 0; j < tp->ElimPair[i].n11; j++) {
	ni = tp->ElimPair[i].list11[j].atmX;
	nj = tp->ElimPair[i].list11[j].atmY;
	excl.map[i][nj] = 0.0;
	excl.map[nj][i] = 0.0;
      }
      for (j = 0; j < tp->ElimPair[i].n12; j++) {
	ni = tp->ElimPair[i].list12[j].atmX;
	nj = tp->ElimPair[i].list12[j].atmY;
	excl.map[i][nj] = 0.0;
	excl.map[nj][i] = 0.0;
      }
      for (j = 0; j < tp->ElimPair[i].n13; j++) {
	ni = tp->ElimPair[i].list13[j].atmX;
	nj = tp->ElimPair[i].list13[j].atmY;
	excl.map[i][nj] = 0.0;
	excl.map[nj][i] = 0.0;
      }
      for (j = 0; j < tp->ElimPair[i].n14; j++) {
	ni = tp->ElimPair[i].list14[j].atmX;
	nj = tp->ElimPair[i].list14[j].atmY;
	excl.map[i][nj] = 1.0 - tp->HLC[i].HC[j].scee;
	excl.map[nj][i] = 1.0 - tp->HLC[i].HC[j].scnb;
      }
    }
  }

  /*** Compute the non-bonded electrostatic and Lennard-Jones kernels ***/
  conf->EEkernel = 0.0;
  conf->LJkernel = 0.0;
  conf->EEnonfit = 0.0;
  conf->LJnonfit = 0.0;
  for (i = 0; i < tp->natom-1; i++) {
    for (j = i+1; j < tp->natom; j++) {
      if (partnb.map[i][j] == 1) {
	if (mp->fitscee == 1) {
	  conf->EEkernel += nrg.map[i][j];
	}
	else {
	  conf->EEnonfit += nrg.map[i][j] * excl.map[i][j];
	}
	if (mp->fitscnb == 1) {
	  conf->LJkernel += nrg.map[j][i];
	}
	else {
	  conf->LJnonfit += nrg.map[j][i] * excl.map[j][i];
	}
      }
      else {
	conf->EEnonfit += nrg.map[i][j] * excl.map[i][j];
	conf->LJnonfit += nrg.map[j][i] * excl.map[j][i];
      }
    }
  }

  /*** Commit tables to the system ***/
  conf->nbnrg = nrg;
  conf->excl = excl;

  /*** Free allocated memory ***/
  DestroyImat(&partnb);
}

/***=======================================================================***/
/*** SumEnergyMM: sum the energy according to a molecular mechanics model. ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   conf:     the system of interest                                    ***/
/***=======================================================================***/
static double SumEnergyMM(mmsys *conf)
{
  int i, j, natom;
  double eQQ, eLJ, etot;

  eQQ = 0.0;
  eLJ = 0.0;
  natom = conf->tp.natom;
  for (i = 0; i < natom-1; i++) {
    for (j = i+1; j < natom; j++) {
      eQQ += conf->nbnrg.map[i][j] * conf->excl.map[i][j];
      eLJ += conf->nbnrg.map[j][i] * conf->excl.map[j][i];
    }
  }
  etot = DSum(conf->bmap.contrib, conf->bmap.nbond);
  etot += DSum(conf->amap.contrib, conf->amap.nangl);
  etot += DSum(conf->hmap.contrib, conf->hmap.ntterm);
  etot += eQQ + eLJ;

  return etot;
}

/***=======================================================================***/
/*** FitMatrixConstraints: add constraints to the fitting matrix of bonds, ***/
/***                       angles, and dihedral potentials.                ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   A:        the fitting matrix, in its original form before the QR    ***/
/***             decomposition                                             ***/
/***   b:        the target vector for the equation Ax = b                 ***/
/***   rstart:   the starting row of the fitting matrix, when we switch    ***/
/***             from data to restraints                                   ***/
/***   loo[key]: flag to indicate that leave-one-out cross-valdiation is   ***/
/***             in effect, and that column indices should therefore be    ***/
/***             taken from a specialized key, "lookey"                    ***/
/***=======================================================================***/
void FitMatrixConstraints(prmset *mp, dmat *A, double* b, int rstart, int loo,
			  int* lookey)
{
  int h, i, nrow, toridx, ninst, vcon;

  /*** Add constraints ***/
  nrow = rstart;
  if (mp->cnstB > 0.0) {
    nrow += mp->nbvar;
  }
  if (mp->cnstA > 0.0) {
    nrow += mp->navar;
  }
  if (mp->cnstH > 0.0) {
    vcon = 0;
    for (h = 0; h < mp->ntor; h++) {
      toridx = (loo == 0) ? mp->torsions[h].fitcol : lookey[h];
      if (toridx < 0) {
        continue;
      }
      ninst = 0;
      for (i = 0; i < rstart; i++) {
        if (fabs(A->map[i][toridx]) > 1.0e-8) {
          ninst += 1;
        }
      }
      if (ninst == 0) {
	ninst = 1;
      }
      if (mp->torsions[h].impr == 0) {
        A->map[vcon+nrow][toridx] = ninst*mp->cnstH;
        b[vcon+nrow] = 0.0;
      }
      else {
        A->map[vcon+nrow][toridx] = ninst*mp->cnstH;
        b[vcon+nrow] = ninst*mp->cnstH*mp->torsions[h].K;
      }
      vcon++;
    }
    nrow += mp->nhvar;
  }
}

/***=======================================================================***/
/*** FillFittingMatrix: fill up the fitting matrix given the set of        ***/
/***                    adjustable parameters, known evergy targets, and   ***/
/***                    molecular mechanics properties of each system.     ***/
/***                    Also loads energy targets into a pre-allocated     ***/
/***                    array.                                             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains history and number of      ***/
/***             fitting points)                                           ***/
/***   A:        the pre-allocated fitting matrix                          ***/
/***   b:        the pre-allocated array for computing Ax = B              ***/
/***=======================================================================***/
static void FillFittingMatrix(prmset *mp, dmat *A, double* b)
{
  int h, i, j, cidx, natom;
  double colsum, ctest;
  xbonddef *bondkey;
  xangldef *anglkey;
  torterm *torkey;
  dmat At;
  dmat *nbnrg, *excl;

  /*** Loop over all systems ***/
  for (h = 0; h < mp->nconf; h++) {

    /*** Start accumulating the MM energy components ***/
    /*** which are not subject to fitting            ***/
    mp->conf[h].nonfitmm = 0.0;

    /*** Bonds ***/
    for (i = 0; i < mp->conf[h].bmap.nbond; i++) {
      bondkey = &mp->bonds[mp->conf[h].bmap.id[i].key];
      cidx = bondkey->fitcol;
      if (cidx < 0) {
	mp->conf[h].nonfitmm += mp->conf[h].bmap.contrib[i];
      }
      else {
	A->map[h][cidx] += mp->conf[h].bmap.kernel[i];
      }
    }

    /*** Angles ***/
    for (i = 0; i < mp->conf[h].amap.nangl; i++) {
      anglkey = &mp->angls[mp->conf[h].amap.id[i].key];
      cidx = anglkey->fitcol;
      if (cidx < 0) {
	mp->conf[h].nonfitmm += mp->conf[h].amap.contrib[i];
      }
      else {
	A->map[h][cidx] += mp->conf[h].amap.kernel[i];
      }
    }

    /*** Torsions ***/
    for (i = 0; i < mp->conf[h].hmap.ntterm; i++) {
      torkey = &mp->torsions[mp->conf[h].hmap.id[i].key];
      cidx = torkey->fitcol;
      if (cidx < 0) {
	mp->conf[h].nonfitmm += mp->conf[h].hmap.contrib[i];
      }
      else {
	A->map[h][cidx] += mp->conf[h].hmap.kernel[i];
      }
    }

    /*** Non-bonded energy ***/
    natom = mp->conf[h].tp.natom;
    nbnrg = &mp->conf[h].nbnrg;
    excl = &mp->conf[h].excl;
    mp->conf[h].nonfitmm += mp->conf[h].EEnonfit;
    mp->conf[h].nonfitmm += mp->conf[h].LJnonfit;
    if (mp->fitscee == 1) {
      cidx = mp->nbvar + mp->navar + mp->nhvar;
      A->map[h][cidx] = mp->conf[h].EEkernel;
    }
    if (mp->fitscnb == 1) {
      cidx = mp->nbvar + mp->navar + mp->nhvar + 1;
      A->map[h][cidx] = mp->conf[h].LJkernel;
    }

    /*** Constant for energy offset between MM and QM models ***/
    A->map[h][mp->nparm + mp->conf[h].GroupNum] = 1.0;

    /*** Energy target ***/
    b[h] = mp->conf[h].enorm - mp->conf[h].nonfitmm;
  }

  /*** Add constraints ***/
  FitMatrixConstraints(mp, A, b, mp->nconf, 0, &i);

  /*** Check the matrix A for zero and highly correlated columns ***/
  mp->nzerocol = 0;
  mp->zerocol = (int*)malloc(A->col*sizeof(int));
  mp->ncorrcol = 0;
  mp->corrcol = (int*)malloc((A->col+1)*A->col*sizeof(int));
  mp->corrval = (double*)malloc(((A->col+1)*A->col/2)*sizeof(double));
  At = CreateDmat(mp->nparm, mp->nconf, 0);
  for (i = 0; i < mp->nparm; i++) {
    for (j = 0; j < mp->nconf; j++) {
      At.map[i][j] = A->map[j][i];
    }
  }
  for (i = 0; i < At.row; i++) {
    colsum = 0.0;
    for (j = 0; j < At.col; j++) {
      colsum += fabs(At.map[i][j]);
    }
    if (colsum < 1.0e-8) {
      mp->zerocol[mp->nzerocol] = i;
      mp->nzerocol += 1;
      for (j = 0; j < mp->nbond; j++) {
	if (mp->bonds[j].fitcol == i) {
	  printf("Zero column %4d: ( BOND %.4s %.4s , %4d instances )\n",
		 mp->bonds[j].fitcol, mp->bonds[j].atype, mp->bonds[j].btype,
                 mp->bonds[j].ninst);
	}
      }
      for (j = 0; j < mp->nangl; j++) {
	if (mp->angls[j].fitcol == i) {
	  printf("Zero column %4d: ( ANGL %.4s %.4s %.4s , %4d instances "
		 ")\n", mp->angls[j].fitcol, mp->angls[j].atype,
                 mp->angls[j].btype, mp->angls[j].ctype, mp->angls[j].ninst);
	}
      }
      for (j = 0; j < mp->ntor; j++) {
	if (mp->torsions[j].fitcol == i) {
	  printf("Zero column %4d: ( DIHE %.4s %.4s %.4s %.4s , %4d "
		 "instances )\n", mp->torsions[j].fitcol,
                 mp->torsions[j].atype, mp->torsions[j].btype,
		 mp->torsions[j].ctype, mp->torsions[j].dtype,
		 mp->torsions[j].ninst);
	  printf("             ( PN = %8.4lf Phase = %8.4lf )\n",
		 mp->torsions[j].pn, mp->torsions[j].phase);
	}
      }
    }
    for (j = i+1; j < At.row; j++) {
      ctest = Pearson(At.map[i], At.map[j], At.col);
      if (fabs(fabs(ctest) - 1.0) < 1.0e-4) {
	mp->corrcol[2*mp->ncorrcol] = i;
	mp->corrcol[2*mp->ncorrcol+1] = j;
	mp->corrval[mp->ncorrcol] = ctest;
	mp->ncorrcol += 1;
      }
    }
  }
}

/***=======================================================================***/
/*** FindColumnTerm: find the term, be it a bond, angle, or dihedral,      ***/
/***                 corresponding to the column of interest.  It is       ***/
/***                 terribly inefficient to search through all bonds,     ***/
/***                 angles, and dihedrals just to find the one that feeds ***/
/***                 into a particular matrix column when all of that data ***/
/***                 could be recorded, but this routine is only called at ***/
/***                 output and therefore saves some complexity in the     ***/
/***                 data structures.                                      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:      fitting control data (contains notes about the fitting     ***/
/***            matrix)                                                    ***/
/***   cc:      the column to match a term against                         ***/
/***   outp:    the output file                                            ***/
/***   order:   returned value, identifying the column as associated with  ***/
/***            a bond (return 2), angle (return 3), or torsion (return 4) ***/
/***   rval:    flag to set the actions taken upon finding the term        ***/
/***            corresponding to the column cc                             ***/
/***=======================================================================***/
static double FindColumnTerm(prmset *mp, int cc, FILE *outp, int *order,
			     int rval)
{
  int i;

  /*** Leading spaces for output in correlated columns output ***/
  if (rval == 0) {
    fprintf(outp, "  ");
  }
  for (i = 0; i < mp->nbond; i++) {
    if (mp->bonds[i].fitcol == cc) {
      *order = 2;
      if (rval <= 0) {
	fprintf(outp, " BOND %.2s %.2s      ", mp->bonds[i].atype,
		mp->bonds[i].btype);
      }
      else if (rval == 1) {
	return  mp->bonds[i].K;
      }
    }
  }
  for (i = 0; i < mp->nangl; i++) {
    if (mp->angls[i].fitcol == cc) {
      *order = 3;
      if (rval <= 0) {
	fprintf(outp, " ANGL %.2s %.2s %.2s   ", mp->angls[i].atype,
		mp->angls[i].btype, mp->angls[i].ctype);
      }
      else if (rval == 1) {
	return mp->angls[i].K;
      }
    }
  }
  for (i = 0; i < mp->ntor; i++) {
    if (mp->torsions[i].fitcol == cc) {
      if (rval <= 0) {
	*order = 4;
	if (mp->torsions[i].impr == 0) {
	  fprintf(outp, " DIHE ");
	}
	else {
	  fprintf(outp, " IMPR ");
	}
	fprintf(outp, "%.2s %.2s %.2s %.2s", mp->torsions[i].atype,
		mp->torsions[i].btype, mp->torsions[i].ctype,
		mp->torsions[i].dtype);
      }
      else if (rval == 1) {
	return mp->torsions[i].K;
      }
    }
  }

  /*** Carriage return for output relating to term glossary ***/
  if (rval == -1) {
    fprintf(outp, "\n");
  }

  /*** By default return zero ***/
  return 0.0;
}

/***=======================================================================***/
/*** CountRestraints: function for totaling up the number of restraints    ***/
/***                  on bonds, angles, and dihedrals.  General restraints ***/
/***                  on all these terms may be placed, as well as         ***/
/***                  specific restraints on individual fitted variables.  ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:  fitting control data                                           ***/
/***=======================================================================***/
static int CountRestraints(prmset *mp)
{
  int ncnst;

  /*** Initialize the counter ***/
  ncnst = 0;

  /*** Bonds ***/
  if (mp->cnstB > 0.0) ncnst += mp->nbvar;

  /*** Angles ***/
  if (mp->cnstA > 0.0) ncnst += mp->navar;

  /*** Dihedrals ***/
  if (mp->cnstH > 0.0) ncnst += mp->nhvar;

  /*** 1:4 fitting terms ***/
  ncnst += mp->fitscnb + mp->fitscee;

  return ncnst;
}

/***=======================================================================***/
/*** DestroyPrmset: detroy a prmset structure and all data in it; this may ***/
/***                not free data associated with conformations which have ***/
/***                been rejected from the fitting set, but those are      ***/
/***                freed as they are rejected.                            ***/
/***=======================================================================***/
static void DestroyPrmset(prmset *mp)
{
  int i;

  free(mp->zerocol);
  free(mp->corrcol);
  free(mp->corrval);
  free(mp->GroupCount);
  for (i = 0; i < mp->nconf; i++) {
    FreeConformation(&mp->conf[i]);
  }
  for (i = 0; i < mp->natom; i++) {
    free(mp->atoms[i].comment);
  }
  for (i = 0; i < mp->nbond; i++) {
    free(mp->bonds[i].comment);
  }
  for (i = 0; i < mp->nangl; i++) {
    free(mp->angls[i].comment);
  }
  for (i = 0; i < mp->ntor; i++) {
    free(mp->torsions[i].comment);
  }
  free(mp->conf);
  free(mp->badj);
  free(mp->aadj);
  free(mp->hadj);
  free(mp->atoms);
  free(mp->bonds);
  free(mp->angls);
  free(mp->torsions);
  free(mp->ititl);
  free(mp->icomm);
}

/***=======================================================================***/
/*** FitParams: the main function for optimizing parameters to fit a given ***/
/***            potential energy surface.                                  ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:  fitting control data, contains topology and coordinate file    ***/
/***        names, as well as keys for correlating parameters across       ***/
/***        systems.  Also contains the frcmod output file name.           ***/
/***   tj:  the trajectory control data, contains the statistics output    ***/
/***        file name                                                      ***/
/***=======================================================================***/
void FitParams(prmset *mp, trajcon *tj)
{
  int i, j, natom;
  double* b;
  double* bcpy;
  dmat A, Acpy;

  /*** Sort the command input information ***/
  CleanCommandInput(mp);

  /*** Read the conformations and topologies ***/
  for (i = 0; i < mp->nconf; i++) {
#ifdef MPI
    ReadSystemConf(&mp->conf[i], i, mp, tj);
#else
    ReadSystemConf(&mp->conf[i], i, mp);
#endif
  }

  /*** Correlate systems with similar topologies and compute enorm ***/
  GroupSystems(mp);

  /*** Read the underlying parameters ***/
  ReadParmFile(mp, tj);
  ReadFrcmodFile(mp, tj);

  /*** Correlate bond, angle, and torsion terms to each topology ***/
  MakeTermKey(mp, 2);
  MakeTermKey(mp, 3);
  MakeTermKey(mp, 4);

  /*** Correlate Lennard-Jones terms to each topology ***/

  /*** Identify terms that are adjustable by ***/
  /*** linear least-squares fitting, and tie ***/
  /*** them to fitting matrix columns        ***/
  FindAdjustableTerms(mp);

  /*** Compute molecular mechanics energies and    ***/
  /*** contributions from various adjustable terms ***/
  for (i = 0; i < mp->nconf; i++) {

    /*** Bonded interactions ***/
    CompBondedMME(&mp->conf[i], mp, tj);

    /*** Nonbonded interactions ***/
    natom = mp->conf[i].tp.natom;
    CompNBMatrix(mp, &mp->conf[i]);

    /*** Sum the original molecular mechanics energy ***/
    mp->conf[i].eorig = SumEnergyMM(&mp->conf[i]);
  }

  /*** Prune bad conformations ***/
  for (i = 0; i < mp->nconf; i++) {
    if (mp->conf[i].PassedEtol == 0) {
      FreeConformation(&mp->conf[i]);
      for (j = i+1; j < mp->nconf; j++) {
	mp->conf[j-1] = mp->conf[j];
      }
      mp->nconf -= 1;
      i--;
    }
  }

  /*** Adjust the target energies to meet the averages ***/
  /*** of the original MM energies for the surviving   ***/
  /*** conformations of each system                    ***/
  CompEnorm(mp);

  /*** Create the fitting matrix and solve ***/
  mp->ncnst = CountRestraints(mp);
  A = CreateDmat(mp->nconf + mp->ncnst, mp->nparm + mp->nunisys, 0);
  b = (double*)malloc((mp->nconf + mp->ncnst)*sizeof(double));
  FillFittingMatrix(mp, &A, b);
  CopyDmat(&Acpy, &A, 0);
  bcpy = CpyDVec(b, mp->nconf);
  AxbQRRxc(A, b, 0);
  BackSub(A, b);

  /*** Report the best parameters and the fit ***/
  ParamReport(mp, &Acpy, bcpy, b, tj);

  /*** Free allocated memory ***/
  free(b);
  free(bcpy);
  DestroyDmat(&A);
  DestroyDmat(&Acpy);
  DestroyPrmset(mp);

  /*** Exit!  We are done. ***/
  exit(1);
}
