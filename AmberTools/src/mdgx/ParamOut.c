#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
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
#include "ParamRead.h"
#include "ParamFit.h"

#include "CrdManipDS.h"

/***=======================================================================***/
/*** EvalParamFit: evaluate the degree to which a parameter set reproduces ***/
/***               the target energies.                                    ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains history and number of      ***/
/***             fitting points)                                           ***/
/***   nrg:      the energies produced by this parameter set               ***/
/***   etrg:     the target energies                                       ***/
/***   ewt:      the weights of each conformation in the fit               ***/
/***=======================================================================***/
static double EvalParamFit(prmset *mp, double* nrg, double* etrg,
                           double* ewt)
{
  int i;
  double dE, wrmsd;

  wrmsd = 0.0;
  for (i = 0; i < mp->nconf; i++) {
    dE = nrg[i] - etrg[i];
    wrmsd += dE*dE*ewt[i];
  }

  return sqrt(wrmsd/mp->nconf);
}

/***=======================================================================***/
/*** ExtractLJ: extract the Lennard-Jones sigma and epsilon values from a  ***/
/***            topology file for a certain atom type                      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tp:     the topology                                                ***/
/***   ljidx:  the Lennard-Jones type index                                ***/
/***   sig:    the sigma value (returned)                                  ***/
/***   eps:    the epsilon value (returned)                                ***/
/***=======================================================================***/
static void ExtractLJ(prmtop *tp, int ljidx, double *sig, double *eps)
{
  double aval, bval;

  if (ljidx >= 0) {
    aval = tp->LJA[tp->NBParmIdx[ljidx*tp->ntypes+ljidx]];
    bval = tp->LJB[tp->NBParmIdx[ljidx*tp->ntypes+ljidx]];
    *sig = pow(aval/bval, 1.0/6.0);
    *eps = 0.25*bval/pow(*sig, 6.0);
  }
  else {
    *sig = 0.0;
    *eps = 0.0;
  }
}

/***=======================================================================***/
/*** PrintParmAtoms: print the masses involved in a frcmod file for this   ***/
/***                   parameter fit.                                      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   outp:     the output frcmod file (pointer)                          ***/
/***=======================================================================***/
static void PrintParmAtoms(prmset *mp, FILE *outp)
{
  int i, j, k, atmid;
  prmtop *tp;

  /*** Determine which atoms appear in any of the systems in this fit ***/
  for (i = 0; i < mp->natom; i++) {
    mp->atoms[i].inreport = mp->reportall;
  }
  for (i = 0; i < mp->nunisys; i++) {
    for (j = 0; j < mp->nconf; j++) {
      if (mp->conf[j].GroupNum == i) {
	tp = &mp->conf[j].tp;
	for (k = 0; k < tp->natom; k++) {
	  atmid = CrossRefAtomType(mp, &tp->AtomTypes[4*k]);
	  mp->atoms[atmid].inreport = 1;
	}
      }
    }
  }

  for (i = 0; i < mp->natom; i++) {
    if (mp->atoms[i].inreport == 1) {
      fprintf(outp, "%.2s               %10.4lf  %10.4lf  %s\n",
	      mp->atoms[i].atype, mp->atoms[i].mass, mp->atoms[i].apol,
	      mp->atoms[i].comment);
    }
  }
  fprintf(outp, "\n");
}

/***=======================================================================***/
/*** TrimParmLine: trims a frcmod line to remove white space amidst atom   ***/
/***               type names.  All other characters on the line retain    ***/
/***               their original positions to preserve formatting.        ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   line:    the line to trim                                           ***/
/***   ntrim:   the number of characters to trim                           ***/
/***=======================================================================***/
static void TrimParmLine(char* line, int ntrim)
{
  int i, j, nlim;

  nlim = ntrim;
  for (i = 0; i < nlim; i++) {
    if (line[i] == ' ') {
      for (j = i; j < nlim; j++) {
        line[j] = line[j+1];
      }
      line[nlim-1] = ' ';
      i--;
      nlim--;
    }
  }
}

/***=======================================================================***/
/*** PrintParmBond: this procedure is much like PrintParmAtoms above.      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   order:    the order of bonded terms to seek out (2 = bonds,         ***/
/***             3 = angles, 4 = dihedrals)                                ***/
/***   outp:     the output file                                           ***/
/***   x:        the vector of fitted parameters                           ***/
/***=======================================================================***/
static void PrintParmBond(prmset *mp, int order, double* x, FILE *outp)
{
  int i, j, ilim;
  double stiffness;
  char tmpline[128];

  if (order == 2) {
    for (i = 0; i < mp->nbond; i++) {
      mp->bonds[i].inreport = mp->reportall;
    }
  }
  else if (order == 3) {
    for (i = 0; i < mp->nangl; i++) {
      mp->angls[i].inreport = mp->reportall;
    }
  }
  else if (order == 4 || order == 5) {
    for (i = 0; i < mp->ntor; i++) {
      mp->torsions[i].inreport = mp->reportall;
    }
  }
  else if (order == 9) {
    for (i = 0; i < mp->nhb1012; i++) {
      mp->hb1012[i].inreport = 1;
    }
  }
  ilim = (order == 2) ? mp->nbond : (order == 3) ? mp->nangl : (order == 9) ?
    mp->nhb1012 : mp->ntor;

  /*** Determine whether these bonds, angles, or dihedrals ***/
  /*** are present in ANY system used in this fit.         ***/
  for (i = 0; i < mp->nconf; i++) {
    if (order == 2) {
      for (j = 0; j < mp->conf[i].bmap.nbond; j++) {
	mp->bonds[mp->conf[i].bmap.id[j].key].inreport = 1;
      }
    }
    if (order == 3) {
      for (j = 0; j < mp->conf[i].amap.nangl; j++) {
	mp->angls[mp->conf[i].amap.id[j].key].inreport = 1;
      }
    }
    if (order == 4 || order == 5) {
      for (j = 0; j < mp->conf[i].hmap.ntterm; j++) {
	mp->torsions[mp->conf[i].hmap.id[j].key].inreport = 1;
      }
    }
  }

  /*** Print bonds, angles, and dihedrals to the new parameter file ***/
  for (i = 0; i < ilim; i++) {
    if (order == 2 && mp->bonds[i].inreport == 1) {
      if (mp->bonds[i].fitcol >= 0) {
	stiffness = x[mp->bonds[i].fitcol];
	sprintf(mp->bonds[i].comment, mp->icomm);
      }
      else {
	stiffness = mp->bonds[i].K;
      }
      sprintf(tmpline, "%.2s-%.2s            %10.4lf  %10.4lf  %s\n",
	      mp->bonds[i].atype, mp->bonds[i].btype, stiffness,
	      mp->bonds[i].l0, mp->bonds[i].comment);
      fprintf(outp, "%s", tmpline);
    }
    if (order == 3 && mp->angls[i].inreport == 1) {
      if (mp->angls[i].fitcol >= 0) {
	stiffness = x[mp->angls[i].fitcol];
	sprintf(mp->angls[i].comment, mp->icomm);
      }
      else {
	stiffness = mp->angls[i].K;
      }
      sprintf(tmpline, "%.2s-%.2s-%.2s         %10.4lf  %10.2lf  %s\n",
              mp->angls[i].atype, mp->angls[i].btype, mp->angls[i].ctype,
              stiffness, (180.0/PI)*mp->angls[i].th0, mp->angls[i].comment);
      fprintf(outp, "%s", tmpline);
    }
    if (order == 4 && mp->torsions[i].impr == 0 &&
	mp->torsions[i].inreport == 1) {
      if (mp->torsions[i].fitcol >= 0) {
	stiffness = x[mp->torsions[i].fitcol];
	sprintf(mp->torsions[i].comment, mp->icomm);
      }
      else {
	stiffness = mp->torsions[i].K;
      }
      sprintf(tmpline, "%.2s-%.2s-%.2s-%.2s   1  %10.5lf %6.1lf %4.1lf  "
	      "%s\n", mp->torsions[i].atype, mp->torsions[i].btype,
              mp->torsions[i].ctype, mp->torsions[i].dtype, stiffness,
	      mp->torsions[i].phase*180.0/PI,
	      mp->torsions[i].singlet * mp->torsions[i].pn,
	      mp->torsions[i].comment);
      fprintf(outp, "%s", tmpline);
    }
    if (order == 5 && mp->torsions[i].impr == 1 &&
	mp->torsions[i].inreport == 1) {
      if (mp->torsions[i].fitcol >= 0) {
	stiffness = x[mp->torsions[i].fitcol];
	sprintf(mp->torsions[i].comment, mp->icomm);
      }
      else {
	stiffness = mp->torsions[i].K;
      }
      sprintf(tmpline, "%.2s-%.2s-%.2s-%.2s      %10.5lf %6.1lf %4.1lf  "
	      "%s\n", mp->torsions[i].atype, mp->torsions[i].btype,
              mp->torsions[i].ctype, mp->torsions[i].dtype, stiffness,
	      mp->torsions[i].phase*180.0/PI,
              mp->torsions[i].singlet * mp->torsions[i].pn,
	      mp->torsions[i].comment);
      fprintf(outp, "%s", tmpline);
    }
    if (order == 9 && mp->hb1012[i].inreport == 1) {

      /*** Currently there is no support for ***/
      /*** fitting H-bond 10-12 potentials.  ***/
      sprintf(tmpline, "  %.2s  %.2s         %10.4lf  %10.4lf  %s\n",
	      mp->hb1012[i].atype, mp->hb1012[i].btype, mp->hb1012[i].Aterm,
	      mp->hb1012[i].Bterm, mp->hb1012[i].comment);
      fprintf(outp, "%s", tmpline);
    }
  }
  fprintf(outp, "\n");
}

/***=======================================================================***/
/*** AccumulateSamplingTable: accumulate the sampling histogram based on   ***/
/***                          all conformations.  This routine is          ***/
/***                          encapsulated to allow it to be used in       ***/
/***                          multiple cases.                              ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   S:        the sampling matrix                                       ***/
/***   order:    the order of bonded terms to seek out (2 = bonds,         ***/
/***             3 = angles, 4 = dihedrals)                                ***/
/***   sysid:    system ID number, -1 for accumulation over all systems    ***/
/***   track:    activates tracking of the exact values of the variable;   ***/
/***             only useable if there is only 1 row in S                  ***/
/***   record:   all values of the variable encountered                    ***/
/***=======================================================================***/
static void AccumulateSamplingTable(prmset *mp, imat *S, int order, int sysid,
				    int track, double* record)
{
  int i, j, dev, nrec;
  double bdisp, adisp, hdisp;
  bidx *tbond;
  aidx *tangl;
  hidx *tdihe;

  /*** Compile sampling histograms ***/
  nrec = 0;
  for (i = 0; i < mp->nconf; i++) {
    if (sysid >= 0 && mp->conf[i].GroupNum != sysid) {
      continue;
    }
    if (order == 2) {
      for (j = 0; j < mp->conf[i].bmap.nbond; j++) {
        tbond = &mp->conf[i].bmap.id[j];
        if (mp->bonds[tbond->key].fitcol < 0) {
          continue;
        }
        bdisp = mp->conf[i].bmap.val[j];
        dev = (bdisp / 0.03) + 9.00001;
        if (dev >= 0 && dev < 18) {
          S->map[mp->bonds[tbond->key].fitcol][dev] += 1;
        }
	if (track == 1) {
	  record[nrec] = bdisp;
	  nrec++;
	}
      }
    }
    else if (order == 3) {
      for (j = 0; j < mp->conf[i].amap.nangl; j++) {
        tangl = &mp->conf[i].amap.id[j];
        if (mp->angls[tangl->key].fitcol < 0) {
          continue;
        }
        adisp = mp->conf[i].amap.val[j];
        dev = (adisp / 0.05) + 9.00001;
        if (dev >= 0 && dev < 18) {
          S->map[mp->angls[tangl->key].fitcol][dev] += 1;
        }
	if (track == 1) {
	  record[nrec] = adisp;
	  nrec++;
	}
      }
    }
    else {
      for (j = 0; j < mp->conf[i].hmap.ntterm; j++) {
        tdihe = &mp->conf[i].hmap.id[j];
        if (mp->torsions[tdihe->key].fitcol < 0) {
          continue;
        }
        hdisp = mp->conf[i].hmap.val[j];
        dev = hdisp * 9.0 / PI + 9.00001;
        if (dev >= 0 && dev < 18) {
          S->map[mp->torsions[tdihe->key].fitcol][dev] += 1;
        }
	if (track == 1) {
	  record[nrec] = hdisp;
	  nrec++;
	}
      }
    }
  }
}

/***=======================================================================***/
/*** CharacterHistogram: print a character-keyed histogram, values ranging ***/
/***                     from 0 to 10.                                     ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   S:      integer matrix containing the histogram data                ***/
/***   ncol:   the number of bins in the histogram                         ***/
/***   outp:   the output file                                             ***/
/***=======================================================================***/
static void CharacterHistogram(int* S, int ncol, FILE *outp)
{
  int i;
  char occmap[16];

  sprintf(occmap, " .:=eoUO0@");
  for (i = 0; i < ncol; i++) {
    if (S[i] < 10) {
      fprintf(outp, "%c", occmap[S[i]]);
    }
    else {
      fprintf(outp, "X");
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
/*** SamplingMatrixOutput: print the output of the sampling matrix, in     ***/
/***                       human-readable format.                          ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   order:    the order of bonded terms to seek out (2 = bonds,         ***/
/***             3 = angles, 4 = dihedrals)                                ***/
/***   outp:     the output file                                           ***/
/***   S:        the sampling matrix                                       ***/
/***   x:        the solution vector                                       ***/
/***=======================================================================***/
static void SamplingMatrixOutput(prmset *mp, int order, FILE *outp, imat *S,
                                 double* x)
{
  int i, j, nvar, ninst;
  double Korig;

  fprintf(outp, "  Term  Amber Types  Count -PI"
          "     +0      PI   Fitted  Original\n");
  fprintf(outp, "  ----  -----------  ----- "
          "------------------  -------- --------\n");
  nvar = 0;
  for (i = 0; i < S->row; i++) {
    if (order == 2) {
      while (mp->bonds[nvar].fitcol < 0) {
        nvar++;
      }
      ninst = mp->bonds[nvar].ninst;
      fprintf(outp, "  BOND  %.2s %.2s       ", mp->bonds[nvar].atype,
              mp->bonds[nvar].btype);
    }
    if (order == 3) {
      while (mp->angls[nvar].fitcol < 0) {
        nvar++;
      }
      ninst = mp->angls[nvar].ninst;
      fprintf(outp, "  ANGL  %.2s %.2s %.2s    ", mp->angls[nvar].atype,
              mp->angls[nvar].btype, mp->angls[nvar].ctype);
    }
    if (order == 4) {
      while (mp->torsions[nvar].fitcol < 0) {
        nvar++;
      }
      ninst = mp->torsions[nvar].ninst;
      if (mp->torsions[nvar].impr == 0) {
        fprintf(outp, "  DIHE  ");
      }
      else {
        fprintf(outp, "  IMPR  ");
      }
      fprintf(outp, "%.2s %.2s %.2s %.2s ", mp->torsions[nvar].atype,
              mp->torsions[nvar].btype, mp->torsions[nvar].ctype,
              mp->torsions[nvar].dtype);
    }
    nvar++;
    fprintf(outp, "%6d ", ninst);
    CharacterHistogram(S->map[i], S->col, outp);
    Korig = FindColumnTerm(mp, i, outp, &j, 1);
    fprintf(outp, "  %8.2lf %8.2lf\n", x[i], Korig);
  }
}

/***=======================================================================***/
/*** ParameterWarnings: runs a series of diagnostics on the parameter and  ***/
/***                    its possible effects on the resulting model.       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   n:        the column of the parameter to detail                     ***/
/***   sysid:    the system ID for context                                 ***/
/***   outp:     the output file                                           ***/
/***=======================================================================***/
static void ParameterWarnings(prmset *mp, int n, int sysid, FILE *outp)
{
  int order;
  double Korig;

  /*** Get the original stiffness and order of this parameter ***/
  Korig = FindColumnTerm(mp, n, outp, &order, 1);

  /*** Make a complete list of every time this parameter is ***/
  /*** encountered in the context of a particular system    ***/

  /*** Carriage return to end the line in the output file ***/
  fprintf(outp, "\n");
}

/***=======================================================================***/
/*** DetailParameter: details the instances in which a parameter occurs in ***/
/***                  the fitting set.                                     ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   n:        the column of the parameter to detail                     ***/
/***   outp:     the output file                                           ***/
/***=======================================================================***/
static void DetailParameter(prmset *mp, int n, FILE *outp)
{
  int i, j, bkey, akey, hkey, ninst, atma, atmb, atmc, atmd;
  int resa, resb, resc, resd, nspc, rescon, bcol, acol, hcol;
  int* systems;
  double* dtmp;
  char *tpatoms, *tpres;
  imat* BondSamp;
  imat* AnglSamp;
  imat* DiheSamp;
  itrack* inst;

  /*** Determine the maximum possible number of instances ***/
  j = 0;
  for (i = 0; i < mp->nconf; i++) {
    j += mp->conf[i].bmap.nbond;
    j += mp->conf[i].amap.nangl;
    j += mp->conf[i].hmap.ntterm;
  }
  inst = (itrack*)malloc(j*sizeof(itrack));

  /*** Check for this term in exactly one example of each system ***/
  ninst = 0;
  systems = (int*)calloc(mp->nunisys, sizeof(int));
  for (i = 0; i < mp->nconf; i++) {
    if (mp->conf[i].GroupNum < 0 || systems[mp->conf[i].GroupNum] == 1) {
      continue;
    }
    tpatoms = mp->conf[i].tp.AtomNames;
    tpres = mp->conf[i].tp.ResNames;
    systems[mp->conf[i].GroupNum] = 1;
    for (j = 0; j < mp->conf[i].bmap.nbond; j++) {
      bkey = mp->conf[i].bmap.id[j].key;
      if (mp->bonds[bkey].fitcol == n) {

        /*** Column n pertains to a bond in this molecule ***/
        inst[ninst].sysid = mp->conf[i].GroupNum;
        inst[ninst].sysno = i;
        inst[ninst].order = 2;
        inst[ninst].tnum = j;
        atma = mp->conf[i].hmap.id[j].a;
        atmb = mp->conf[i].hmap.id[j].b;
        sprintf(inst[ninst].atom, "%.4s%.4s", &tpatoms[4*atma],
                &tpatoms[4*atmb]);
        sprintf(inst[ninst].res, "%.4s%.4s", &tpres[4*atma], &tpres[4*atmb]);
        ninst++;
      }
    }
    for (j = 0; j < mp->conf[i].amap.nangl; j++) {
      akey = mp->conf[i].amap.id[j].key;
      if (mp->angls[akey].fitcol == n) {

        /*** Column n pertains to an angle in this molecule ***/
        inst[ninst].sysid = mp->conf[i].GroupNum;
        inst[ninst].sysno = i;
        inst[ninst].order = 3;
        inst[ninst].tnum = j;
        atma = mp->conf[i].hmap.id[j].a;
        atmb = mp->conf[i].hmap.id[j].b;
        atmc = mp->conf[i].hmap.id[j].c;
        sprintf(inst[ninst].atom, "%.4s%.4s%.4s", &tpatoms[4*atma],
                &tpatoms[4*atmb], &tpatoms[4*atmc]);
        sprintf(inst[ninst].res, "%.4s%.4s%.4s", &tpres[4*atma],
                &tpres[4*atmb], &tpres[4*atmc]);
        ninst++;
      }
    }
    for (j = 0; j < mp->conf[i].hmap.ntterm; j++) {
      hkey = mp->conf[i].hmap.id[j].key;
      if (mp->torsions[hkey].fitcol == n) {

        /*** Column n pertains to a dihedral in this molecule ***/
        inst[ninst].sysid = mp->conf[i].GroupNum;
        inst[ninst].sysno = i;
        inst[ninst].order = 4;
        inst[ninst].tnum = j;
        atma = mp->conf[i].hmap.id[j].a;
        atmb = mp->conf[i].hmap.id[j].b;
        atmc = mp->conf[i].hmap.id[j].c;
        atmd = mp->conf[i].hmap.id[j].d;
        resa = LocateResID(&mp->conf[i].tp, atma, 0, mp->conf[i].tp.nres);
        resb = LocateResID(&mp->conf[i].tp, atmb, 0, mp->conf[i].tp.nres);
        resc = LocateResID(&mp->conf[i].tp, atmc, 0, mp->conf[i].tp.nres);
        resd = LocateResID(&mp->conf[i].tp, atmd, 0, mp->conf[i].tp.nres);
        sprintf(inst[ninst].atom, "%.4s%.4s%.4s%.4s", &tpatoms[4*atma],
                &tpatoms[4*atmb], &tpatoms[4*atmc], &tpatoms[4*atmd]);
        sprintf(inst[ninst].res, "%.4s%.4s%.4s%.4s", &tpres[4*resa],
                &tpres[4*resb], &tpres[4*resc], &tpres[4*resd]);
        ninst++;
      }
    }
  }

  /*** Create sampling tables for conformations  ***/
  /*** sampled by each individual system.  dtmp  ***/
  /*** is set just to avoid compiler complaints. ***/
  if (mp->nbvar > 0) BondSamp = (imat*)malloc(mp->nunisys*sizeof(imat));
  if (mp->navar > 0) AnglSamp = (imat*)malloc(mp->nunisys*sizeof(imat));
  if (mp->nhvar > 0) DiheSamp = (imat*)malloc(mp->nunisys*sizeof(imat));
  dtmp = mp->corrval;
  for (i = 0; i < mp->nunisys; i++) {
    if (mp->nbvar > 0) {
      BondSamp[i] = CreateImat(mp->nbvar, 18);
      AccumulateSamplingTable(mp, &BondSamp[i], 2, i, 0, dtmp);
    }
    if (mp->navar > 0) {
      AnglSamp[i] = CreateImat(mp->navar, 18);
      AccumulateSamplingTable(mp, &AnglSamp[i], 3, i, 0, dtmp);
    }
    if (mp->nhvar > 0) {
      DiheSamp[i] = CreateImat(mp->nhvar, 18);
      AccumulateSamplingTable(mp, &DiheSamp[i], 4, i, 0, dtmp);
    }
  }

  /*** Recount the instances ***/
  for (i = 0; i < ninst; i++) {
    fprintf(outp, "  ");
    fprintf(outp, "%.4s(", inst[i].res);
    rescon = 0;
    nspc = 7;
    for (j = 0; j < inst[i].order; j++) {
      if (strncmp(&inst[i].res[4*j], &inst[i].res[4*rescon], 4) == 0) {
        fprintf(outp, "%.4s", &inst[i].atom[4*j]);
        nspc += 4;
      }
      else {
        fprintf(outp, ") %.4s(%.4s", &inst[i].res[4*j], &inst[i].atom[4*j]);
        rescon = j;
        nspc += 11;
      }
      if (j < inst[i].order-1 &&
          strncmp(&inst[i].res[4*(j+1)], &inst[i].res[4*rescon], 4) == 0) {
        fprintf(outp, " ");
        nspc += 1;
      }
    }
    fprintf(outp, ")");
    nspc += 1;
    for (j = nspc; j < 41; j++) {
      fprintf(outp, " ");
    }
    fprintf(outp, "%5d  ", mp->GroupCount[inst[i].sysid]);
    if (inst[i].order == 2) {
      bkey = mp->conf[inst[i].sysno].bmap.id[inst[i].tnum].key;
      bcol = mp->bonds[bkey].fitcol;
      CharacterHistogram(BondSamp[inst[i].sysid].map[bcol], 18, outp);
    }
    else if (inst[i].order == 3) {
      akey = mp->conf[inst[i].sysno].amap.id[inst[i].tnum].key;
      acol = mp->angls[akey].fitcol;
      CharacterHistogram(AnglSamp[inst[i].sysid].map[acol], 18, outp);
    }
    else if (inst[i].order == 4) {
      hkey = mp->conf[inst[i].sysno].hmap.id[inst[i].tnum].key;
      hcol = mp->torsions[hkey].fitcol;
      CharacterHistogram(DiheSamp[inst[i].sysid].map[hcol], 18, outp);
    }
    ParameterWarnings(mp, n, inst[i].sysid, outp);
  }

  /*** Free allocated memory ***/
  for (i = 0; i < mp->nunisys; i++) {
    if (mp->nbvar > 0) DestroyImat(&BondSamp[i]);
    if (mp->navar > 0) DestroyImat(&AnglSamp[i]);
    if (mp->nhvar > 0) DestroyImat(&DiheSamp[i]);
  }
  if (mp->nbvar > 0) free(BondSamp);
  if (mp->navar > 0) free(AnglSamp);
  if (mp->nhvar > 0) free(DiheSamp);
  free(systems);
  free(inst);
}

/***=======================================================================***/
/*** PrintSamplingTable: print a table to describe the sampling of each    ***/
/***                     fitted parameter across all conformations.        ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   order:    the order of bonded terms to seek out (2 = bonds,         ***/
/***             3 = angles, 4 = dihedrals)                                ***/
/***   x:        the solution vector (to print alongside sampling data)    ***/
/***   outp:     the output file                                           ***/
/***=======================================================================***/
static void PrintSamplingTable(prmset *mp, int order, double* x, FILE *outp)
{
  double* dtmp;
  imat samp;

  /*** Allocate the histogram matrix ***/
  samp = (order == 2) ? CreateImat(mp->nbvar, 18) :
    (order == 3) ? CreateImat(mp->navar, 18) : CreateImat(mp->nhvar, 18);

  /*** Accumulate the histogram ***/
  AccumulateSamplingTable(mp, &samp, order, -1, 0, dtmp);

  /*** Print the histogram ***/
  if (order == 2 && mp->nbvar > 0) {
    fprintf(outp, " Bond sampling:\n");
    SamplingMatrixOutput(mp, 2, outp, &samp, x);
  }
  if (order == 3 && mp->navar > 0) {
    fprintf(outp, " Angle sampling:\n");
    SamplingMatrixOutput(mp, 3, outp, &samp, x);
  }
  if (order == 4 && mp->nhvar > 0) {
    fprintf(outp, " Torsion term sampling:\n"
            "                             Sample, 18 bins \n");
    SamplingMatrixOutput(mp, 4, outp, &samp, x);
  }

  /*** Free allocated memory ***/
  DestroyImat(&samp);
}

/***=======================================================================***/
/*** PrintMatrixAnalysis: print an analysis of the fitting matrix,         ***/
/***                      identifying highly correlated columns as well as ***/
/***                      columns with no fitting data in them.            ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:      fitting control data (contains notes about the fitting     ***/
/***            matrix)                                                    ***/
/***   x:       the solution vector                                        ***/
/***   outp:    the output file                                            ***/
/***=======================================================================***/
static void PrintMatrixAnalysis(prmset *mp, double* x, FILE *outp)
{
  int i, order;

  /*** Zero-data columns? ***/
  if (mp->nzerocol == 0) {
    fprintf(outp, " - No columns with zero data detected\n");
  }
  else {
    fprintf(outp, " - %d columns with zero data detected, corresponding to "
            "parameters:\n", mp->nzerocol);
    for (i = 0; i < mp->nzerocol; i++) {
      FindColumnTerm(mp, mp->zerocol[i], outp, &order, 0);
      fprintf(outp, "\n");
    }
  }

  /*** Correlated columns? ***/
  if (mp->ncorrcol == 0) {
    fprintf(outp, " - No correlated columns detected\n");
  }
  else {
    fprintf(outp, " - %d highly correlated column pairs detected, "
            "corresponding to parameters:\n\n", mp->ncorrcol);
    fprintf(outp, "     First term,        Second term,     Fitted   Fitted   "
            "   Pearson\n");
    fprintf(outp, "   Type    Atoms      Type    Atoms      value 1  value 2  "
            " Correlation\n");
    fprintf(outp, "   ---- -----------   ---- -----------   -------- -------- "
            " -----------\n");
    for (i = 0; i < mp->ncorrcol; i++) {
      FindColumnTerm(mp, mp->corrcol[2*i], outp, &order, 0);
      FindColumnTerm(mp, mp->corrcol[2*i+1], outp, &order, 0);
      fprintf(outp, "   %8.2lf %8.2lf  %11.8lf\n",
              x[mp->corrcol[2*i]], x[mp->corrcol[2*i+1]], mp->corrval[i]);
    }
  }
}

/***=======================================================================***/
/*** ParameterDescriptions: this routine will print information to put all ***/
/***                        fitted parameters in context.  It details the  ***/
/***                        residues and atoms that each parameter affects ***/
/***                        and gives counts on the number of instances of ***/
/***                        each occurrence.  This routine will also check ***/
/***                        for criteria that may indicate a parameter has ***/
/***                        been poorly fitted.  Output is tabulated in a  ***/
/***                        special section of the mdout file.             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   outp:     pointer to the mdout file                                 ***/
/***=======================================================================***/
static void ParameterDescription(prmset *mp, FILE *outp)
{
  int i, order;

  PrintVADesc(0, "(4.)", 4, " ", 1, "Context of each parameter.  Instances of "
              "each parameter found in the fitting conformations are "
              "enumerated below.\n", 74, 0, outp);
  fprintf(outp, " [ Parameter ]                                  Sampling "
          "Histogram\n"
          " Residue (Atom Names)                    Count  -PI     +0      PI "
          " Warnings\n"
          " --------------------------------------- -----  ------------------ "
          " --------\n");
  for (i = 0; i < mp->nparm; i++) {
    FindColumnTerm(mp, i, outp, &order, -1);
    DetailParameter(mp, i, outp);
  }
}

/***=======================================================================***/
/*** CheckMinima: this function computes the minima of each fourier series ***/
/***              and determines whether any of them are significantly     ***/
/***              lower than the points at which the fourier series is     ***/
/***              sampled in the training data.                            ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   tj:       trajectory control information (for the input file name)  ***/
/***   mdout:    pointer to the mdout file                                 ***/
/***   x:        the solution vector from the fit                          ***/
/***=======================================================================***/
static void CheckMinima(prmset *mp, trajcon *tj, double* x)
{
  int i, j, ntt, nser, totsamp;
  int* examined;
  int* sterm;
  double hangle, sangle, tfunc;
  FILE *outp;
  time_t ct;
  torterm *tortmp;

  /*** Checks are done here to avoid cluttering ParamReport() ***/
  if (mp->sr[0] == '\0' || mp->nhvar == 0) {
    return;
  }

  /*** Array to flag torsions that have ***/
  /*** been included in fourier series  ***/
  examined = (int*)calloc(mp->ntor, sizeof(int));

  /*** Array to store the torsion terms in a series ***/
  sterm = (int*)malloc(mp->ntor*sizeof(int));

  /*** Prepare a file to plot each series ***/
  outp = fopen(mp->sr, "w");
  ct = time(&ct);
  fprintf(outp, "%%\n%% Printed by mdgx on %s", asctime(localtime(&ct)));
  fprintf(outp, "%% This MatLab-readable results file plots torsion "
	  "potentials fitted\n%% according to input found in %s.\n%%\n\n",
	  tj->inpname);

  /*** Loop over all fitted torsion parameters ***/
  nser = 0;
  for (i = 0; i < mp->ntor; i++) {
    if (examined[i] == 1 || mp->torsions[i].fitcol < 0) {
      continue;
    }

    /*** Mark this torsion term ***/
    examined[i] = 1;    
    sterm[0] = i;
    ntt = 1;
    totsamp = mp->torsions[i].ninst;

    /*** Determine whether this torsion is part of a series ***/
    for (j = i+1; j < mp->ntor; j++) {
      if (TypeCompare(mp->torsions[i].atype, mp->torsions[i].btype,
		      mp->torsions[i].ctype, mp->torsions[i].dtype,
		      mp->torsions[j].atype, mp->torsions[j].btype,
		      mp->torsions[j].ctype, mp->torsions[j].dtype,
		      4, 1) == 2) {
	examined[j] = 1;
	sterm[ntt] = j;
	ntt++;
	totsamp += mp->torsions[j].ninst;
      }
    }

    /*** Print the output file ***/
    nser++;
    fprintf(outp, "%% Series for atom types %.4s %.4s %.4s %.4s\n"
	    "%% %d samples, %d Fourier terms\n",
	    mp->torsions[i].atype, mp->torsions[i].btype,
	    mp->torsions[i].ctype, mp->torsions[i].dtype, totsamp, ntt);
    fprintf(outp, "FSeries%d = [\n", nser);
    for (hangle = 0.0; hangle < 2.0*PI; hangle += 0.001) {
      tfunc = 0.0;
      for (j = 0; j < ntt; j++) {
	tortmp = &mp->torsions[sterm[j]];
	sangle = tortmp->pn*hangle - tortmp->phase;
	tfunc += tortmp->fitcol * (1.0 + cos(sangle));
      }
      fprintf(outp, "%9.4lf\n", tfunc);
    }
    fprintf(outp, "];\n");
  }

  /*** Close output file ***/
  fclose(outp);

  /*** Free allocated memory ***/
  free(examined);
  free(sterm);
}

/***=======================================================================***/
/*** SystemLOO: leave-one-out cross-validation of the matrix by excluding  ***/
/***            systems, one by one, from the fit.  If certain parameters  ***/
/***            are found only in one system then those parameters will    ***/
/***            also be excluded from the fit.                             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   A:        the fitting matrix, in its original form before the QR    ***/
/***             decomposition                                             ***/
/***   b:        the original target vector, Ax = b                        ***/
/***   x:        the original solution, containing all fitted parameters   ***/
/***             based on all systems                                      ***/
/***=======================================================================***/
static void SystemLOO(prmset *mp, dmat *A, double* b, double* x)
{
  int h, i, j, nc, nr, found;
  int* colcorr;
  dmat An;
  double* bn;

  /*** Column correspondence key (in case  ***/
  /*** some parameters must be eliminated) ***/
  colcorr = (int*)malloc(mp->nparm*sizeof(int));

  /*** Loop over all systems and make new fitting matrices ***/
  for (h = 0; h < mp->nunisys; h++) {

    /*** Assuming no eliminations, it is simple ***/
    /*** to determine the correspondence.  If   ***/
    /*** there are elimiations, this loop will  ***/
    /*** find them.  The result is colcorr      ***/
    /*** stating which columns of the original  ***/
    /*** fitting matrix the columns of the new  ***/
    /*** fitting matrix correspond to.          ***/
    nc = 0;
    for (i = 0; i < mp->nparm; i++) {
      colcorr[nc] = i;
      found = 0;
      for (j = 0; j < mp->nconf; j++) {
	if (mp->conf[i].GroupNum != h && fabs(A->map[j][i]) > 1.0e-8) {
	  found = 1;
	  break;
	}
      }
      if (found == 1) {
	nc++;
      }
    }

    /*** Construct the new fitting matrix ***/
    An = CreateDmat(A->row - mp->GroupCount[h], nc, 0);
    nr = 0;
    for (i = 0; i < mp->nconf; i++) {
      if (mp->conf[i].GroupNum == h) {
	continue;
      }
      for (j = 0; j < nc; j++) {
	An.map[nr][j] = A->map[i][colcorr[j]];
      }
      bn[nr] = b[i];
      nr++;
    }
  }
}

/***=======================================================================***/
/*** CountFittedTerms: count the number of fitted terms in a system.       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   sysno:    the system number                                         ***/
/***   mp:       the master parameter and fitting set                      ***/
/***=======================================================================***/
static int CountFittedTerms(int sysno, prmset *mp)
{
  int i, exid, nparm;
  prmtop *tp;
  mmsys *exconf;

  /*** First, find one example of this system ***/
  exid = -1;
  for (i = 0; i < mp->nconf; i++) {
    if (mp->conf[i].GroupNum == sysno) {
      exid = i;
      break;
    }
  }
  if (exid == -1) {
    printf("CountFittedTerms >> Error.  System %d does not exist.\n", sysno);
    exit(1);
  }
  exconf = &mp->conf[exid];
  tp = &exconf->tp;

  /*** Bonds, angles, and dihedrals ***/
  nparm = 0;
  for (i = 0; i < exconf->bmap.nbond; i++) {
    if (mp->bonds[exconf->bmap.id[i].key].fitcol >= 0) {
      nparm++;
    }
  }
  for (i = 0; i < exconf->amap.nangl; i++) {
    if (mp->angls[exconf->amap.id[i].key].fitcol > 0) {
      nparm++;
    }
  }
  for (i = 0; i < exconf->hmap.ntterm; i++) {
    if (mp->torsions[exconf->hmap.id[i].key].fitcol > 0) {
      nparm++;
    }
  }

  return nparm;
}

/***=======================================================================***/
/*** PrintAccuracy: report the accuracy of this model in relation to the   ***/
/***                training set.                                          ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   A:        the fitting matrix, in its original form before the QR    ***/
/***             decomposition                                             ***/
/***   b:        the target vector, Ax = b                                 ***/
/***   x:        the solution vector, containing all fitted parameters     ***/
/***   tj:       trajectory control information                            ***/
/***   outp:     the output file                                           ***/
/***=======================================================================***/
void PrintAccuracy(prmset *mp, dmat *A, double* b, double* x, trajcon *tj,
		   FILE *outp)
{
  int i, j, k, m, ncfg, fsys, natom;
  int* nsysvar;
  double fitnrg, eQQ, eLJ;
  double *atmp;
  double* enorm;
  double* eorig;
  double* efin;
  double* eelec;
  double* elj;
  double* ebond;
  double* eangl;
  double* edihe;
  double* efit;
  FILE *mscript;
  time_t ct;

  /*** If a MatLab readable display is requested, print it ***/
  if (mp->ao[0] != '\0') {
    mscript = FOpenSafe(mp->ao, tj->OverwriteOutput);
    ct = time(&ct);
    fprintf(mscript, "%% Model accuracy report for fitting ordered by\n"
	    "%% %s.\n", tj->inpname);
    fprintf(mscript, "%% Written on %s\n", asctime(localtime(&ct)));
  }
  fprintf(outp, 
	  "                 RMS Error         Correlation    Fitted   Fitted  "
	  " Model\n"
          "    System     Orig    Fitted    Orig    Fitted   Terms    Energy  "
	  " Diff.\n"
          " -----------  -------  -------  -------  -------  ------  -------  "
	  "-------\n");

  /*** Count the number of fitted terms in each system ***/
  nsysvar = (int*)malloc(mp->nunisys*sizeof(int));
  for (i = 0; i < mp->nunisys; i++) {
    nsysvar[i] = CountFittedTerms(i, mp);
  }

  /*** Compute the accuracy of the model ***/
  enorm = (double*)malloc(mp->nconf*sizeof(double));
  eorig = (double*)malloc(mp->nconf*sizeof(double));
  efin = (double*)malloc(mp->nconf*sizeof(double));
  eelec = (double*)malloc(mp->nconf*sizeof(double));
  elj = (double*)malloc(mp->nconf*sizeof(double));
  ebond = (double*)malloc(mp->nconf*sizeof(double));
  eangl = (double*)malloc(mp->nconf*sizeof(double));
  edihe = (double*)malloc(mp->nconf*sizeof(double));
  efit = (double*)calloc(mp->nconf, sizeof(double));
  for (i = 0; i < mp->nunisys; i++) {
    ncfg = 0;
    for (j = 0; j < mp->nconf; j++) {
      if (mp->conf[j].GroupNum != i) {
        continue;
      }
      if (mp->ao[0] != '\0' && ncfg == 0) {
	fprintf(mscript, "%% System %s\ndata%d = [\n", mp->conf[j].tp.source,
		i);
      }
      fsys = j;
      fitnrg = 0.0;
      atmp = A->map[j];
      for (k = 0; k < A->col; k++) {
        fitnrg += atmp[k] * x[k];
      }
      ebond[ncfg] = DSum(mp->conf[j].bmap.contrib, mp->conf[j].bmap.nbond);
      eangl[ncfg] = DSum(mp->conf[j].amap.contrib, mp->conf[j].amap.nangl);
      edihe[ncfg] = DSum(mp->conf[j].hmap.contrib, mp->conf[j].hmap.ntterm);
      natom = mp->conf[j].tp.natom;
      eLJ = 0.0;
      eQQ = 0.0;
      for (k = 0; k < natom-1; k++) {
	for (m = k+1; m < natom; m++) {
	  eQQ += mp->conf[j].nbnrg.map[k][m] * mp->conf[j].excl.map[k][m];
	  eLJ += mp->conf[j].nbnrg.map[m][k] * mp->conf[j].excl.map[m][k];
	}
      }
      eelec[ncfg] = eQQ;
      elj[ncfg] = eLJ;
      efin[ncfg] = mp->conf[j].nonfitmm + fitnrg;
      enorm[ncfg] = mp->conf[j].enorm;
      eorig[ncfg] = mp->conf[j].eorig;
      efit[i] += fitnrg;
      ncfg++;
    }
    efit[i] /= ncfg;
    if (mp->ao[0] != '\0') {
      for (j = 0; j < ncfg; j++) {
	fprintf(mscript, "%9.4lf %9.4lf %9.4lf\n", enorm[j], eorig[j],
		efin[j]);
      }
      fprintf(mscript, "];\n\ncomponents%d = [", i);
      fprintf(mscript, "  %% elec      LJ        bond      angl      dihe\n");
      for (j = 0; j < ncfg; j++) {
	fprintf(mscript, "%9.4lf %9.4lf %9.4lf %9.4lf %9.4lf\n", eelec[j],
		elj[j], ebond[j], eangl[j], edihe[j]);
      }
      fprintf(mscript, "];\nfigure; hold on;\n"
	      "plot(data%d(:,1), data%d(:,2), 'k.');\n"
	      "plot(data%d(:,1), data%d(:,3), 'r.');\n", i, i, i, i);
      fprintf(mscript, "legend('Original', 'Fitted');\n"
	      "xlabel('Target energy');\nylabel('Model energy');\n");
      fprintf(mscript, "title('Model results for %s');\n",
	      mp->conf[fsys].tp.source);
      fprintf(mscript, "n = input('Ready?');\n\n");
    }
    fprintf(outp, " %.11s  %7.4lf  %7.4lf  %7.4lf  %7.4lf  %6d  %7.3lf  "
	    "%7.3lf\n", mp->conf[fsys].tp.source, VecRMSD(eorig, enorm, ncfg),
            VecRMSD(efin, enorm, ncfg), Pearson(eorig, enorm, ncfg),
            Pearson(efin, enorm, ncfg), nsysvar[i], efit[i],
	    VecRMSD(efin, eorig, ncfg));
  }

  /*** Complete the MatLab output if requested ***/
  if (mp->ao[0] != '\0') {  
    fclose(mscript);
  }

  /*** Free allocated memory ***/
  free(enorm);
  free(eorig);
  free(efin);
  free(ebond);
  free(eangl);
  free(edihe);
  free(efit);
  free(elj);
  free(eelec);
  free(nsysvar);
}

/***=======================================================================***/
/*** ParamReport: report the best parameters and their fit to the target   ***/
/***              data.  New parameters are reported in a frcmod-format    ***/
/***              file, while statistics from the run are reported in a    ***/
/***              text file.                                               ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains atom types and examples of ***/
/***             each atom in various topologies)                          ***/
/***   A:        the fitting matrix, in its original form before the QR    ***/
/***             decomposition                                             ***/
/***   b:        the target vector, Ax = b                                 ***/
/***   x:        the solution vector, containing all fitted parameters     ***/
/***   tj:       trajectory control information                            ***/
/***=======================================================================***/
void ParamReport(prmset *mp, dmat *A, double* b, double* x, trajcon *tj)
{
  int i, j;
  FILE *outp;
  time_t ct;

  /*** Input file header ***/
  ct = time(&ct);
  outp = FOpenSafe(tj->outbase, tj->OverwriteOutput);
  PrintSplash(outp);
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

  /*** Print the energies according to the derived model ***/
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(1.)", 4, " ", 1, "Energies of each system, according to "
              "the fitted parameters.  Units on errors are kcal/mol.  The "
              "original model's energies are compared to target energies "
              "after applying an offset to make the average of each set of "
              "energies for any given system equal.\n", 74, 0, outp);
  PrintAccuracy(mp, A, b, x, tj, outp);
  HorizontalRule(outp, 1);

  /*** Print the sampling in each fitted ***/
  /*** parameter, across all systems     ***/
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(2.)", 4, " ", 1, "Parameter sampling.  Bonds and angles "
              "are binned over 18 intervals based on deviations, positive or "
              "negative, from the ideal bond length.  Each interval "
              "represents 0.03 Angstrom or 0.05 radian deviation for bonds or "
              "angles, respectively.  Dihedrals and impropers are binned at "
              "10-degree intervals, ranging from zero to 180 (the range of "
              "the arccosine function).\n", 74, 0, outp);
  fprintf(outp, " Sampling histogram key:  (zero)   . : = e o U O 0 @ X "
          "( > ten)\n\n");
  PrintSamplingTable(mp, 2, x, outp);
  PrintSamplingTable(mp, 3, x, outp);
  PrintSamplingTable(mp, 4, x, outp);
  HorizontalRule(outp, 1);
  HorizontalRule(outp, 1);
  fprintf(outp, "(3.) Fitting matrix analysis:\n\n");
  PrintMatrixAnalysis(mp, x, outp);
  HorizontalRule(outp, 1);
  HorizontalRule(outp, 1);
  ParameterDescription(mp, outp);
  HorizontalRule(outp, 1);
  fclose(outp);

  /*** FIX ME!!!  The CheckMinima function may not work so well. ***/
  CheckMinima(mp, tj, x);

  /*** Print the frcmod file for this parameter set ***/
  outp = FOpenSafe(tj->dumpname, tj->OverwriteOutput);
  fprintf(outp, "%s\n", mp->ititl);
  PrintParmAtoms(mp, outp);
  if (mp->Hydrophilics.row > 0) {
    for (i = 0; i < mp->Hydrophilics.row; i++) {
      fprintf(outp, "%.4s", mp->Hydrophilics.map[i]);
    }
    fprintf(outp, "\n");
  }
  PrintParmBond(mp, 2, x, outp);
  PrintParmBond(mp, 3, x, outp);
  PrintParmBond(mp, 4, x, outp);
  PrintParmBond(mp, 5, x, outp);
  PrintParmBond(mp, 9, x, outp);
  for (i = 0; i < mp->neqgroups; i++) {
    for (j = 0; j < mp->eqgroups[i].natom; j++) {
      fprintf(outp, "%.4s", &mp->eqgroups[i].types[4*j]);
    }
    fprintf(outp, "\n");
  }
  fprintf(outp, "\nMOD4      RE\n");
  for (i = 0; i < mp->natom; i++) {
    fprintf(outp, "  %.2s             %10.6lf  %10.6lf  %s\n",
	    mp->atoms[i].atype, mp->atoms[i].ljsig, mp->atoms[i].ljeps,
	    mp->atoms[i].comment);
  }
  fprintf(outp, "\nEND\n");
  fclose(outp);
}
