#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pdbRead.h"
#include "vector.h"
#include "crdmanip.h"
#include "stringDefs.h"
#include "matrix.h"
#include "myconstants.h"
#include "ran2.h"
#include "grid.h"

/***=======================================================================***/
/*** PFrac: scan the entire grid and determine the "protein fraction," the ***/
/***        portion of the grid that is masked off within RP of some       ***/
/***        protein atom.                                                  ***/
/***=======================================================================***/
double PFrac(bgrid *tg)
{
  int i, j, k;
  double dx;

  dx = 0.0;
  for (i = 0; i < tg->lx; i++) {
    for (j = 0; j < tg->ly; j++) {
      for (k = 0; k < tg->lz; k++) {
	dx += ReadBoolGrid(*tg, i, j, k);
      }
    }
  }
  dx /= tg->lx;
  dx /= tg->ly;
  dx /= tg->lz;

  return dx;
}

/***=======================================================================***/
/*** RunBGrid: this routine iteratively runs the grid coloring algorithm   ***/
/***           for the primary image and all of its 26 nearest neighbors.  ***/
/***=======================================================================***/
void RunBGrid(bgrid tg, double* crds, int P, dmat U, dmat invU, double RP,
	      int update_user)
{

  int i, j, k, m, m3;
  double* crdcpy;

  crdcpy = CpyDVec(crds, 3*P);
  BiColorGrid(tg, crds, P, U, invU, RP, update_user);
  for (i = -1; i <= 1; i++) {
    for (j = -1; j <= 1; j++) {
      for (k = -1; k <= 1; k++) {
	if (i != 0 || j != 0 || k != 0) {
	  for (m = 0; m < P; m++) {
	    m3 = 3*m;
	    crdcpy[m3] = crds[m3] + i;
	    crdcpy[m3+1] = crds[m3+1] + j;
	    crdcpy[m3+2] = crds[m3+2] + k;
	  }
	  BiColorGrid(tg, crdcpy, P, U, invU, RP, 0);
	}
      }
    }
  }

  free(crdcpy);
}

/***=======================================================================***/
/*** main                                                                  ***/
/***=======================================================================***/
int main(int argc, char* argv[])
{
  int h, i, j, N, P, placed, wx, wy, wz, nfail, nplaced, RCS, noOutPut;
  long counter;
  double RW, RP, R2P, R2W, gspc;
  double omx, omy, omz, dx, dy, dz;
  double boxd[6], ccm[3], acm[3], tvec[3];
  double* acpy;
  char tag[MAXNAME], syscall[MAXHEADER], a2bpath[MAXNAME];
  pdb c, a, o;
  dmat U, invU, M;
  bgrid tg;

  /*** Options ***/
  if (argc == 1) {
    printf("\nAddToBox >> A program for adding solvent molecules to a crystal "
	   "cell.\n\nOptions:\n"
	   "  -c    : the molecule cell (PDB format)\n"
	   "  -a    : the molecule to add\n"
	   "  -na   : the number of copies to add\n"
           "  -P    : the upper limit of protein atoms\n"
	   "  -o    : output file (PDB format)\n"
           "  -X    : X dimension of the box\n"
           "  -Y    : Y dimension of the box\n"
           "  -Z    : Z dimension of the box\n"
           "  -al   : a dimension of the box\n"
           "  -bt   : b dimension of the box\n"
           "  -gm   : g dimension of the box\n"
           "  -RW   : Clipping radius for solvent atoms\n"
	   "  -RP   : Clipping radius for protein atoms\n"
	   "  -IG   : Random number seed\n"
	   "  -NO   : flag for no PDB output (stops after determining the "
	   "protein\n          fraction of the box)\n"
	   "  -G    : Grid spacing for search (default 0.2)\n"
	   "  -V    : Recursively call AddToBox until all residues have been "
	   "added.\n"
	   "          (Default 0 ; any other setting activates recursion)\n"
	   "  -path : Path for AddToBox program on subsequent calls\n"
	   "          (default ${AMBERHOME}/bin/AddToBox)\n\n");
    exit(1);
  }

  /*** Initialize random number generator ***/
  counter = -1;
  ran2(&counter);

  /*** Input ***/
  c.source[0] = '\0';
  a.source[0] = '\0';
  o.source[0] = '\0';
  sprintf(a2bpath, "${AMBERHOME}/bin/AddToBox");
  boxd[0] = -1.0;
  boxd[1] = -1.0;
  boxd[2] = -1.0;
  boxd[3] = 90.0;
  boxd[4] = 90.0;
  boxd[5] = 90.0;
  RW = 1.0;
  RP = 5.0;
  RCS = 0;
  P = 0;
  noOutPut = 0;
  gspc = 0.2;
  for (i = 0; i < argc-2; i += 2) {
    strcpy(tag, *++argv);
    if (strcmp(tag, "-c") == 0) {
      strcpy(c.source, *++argv);
    }
    else if (strcmp(tag, "-a") == 0) {
      strcpy(a.source, *++argv);
    }
    else if (strcmp(tag, "-o") == 0) {
      strcpy(o.source, *++argv);
    }
    else if (strcmp(tag, "-na") == 0) {
      N = atoi(*++argv);
    }
    else if (strcmp(tag, "-P") == 0) {
      P = atoi(*++argv);
    }
    else if (strcmp(tag, "-X") == 0) {
      boxd[0] = atof(*++argv);
    }
    else if (strcmp(tag, "-Y") == 0) {
      boxd[1] = atof(*++argv);
    }
    else if (strcmp(tag, "-Z") == 0) {
      boxd[2] = atof(*++argv);
    }
    else if (strcmp(tag, "-al") == 0) {
      boxd[3] = atof(*++argv);
    }
    else if (strcmp(tag, "-bt") == 0) {
      boxd[4] = atof(*++argv);
    }
    else if (strcmp(tag, "-gm") == 0) {
      boxd[5] = atof(*++argv);
    }
    else if (strcmp(tag, "-RW") == 0) {
      RW = atof(*++argv);
    }
    else if (strcmp(tag, "-RP") == 0) {
      RP = atof(*++argv);
    }
    else if (strcmp(tag, "-G") == 0) {
      gspc = atof(*++argv);
    }
    else if (strcmp(tag, "-NO") == 0) {
      noOutPut = atoi(*++argv);
    }
    else if (strcmp(tag, "-V") == 0) {
      RCS = atoi(*++argv);
    }
    else if (strcmp(tag, "-IG") == 0) {
      counter = atoi(*++argv);
      ran2(&counter);
    }
    else if (strcmp(tag, "-path") == 0) {
      strcpy(a2bpath, *++argv);
    }
    else {
      printf("AddToBox >> Error.  Identifier %s unkonwn.\n\n", tag);
      exit(1); 
    }
  }

  /*** Check ***/
  if (c.source[0] == '\0') {
    printf("AddToBox >> Error.  Protein file not specified!\n");
    exit(1);
  }
  if (a.source[0] == '\0') {
    printf("AddToBox >> Error.  Additional molecule file not specified!\n");
    exit(1);
  }
  if (o.source[0] == '\0') {
    printf("AddToBox >> Error.  Output file not specified!\n");
    exit(1);
  }
  for (i = 0; i < 3; i++) {
    if (boxd[i] < 0.0) {
      printf("AddToBox >> Error!  Box dimension wrong!\n");
      printf("AddToBox >> %8.3lf %8.3lf %8.3lf\n", boxd[0], boxd[1], boxd[2]);
      exit(1);
    }
    boxd[i+3] *= PI/180.0;
  }
  R2W = RW*RW;
  R2P = RP*RP;

  /*** Get PDB inputs ***/
  GetPDB(&c, 1, 0);
  GetPDB(&a, 1, 0);

  /*** Prepare output ***/
  o.n_atoms = c.n_atoms + N*a.n_atoms;
  o.n_res = c.n_res + N*a.n_res;
  o.n_ter = c.n_ter + N*(a.n_ter+1);
  o.pqr = 0;
  o.crds = (double*)malloc(3*o.n_atoms*sizeof(double));
  o.atom_names = (char*)malloc(4*o.n_atoms*sizeof(char));
  o.atom_nums = (int*)malloc(o.n_atoms*sizeof(int));
  o.res_names = (char*)malloc(4*o.n_atoms*sizeof(char));
  o.res_nums = (int*)malloc(o.n_atoms*sizeof(int));
  o.chain = (char*)malloc(o.n_atoms*sizeof(char));
  o.termini = (int*)malloc(o.n_ter*sizeof(int));
  o.res_lims= (int*)malloc((o.n_res+1)*sizeof(int));
  o.res_lims[0] = 0;
  for (i = 0; i < c.n_res; i++) {
    o.res_lims[i] = c.res_lims[i];
  }
  for (i = 0; i < c.n_ter; i++) {
    o.termini[i] = c.termini[i];
  }
  for (i = 0; i < N; i++) {
    for (j = 0; j < a.n_res; j++) {
      o.res_lims[c.n_res+i*a.n_res+j] =
	c.n_atoms + i*a.n_atoms + a.res_lims[j];
    }
    for (j = 0; j < a.n_ter; j++) {
      o.termini[c.n_ter+i*(a.n_ter+1)+j] =
	c.n_atoms + i*a.n_atoms + a.termini[j];
    }
    o.termini[c.n_ter+i*(a.n_ter+1)+a.n_ter] = c.n_atoms + i*a.n_atoms;
  }
  o.res_lims[o.n_res] = o.n_atoms;
  CopyAtoms(&c, 0, c.n_atoms, &o, 0, 0);

  /*** Transform into crystal lattice space ***/
  U = CreateDmat(3, 3);
  invU = CreateDmat(3, 3);
  M = CreateDmat(3, 3);
  CmpXfrm(boxd, U, invU);
  RotateCrd(c.crds, c.n_atoms, U);

  /*** Center the protein and move solvent accordingly ***/
  FindCrdCenter(c.crds, c.crds, 0, c.n_atoms, ccm);
  TransCrd(c.crds, c.n_atoms, ccm, -1.0);
  FindCrdCenter(a.crds, a.crds, 0, a.n_atoms, acm);
  TransCrd(a.crds, a.n_atoms, acm, -1.0);

  /*** Re-image the protein ***/
  ReImage(c.crds, c.n_atoms);

  /*** Now, make the grid ***/
  tg = CreateBoolGrid(boxd[0]/gspc+1, boxd[1]/gspc+1, boxd[2]/gspc+1,
		      -0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  RunBGrid(tg, c.crds, P, U, invU, RP, 1);

  /*** Report the amount of "protein" volume ***/
  if (noOutPut == 1) {
    printf("\nAddToBox >> Protein fraction: %.4lf\n", PFrac(&tg));
    exit(1);
  }

  RunBGrid(tg, &c.crds[3*P], c.n_atoms-P, U, invU, RW, 1);

  /*** Now examine all residues ***/
  nfail = 0;
  for (i = 0; i < N; i++) {

    fprintf(stderr, "\rAddToBox >> Examining residue %7d", i);
    fflush(stderr);
                                     
    /*** Pointer to oxygen coordinates ***/
    placed = 0;
    while (placed == 0) {

      /*** Rotate randomly ***/
      FindCrdCenter(a.crds, a.crds, 0, a.n_atoms, acm);
      TransCrd(a.crds, a.n_atoms, acm, -1.0);
      omx = ran2(&counter)*PI;
      omy = ran2(&counter)*PI;
      omz = ran2(&counter)*PI;
      BeardRotMat(omx, omy, omz, M);
      RotateCrd(a.crds, a.n_atoms, M);
      dx = ran2(&counter) - 0.5;
      dy = ran2(&counter) - 0.5;
      dz = ran2(&counter) - 0.5;
      tvec[0] = invU.data[0]*dx + invU.data[1]*dy + invU.data[2]*dz;
      tvec[1] = invU.data[3]*dx + invU.data[4]*dy + invU.data[5]*dz;
      tvec[2] = invU.data[6]*dx + invU.data[7]*dy + invU.data[8]*dz;
      TransCrd(a.crds, a.n_atoms, tvec, 1.0);
      acpy = CpyDVec(a.crds, 3*a.n_atoms);
      RotateCrd(acpy, a.n_atoms, U);
      ReImage(acpy, a.n_atoms); 
      placed = 1;
      for (h = 0; h < a.n_atoms; h++) {

	/*** Find the grid cell ***/
	wx = (acpy[3*h] - tg.ox)*tg.invgx;
	wy = (acpy[3*h+1] - tg.oy)*tg.invgy;
	wz = (acpy[3*h+2] - tg.oz)*tg.invgz;
	if (ReadBoolGrid(tg, wx, wy, wz) == 1) {
	  placed = 0;
	}
      }
      if (placed == 1) {

	/*** Reset failure count ***/
	nfail = 0;

	/*** Update grid ***/
	RunBGrid(tg, acpy, a.n_atoms, U, invU, RW, 0);

	/*** Accumulate molecules ***/
	for (h = 0; h < a.n_atoms; h++) {
	  a.res_nums[h] = i;
	}
	CopyAtoms(&a, 0, a.n_atoms, &o, c.n_atoms + i*a.n_atoms, 0);
      }
      else {
	nfail++;
	FindCrdCenter(a.crds, a.crds, 0, a.n_atoms, acm);
	TransCrd(a.crds, a.n_atoms, acm, -1.0);
	BeardRotMat(-omx, -omy, -omz, M);
	RotateCrd(a.crds, a.n_atoms, M);
      }
      free(acpy);

      /*** If we can't place it, break ***/
      if (nfail > 1000000) {
	break;
      }
    }

    /*** If we can not place anymore residues, break ***/
    if (nfail > 1000000) {
      o.n_atoms = c.n_atoms + i*a.n_atoms;
      break;
    }
  }
  nplaced = i;

  /*** Get back into real space ***/
  for (i = 0; i < 3; i++) {
    acm[i] = DotP(&invU.data[3*i], ccm, 3);
  }
  TransCrd(&o.crds[3*c.n_atoms], o.n_atoms-c.n_atoms, acm, 1.0);

  /*** Make sure residue numbers are all right ***/
  for (i = 0; i < o.n_atoms; i++) {
    if (o.res_nums[i] >= 10000) {
      o.res_nums[i] = o.res_nums[i] % 10000;
    }
  }
  ModPdbRA(&o);
  PutPDB(&o, o.source, "STANDARD", "HEADER  ", 0);

  /*** Free allocated memory ***/
  DestroyDmat(&U);
  DestroyDmat(&invU);
  DestroyDmat(&M);
  free(tg.data);

  printf("\nAddToBox >> Added %d residues.\n", nplaced);
  if (nplaced == N) {
    printf("AddToBox >> Done.\n\n");
  }
  else if (RCS == 1) {
    printf("AddToBox >> Residue goal not achieved.\n"
	   "AddToBox >> Re-attempting with solvent clipping radius %8.4lf.\n",
	   RW*0.90);
    sprintf(syscall, "%s -c %s -o %s -RW %lf -RP %lf -IG %d -G %lf -a %s -na "
	    "%d -X %lf -Y %lf -Z %lf -al %lf -bt %lf -gm %lf -P %d -V 1",
	    a2bpath, o.source, o.source, RW*0.90, RP,
	    (int)(ran2(&counter)*1.0e6), gspc, a.source, N - nplaced, boxd[0],
	    boxd[1], boxd[2], boxd[3]*180.0/PI, boxd[4]*180.0/PI,
	    boxd[5]*180.0/PI, P);
    system(syscall);
  }

  return 0;
}
