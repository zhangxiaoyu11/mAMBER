#include "stringDefs.h"
#include "matrix.h"

#ifndef PDBREAD_STRUCTS
#define PDBREAD_STRUCTS

/***=======================================================================***/
/*** PDBMOL: the format for molecules, based on the standard PDB format.   ***/
/***=======================================================================***/
struct pdbmol {
  int n_atoms;
  int n_res;
  int n_ter;
  int pqr;
  int aug;
  int AnisouRec;
  int* res_nums;
  int* atom_nums;
  int* res_lims;
  int* termini;
  int* z;
  int* anisou;
  char* chain;
  char* atom_names;
  char* res_names;
  double* crds;
  double* charges;
  double* radii;
  double* mass;
  char source[MAXNAME];
};
typedef struct pdbmol pdb;

/***=======================================================================***/
/*** SYMSTRUCT: the symmetries of the crystal space group, extracted from  ***/
/***            the PDB file.                                              ***/
/***=======================================================================***/
struct symstruct {
  int nsym;
  dmat* tmat;
  double* tvec;
  double boxd[6];
};
typedef struct symstruct symT;

#endif

#ifndef PDBREAD_FUNCS
#define PDBREAD_FUNCS

void GetPDB(pdb *thispdb, int verbosity, int pqr_format);

void PutPDB(pdb *thispdb, char* filename, char* pdb_form,
	    char* custom_header, int verbosity);

void FreePDB(pdb *thispdb);

char* sscan_adv(char *ctmp, char* cword);

void touchup(char* a, int asize);

void remove_whitespace(char* a, int asize);

void append_whitespace(char* a, int asize);

void alphanumeric(char *a);

int isalphanum(char a);

void realign_name(char* atom_name);

void CopyAtoms(pdb *add, int astart, int aend, pdb *base, int bstart, int chn);

void DeleteResidue(pdb *mol, int r);

void ModPdbRA(pdb *w);

void AugmentPDB(pdb *thispdb);

symT LoadSymmetryData(char* source);

pdb ExpandPDBSym(pdb *p, symT *S, int stride);

pdb TilePDB(pdb *p, double* boxd, int* irep, int stride);

#endif
