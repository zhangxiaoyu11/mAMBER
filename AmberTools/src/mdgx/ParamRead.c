#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Matrix.h"
#include "Topology.h"
#include "Parse.h"

#include "ParamFitDS.h"

/***=======================================================================***/
/*** Type4Char: fill a type definition string up to 4 characters, padding  ***/
/***            it with blank spaces as necessary.                         ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   a:     the type string to write                                     ***/
/***   t:     the type string input                                        ***/
/***=======================================================================***/
static void Type4Char(char* a, char* t)
{
  int i, pivot;

  strncpy(a, t, 4);
  pivot = 0;
  for (i = 0; i < 4; i++) {
    if (a[i] == '\0') {
      pivot = 1;
    }
    if (pivot == 1) {
      a[i] = ' ';
    }
  }
  a[4] = '\0';
}

/***=======================================================================***/
/*** CrossRefAtomType: function to find the index of a given atom type     ***/
/***                   among all of those stored in a parameter set.       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:             the parameter set                                   ***/
/***   aname:          the name of the atom type                           ***/
/***=======================================================================***/
int CrossRefAtomType(prmset *mp, char* aname)
{
  int i;
  char typeA[8];

  Type4Char(typeA, aname);
  for (i = 0; i < mp->natom; i++) {
    if (strncmp(mp->atoms[i].atype, typeA, 4) == 0) {
      return i;
    }
  }

  return -1;
}

/***=======================================================================***/
/*** CrossRefBondType: function to find the index of a given atom type     ***/
/***                   among all of those stored in a parameter set.       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:             the parameter set                                   ***/
/***   [a,b]name:      the name of the atom type                           ***/
/***=======================================================================***/
static int CrossRefBondType(prmset *mp, char* aname, char* bname)
{
  int i;
  char typeA[8], typeB[8];

  Type4Char(typeA, aname);
  Type4Char(typeB, bname);
  for (i = 0; i < mp->nbond; i++) {
    if ((strncmp(mp->bonds[i].atype, typeA, 4) == 0 &&
	 strncmp(mp->bonds[i].btype, typeB, 4) == 0) ||
	(strncmp(mp->bonds[i].atype, typeB, 4) == 0 &&
	 strncmp(mp->bonds[i].btype, typeA, 4) == 0)) {
      return i;
    }
  }

  return -1;
}

/***=======================================================================***/
/*** CrossRefHBondType: function to find the index of a given hydrogen     ***/
/***                    bond potential among all of those store in a       ***/
/***                    parameter set.                                     ***/
/***=======================================================================***/
static int CrossRefHBondType(prmset *mp, char* aname, char* bname)
{
  int i;
  char typeA[8], typeB[8];

  Type4Char(typeA, aname);
  Type4Char(typeB, bname);
  for (i = 0; i < mp->nhb1012; i++) {
    if ((strncmp(mp->hb1012[i].atype, typeA, 4) == 0 &&
         strncmp(mp->hb1012[i].btype, typeB, 4) == 0) ||
        (strncmp(mp->hb1012[i].atype, typeB, 4) == 0 &&
         strncmp(mp->hb1012[i].btype, typeA, 4) == 0)) {
      return i;
    }
  }

  return -1;
}

/***=======================================================================***/
/*** CrossRefAnglType: function to find the index of a given atom type     ***/
/***                   among all of those stored in a parameter set.       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:             the parameter set                                   ***/
/***   [a,b,c]name:    the name of the atom type                           ***/
/***=======================================================================***/
static int CrossRefAnglType(prmset *mp, char* aname, char* bname, char* cname)
{
  int i;
  char typeA[8], typeB[8], typeC[8];

  Type4Char(typeA, aname);
  Type4Char(typeB, bname);
  Type4Char(typeC, cname);
  for (i = 0; i < mp->nangl; i++) {
    if ((strncmp(mp->angls[i].atype, typeA, 4) == 0 &&
	 strncmp(mp->angls[i].btype, typeB, 4) == 0 && 
	 strncmp(mp->angls[i].ctype, typeC, 4) == 0) ||
	(strncmp(mp->angls[i].atype, typeC, 4) == 0 &&
	 strncmp(mp->angls[i].btype, typeB, 4) == 0 &&
	 strncmp(mp->angls[i].ctype, typeA, 4) == 0)) {
      return i;
    }
  }

  return -1;
}

/***=======================================================================***/
/*** CrossRefDiheType: function to find the index of a given atom type     ***/
/***                   among all of those stored in a parameter set.       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:             the parameter set                                   ***/
/***   [a,b,c,d]name:  the name of the atom type                           ***/
/***=======================================================================***/
static int CrossRefDiheType(prmset *mp, cmat lwords, int order)
{
  int i;
  double K, phase, pn;
  char typeA[8], typeB[8], typeC[8], typeD[8];

  Type4Char(typeA, lwords.map[0]);
  Type4Char(typeB, lwords.map[1]);
  Type4Char(typeC, lwords.map[2]);
  Type4Char(typeD, lwords.map[3]);
  for (i = 0; i < mp->ntor; i++) {
    if ((strncmp(mp->torsions[i].atype, typeA, 4) == 0 && 
	 strncmp(mp->torsions[i].btype, typeB, 4) == 0 &&
	 strncmp(mp->torsions[i].ctype, typeC, 4) == 0 &&
	 strncmp(mp->torsions[i].dtype, typeD, 4) == 0) ||
	(strncmp(mp->torsions[i].atype, typeD, 4) == 0 &&
	 strncmp(mp->torsions[i].btype, typeC, 4) == 0 &&
	 strncmp(mp->torsions[i].ctype, typeB, 4) == 0 &&
	 strncmp(mp->torsions[i].dtype, typeA, 4) == 0)) {

      /*** The types match, but what about other     ***/
      /*** aspects of this fourier term declaration? ***/
      if (order == 4) {
	K = atof(lwords.map[5]) / atof(lwords.map[4]);
	phase = atof(lwords.map[6])*PI/180.0;
	pn = atof(lwords.map[7]);
      }
      else if (order == 5) {
	K = atof(lwords.map[4]);
	phase = atof(lwords.map[5])*PI/180.0;
	pn = atof(lwords.map[6]);
      }
      if (fabs(mp->torsions[i].K - K) < 1.0e-4 &&
	  fabs(mp->torsions[i].phase - phase) < 1.0e-4 &&
	  fabs(mp->torsions[i].pn - pn) < 1.0e-4) {
	return i;
      }
    }
  }

    return -1;
}

/***=======================================================================***/
/*** ParmFileComment: record a comment from a parameter file.              ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   lwords:    the list of words which may form a comment               ***/
/***   wstart:    the starting word                                        ***/
/***=======================================================================***/
static char* ParmFileComment(cmat *lwords, int wstart)
{
  int i, j;
  char* commtext;

  commtext = (char*)malloc(MAXLINE*sizeof(char));
  commtext[0] = '\0';
  if (wstart >= lwords->row) {
    return commtext;
  }
  j = 0;
  for (i = wstart; i < lwords->row; i++) {
    sprintf(&commtext[j], "%s ", lwords->map[i]);
    j = strlen(commtext);
  }

  return commtext;
}

/***=======================================================================***/
/*** RunOfWords: check to see that there is a run of atom types or numbers ***/
/***             in selected positions of a word list.                     ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   lwords:    the list of words                                        ***/
/***   istart:    the start position of the atom types run                 ***/
/***   iend:      the final position of the atom types run                 ***/
/***   style:     set to 0 for atom types, 1 for numbers                   ***/
/***=======================================================================***/
static int RunOfWords(cmat *lwords, int istart, int iend, int style)
{
  int i;

  for (i = istart; i <= iend; i++) {
    if (i >= 0 && i < lwords->row) {
      if (style == 0 && WordIsAtomType(lwords->map[i]) == 0) {
	return 0;
      }
      if (style == 1 && WordIsNumber(lwords->map[i]) == 0) {
	return 0;
      }
    }
  }

  return 1;
}

/***=======================================================================***/
/*** FormatTest: function to test the format of a line for matches to each ***/
/***             type of parameter file data.                              ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   lwords:    the line, converted to words                             ***/
/***   order:     integer denoting the type of input needed                ***/
/***=======================================================================***/
static int FormatTest(cmat *lwords, int order)
{
  /*** Atoms ***/
  if (order == 1 || order == 8) {
    if (lwords->row < 2) {
      return 0;
    }
    if (WordIsAtomType(lwords->map[0]) == 1 &&
	WordIsNumber(lwords->map[1]) == 1) {
      if (order == 8 && WordIsNumber(lwords->map[2]) == 0) {
	return 0;
      }
      return 1;
    }
  }

  /*** Bonds, hydrogen bonds ***/
  if (order == 2 || order == 9) {
    if (lwords->row < 4) {
      return 0;
    }
    if (RunOfWords(lwords, 0, 1, 0) == 1 && RunOfWords(lwords, 2, 3, 1) == 1) {
      return 1;
    }
  }

  /*** Angles ***/
  if (order == 3) {
    if (lwords->row < 5) {
      return 0;
    }
    if (RunOfWords(lwords, 0, 2, 0) == 1 && RunOfWords(lwords, 3, 4, 1) == 1) {
      return 1;
    }
  }

  /*** Dihedrals ***/
  if (order == 4) {
    if (lwords->row < 8) {
      return 0;
    }
    if (RunOfWords(lwords, 0, 3, 0) == 1 && RunOfWords(lwords, 4, 7, 1) == 1) {
      return 1;
    }
  }

  /*** Improper Dihedrals ***/
  if (order == 5) {
    if (lwords->row < 7) {
      return 0;
    }
    if (RunOfWords(lwords, 0, 3, 0) == 1 && RunOfWords(lwords, 4, 6, 1) == 1) {
      return 1;
    }
  }

  /*** Line of atom types ***/
  if (order == 6 || order == 7) {
    if (lwords->row == 0) {
      return 0;
    }
    if (RunOfWords(lwords, 0, lwords->row-1, 0) == 1) {
      return 1;
    }
  }

  /*** Nope, the line didn't conform ***/
  return 0;
}

/***=======================================================================***/
/*** ScanDeclarations: scan through the parameter or forcemod file for     ***/
/***                   lines that conform to the required format.          ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:      the parameter set                                          ***/
/***   cif:     the parameter input file, converted to a character matrix  ***/
/***   order:   the order of the comparison                                ***/
/***   nentry:  the number of entries already logged, to which this scan   ***/
/***            will add                                                   ***/
/***   strline: the line on which to start the scan                        ***/
/***=======================================================================***/
static int ScanDeclarations(prmset *mp, cmat *cif, int order, int nentry,
			    int strline, char* tercode)
{
  int i, j, k, m, llen, ndash, represent, entryID, newadd;
  int segter, formatmatch, matchfound, maxentry;
  char line[MAXLINE], typetmp[8];
  cmat lwords;

  /*** Allocate memory ***/
  maxentry = nentry + 32;
  if (nentry == 0) {
    if (order == 1) {
      mp->atoms = (xatomdef*)malloc(maxentry*sizeof(xatomdef));
    }
    else if (order == 2) {
      mp->bonds = (xbonddef*)malloc(maxentry*sizeof(xbonddef));
    }
    else if (order == 3) {
      mp->angls = (xangldef*)malloc(maxentry*sizeof(xangldef));
    }
    else if (order == 4) {
      mp->torsions = (torterm*)malloc(maxentry*sizeof(torterm));
    }
    else if (order == 6) {
      mp->Hydrophilics = CreateCmat(1, 5);
      mp->Hydrophilics.row = 0;
    }
    else if (order == 7) {
      mp->eqgroups = (eagrp*)malloc(maxentry*sizeof(eagrp));
    }
    else if (order == 9) {
      mp->hb1012 = (xhb1012def*)malloc(maxentry*sizeof(xhb1012def));
    }
  }
  else {
    if (order == 1) {
      mp->atoms = (xatomdef*)realloc(mp->atoms, maxentry*sizeof(xatomdef));
    }
    else if (order == 2) {
      mp->bonds = (xbonddef*)realloc(mp->bonds, maxentry*sizeof(xbonddef));
    }
    else if (order == 3) {
      mp->angls = (xangldef*)realloc(mp->angls, maxentry*sizeof(xangldef));
    }
    else if (order == 4 || order == 5) {
      mp->torsions = (torterm*)realloc(mp->torsions, maxentry*sizeof(torterm));
    }
    else if (order == 9) {
      mp->hb1012 = (xhb1012def*)realloc(mp->hb1012,
					maxentry*sizeof(xhb1012def));
    }
  }

  /*** Loop through the available portions of the input file ***/
  matchfound = 0;
  segter = 0;
  for (i = strline; i < cif->row; i++) {

    /*** Copy the line so it can be modified ***/
    strcpy(line, cif->map[i]);
    RemoveWhiteSpace(line, 256);
    llen = strlen(line);

    /*** Remove dashes ***/
    if (order >= 2 && order <= 5) {
      ndash = 0;
      for (j = 0; j < llen; j++) {
	if (line[j] == '-') {
	  ndash++;
	  line[j] = ' ';
	}

	/*** Break if the maximum number of ***/
	/*** dashes have been encountered   ***/
	if (ndash == order-1) {
	  break;
	}
      }
    }

    /*** Scan the line into words ***/
    lwords = ParseWords(line);

    /*** Determine whether the format matches ***/
    formatmatch = FormatTest(&lwords, order);

    /*** If the format does match, proceed accordingly ***/
    if (formatmatch == 1) {
      matchfound = 1;             
      if (order == 1 || order == 8) {
	entryID = CrossRefAtomType(mp, lwords.map[0]);
      }
      else if (order == 2) {
	entryID = CrossRefBondType(mp, lwords.map[0], lwords.map[1]);
      }
      else if (order == 3) {
	entryID = CrossRefAnglType(mp, lwords.map[0], lwords.map[1],
				   lwords.map[2]);
      }
      else if (order == 4 || order == 5) {
	entryID = CrossRefDiheType(mp, lwords, order);
      }
      else if (order == 9) {
	entryID = CrossRefHBondType(mp, lwords.map[0], lwords.map[1]);
      }
      else {
	entryID = -1;
      }
      if (entryID == -1) {
	newadd = 1;
	entryID = nentry;
      }
      else {
	newadd = 0;
      }
      if (order == 1) {
	Type4Char(mp->atoms[entryID].atype, lwords.map[0]);
	mp->atoms[entryID].mass = atof(lwords.map[1]);
	if (lwords.row > 2) {
	  if (WordIsNumber(lwords.map[2]) == 1) {
	    mp->atoms[entryID].apol = atof(lwords.map[2]);
	    mp->atoms[entryID].comment = ParmFileComment(&lwords, 3);
	  }
	  else {
	    mp->atoms[entryID].apol = 0.0;
	    mp->atoms[entryID].comment = ParmFileComment(&lwords, 2);
	  }
	}
	else {
	  mp->atoms[entryID].apol = 0.0;
	  mp->atoms[entryID].comment = ParmFileComment(&lwords, 2);
	}
	mp->atoms[entryID].ljsig = 0.0;
	mp->atoms[entryID].ljeps = 0.0;
      }
      else if (order == 2) {
        Type4Char(mp->bonds[entryID].atype, lwords.map[0]);
        Type4Char(mp->bonds[entryID].btype, lwords.map[1]);
	mp->bonds[entryID].K = atof(lwords.map[2]);
	mp->bonds[entryID].l0 = atof(lwords.map[3]);
	mp->bonds[entryID].comment = ParmFileComment(&lwords, 4);
      }
      else if (order == 3) {
        Type4Char(mp->angls[entryID].atype, lwords.map[0]);
        Type4Char(mp->angls[entryID].btype, lwords.map[1]);
        Type4Char(mp->angls[entryID].ctype, lwords.map[2]);
	mp->angls[entryID].K = atof(lwords.map[3]);
	mp->angls[entryID].th0 = atof(lwords.map[4])*PI/180.0;
	mp->angls[entryID].comment = ParmFileComment(&lwords, 5);
      }
      else if (order == 4 || order == 5) {
        Type4Char(mp->torsions[entryID].atype, lwords.map[0]);
        Type4Char(mp->torsions[entryID].btype, lwords.map[1]);
        Type4Char(mp->torsions[entryID].ctype, lwords.map[2]);
        Type4Char(mp->torsions[entryID].dtype, lwords.map[3]);
	if (order == 4) {
	  mp->torsions[entryID].K = atof(lwords.map[5]) / atof(lwords.map[4]);
	  mp->torsions[entryID].phase = atof(lwords.map[6])*PI/180.0;
	  mp->torsions[entryID].pn = atof(lwords.map[7]);
	  mp->torsions[entryID].impr = 0;
	  mp->torsions[entryID].comment = ParmFileComment(&lwords, 8);
	}
	else if (order == 5) {
	  mp->torsions[entryID].K = atof(lwords.map[4]);
	  mp->torsions[entryID].phase = atof(lwords.map[5])*PI/180.0;
	  mp->torsions[entryID].pn = atof(lwords.map[6]);
	  mp->torsions[entryID].impr = 1;
	  mp->torsions[entryID].comment = ParmFileComment(&lwords, 7);
	}
	if (mp->torsions[entryID].pn < 0.0) {
	  mp->torsions[entryID].singlet = -1;
	  mp->torsions[entryID].pn = fabs(mp->torsions[entryID].pn);
	}
	else {
	  mp->torsions[entryID].singlet = 1;
	}
      }
      else if (order == 6) {
	k = mp->Hydrophilics.row;
	mp->Hydrophilics = ReallocCmat(&mp->Hydrophilics, k + lwords.row, 5);
	for (j = 0; j < lwords.row; j++) {
	  Type4Char(mp->Hydrophilics.map[k+j], lwords.map[j]);
	}
      }
      else if (order == 7) {
	mp->eqgroups = (eagrp*)realloc(mp->eqgroups,
				       (entryID+1)*sizeof(eagrp));
	mp->eqgroups[entryID].natom = lwords.row;
	mp->eqgroups[entryID].types =
	  (char*)malloc((4*lwords.row+1)*sizeof(char));
        for (j = 0; j < lwords.row; j++) {
          Type4Char(&mp->eqgroups[entryID].types[4*j], lwords.map[j]);
        }
      }
      else if (order == 8) {

	/*** In this case, the objective is not to add ***/
	/*** to a growing list, but to contribute more ***/
	/*** information to what already exists.       ***/             
	mp->atoms[entryID].ljsig = atof(lwords.map[1]);
	mp->atoms[entryID].ljeps = atof(lwords.map[2]);
	Type4Char(typetmp, lwords.map[0]);
	for (j = 0; j < mp->neqgroups; j++) {
	  represent = 0;
	  for (k = 0; k < mp->eqgroups[j].natom; k++) {
	    if (strncmp(&mp->eqgroups[j].types[4*k], typetmp, 4) == 0) {
	      represent = 1;
	    }
	  }
	  if (represent == 1) {
	    for (k = 0; k < mp->eqgroups[j].natom; k++) {
	      m = CrossRefAtomType(mp, &mp->eqgroups[j].types[4*k]);
	      if (m == -1) {
		printf("ScanDeclarations >> Error.  Cross-reference to atom "
		       "type %.4s failed.\n", &mp->eqgroups[j].types[4*k]);
		exit(1);
	      }
	      mp->atoms[m].ljsig = atof(lwords.map[1]);
	      mp->atoms[m].ljeps = atof(lwords.map[2]);
	    }
	  }
	}
      }
      else if (order == 9) {
        Type4Char(mp->hb1012[entryID].atype, lwords.map[0]);
        Type4Char(mp->hb1012[entryID].btype, lwords.map[1]);
	mp->hb1012[entryID].Aterm = atof(lwords.map[2]);
	mp->hb1012[entryID].Bterm = atof(lwords.map[3]);
	mp->hb1012[entryID].comment = ParmFileComment(&lwords, 4);
      }

      /*** Increment the number of entries ***/
      if (newadd == 1) {
	nentry += 1;
	if (order == 1) {
	  mp->natom = nentry;
	}
	else if (order == 2) {
	  mp->nbond = nentry;
	}
	else if (order == 3) {
	  mp->nangl = nentry;
	}
	else if (order == 4 || order == 5) {
	  mp->ntor = nentry;
	}
	else if (order == 7) {
	  mp->neqgroups = nentry;
	}
	else if (order == 9) {
	  mp->nhb1012 = nentry;
	}
      }
      if (nentry >= maxentry) {
	maxentry += 32;
	if (order == 1) {
	  mp->atoms = (xatomdef*)realloc(mp->atoms, maxentry*sizeof(xatomdef));
	}
	else if (order == 2) {
	  mp->bonds = (xbonddef*)realloc(mp->bonds, maxentry*sizeof(xbonddef));
	}
	else if (order == 3) {
	  mp->angls = (xangldef*)realloc(mp->angls, maxentry*sizeof(xangldef));
	}
	else if (order == 4 || order == 5) {
	  mp->torsions = (torterm*)realloc(mp->torsions,
					   maxentry*sizeof(torterm));
	}
      }
    }
    else if (matchfound == 1) {

      /*** Consider the possibility that this is a segment terminator ***/
      if (strcmp(tercode, "BLANKLINE") == 0 && line[0] == '\n') {
	segter = 1;
      }
      if (strcmp(tercode, "ONELINE") == 0) {
	segter = 1;
      }
      if (strcmp(tercode, "END") == 0 &&
	  (strcmp(line, "END") == 0 || strcmp(line, "end") == 0)) {
	segter = 1;
      }
    }

    /*** Free allocated memory ***/
    DestroyCmat(&lwords);

    /*** Terminate the search ***/
    if (segter == 1) {
      break;
    }
  }

  /*** What about the next line of the file? ***/
  if (matchfound == 0) {
    return strline;
  }
  else {
    return i;
  }
}

/***=======================================================================***/
/*** ReadParmFile: read a force field parameter file, to help identify     ***/
/***               wildcard parameters and thus pare down the number of    ***/
/***               variables to fit.                                       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:      the fitting data (contains a list of all systems)          ***/
/***   tj:      trajectory control data (contains the frcmod file name)    ***/
/***=======================================================================***/
void ReadParmFile(prmset *mp, trajcon *tj)
{
  int nextline;
  cmat cfi;

  /*** Test file existence ***/
  if (tj->parmfile[0] == '\0') {
    printf("ReadParmFile >> Error.  Force field parameter file not specified."
           "\n");
    exit(1);
  }
  cfi = Ascii2Mem(tj->parmfile, 256, "Force field parameter file not found.");

  /*** Read the atom declarations ***/
  mp->natom = 0;
  nextline = ScanDeclarations(mp, &cfi, 1, mp->natom, 1, "BLANKLINE");

  /*** Record the hydrophilic atoms ***/
  nextline = ScanDeclarations(mp, &cfi, 6, 0, nextline, "ONELINE");

  /*** Read bond, angle, dihedral, and improper declarations ***/
  mp->nbond = 0;
  nextline = ScanDeclarations(mp, &cfi, 2, mp->nbond, nextline, "BLANKLINE");
  mp->nangl = 0;
  nextline = ScanDeclarations(mp, &cfi, 3, mp->nangl, nextline, "BLANKLINE");
  mp->ntor = 0;
  nextline = ScanDeclarations(mp, &cfi, 4, mp->ntor, nextline, "BLANKLINE");
  nextline = ScanDeclarations(mp, &cfi, 5, mp->ntor, nextline, "BLANKLINE");
  mp->nhb1012 = 0;
  nextline = ScanDeclarations(mp, &cfi, 9, mp->nhb1012, nextline, "BLANKLINE");

  /*** Read atom equivalencies ***/
  nextline = ScanDeclarations(mp, &cfi, 7, 0, nextline, "BLANKLINE");

  /*** Read non-bonded parameters and     ***/
  /*** cross-reference with equivalencies ***/
  nextline = ScanDeclarations(mp, &cfi, 8, 0, nextline, "BLANKLINE");
}

/***=======================================================================***/
/*** ReadFrcmodFile: read a force field modification file, to add to the   ***/
/***                 parameters already detected by ReadParmFile() above.  ***/
/***=======================================================================***/
void ReadFrcmodFile(prmset *mp, trajcon *tj)
{
  int nextline;
  cmat cfi;

  /*** Test file existence ***/
  if (tj->fmodfile[0] == '\0') {
    return;
  }
  cfi = Ascii2Mem(tj->fmodfile, 256,
		  "Force field modification file not found.");

  /*** Read bond, angle, dihedral, and improper declarations ***/
  nextline = ScanDeclarations(mp, &cfi, 1, mp->natom, 1, "BLANKLINE");
  nextline = ScanDeclarations(mp, &cfi, 2, mp->nbond, nextline, "BLANKLINE");
  nextline = ScanDeclarations(mp, &cfi, 3, mp->nangl, nextline, "BLANKLINE");
  nextline = ScanDeclarations(mp, &cfi, 4, mp->ntor, nextline, "BLANKLINE");
  nextline = ScanDeclarations(mp, &cfi, 5, mp->ntor, nextline, "BLANKLINE");
  nextline = ScanDeclarations(mp, &cfi, 9, mp->nhb1012, nextline, "BLANKLINE");

  /*** Read modified nobonded parameters ***/
  nextline = ScanDeclarations(mp, &cfi, 8, 0, nextline, "BLANKLINE");
}
