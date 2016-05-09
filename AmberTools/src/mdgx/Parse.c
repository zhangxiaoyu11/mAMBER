#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Matrix.h"
#include "Constants.h"
#include "Parse.h"
#include "Macros.h"
#include "mdgxVector.h"
#include "ptrajmask.h"

#include "TopologyDS.h"
#include "CrdManipDS.h"
#include "ParamFitDS.h"

/***=======================================================================***/
/*** CountWords: this function counts the number of words appearing on a   ***/
/***             line.  To store all words, sequentially and individually, ***/
/***             use the ParseWords function in this same library.         ***/
/***=======================================================================***/
int CountWords(char* line)
{
  int i, nword, slen, onword, inquotes;

  nword = 0;
  slen = strlen(line);
  onword = 0;
  inquotes = 0;
  for (i = 0; i < slen; i++) {
    if (line[i] == 34 || line[i] == 39) {
      inquotes = 1 - inquotes;
    }
    if (line[i] != ' ' && onword == 0 && line[i] != '\n') {
      onword = 1;
      nword++;
    }
    else if (inquotes == 0 && line[i] == ' ' && onword == 1) {
      onword = 0;
    }
  }

  return nword;
}

/***=======================================================================***/
/*** ToUpper: functino for converting a character to upper case            ***/
/***=======================================================================***/
char ToUpper(char c)
{
  int buffer;

  if (c >= 'a' && c <= 'z') {
    buffer = c - 'a';
    return 'A' + buffer;
  }
  else {
    return c;
  }
}

/***=======================================================================***/
/*** WordIsNumber: function to determine whether a word really looks like  ***/
/***               a number, and then return 0 if false or 1 if true.      ***/
/***=======================================================================***/
int WordIsNumber(char* word)
{
  int i, wlen, hitdot;

  if ((word[0] < '0' || word[0] > '9') && word[0] != '-' && word[0] != '.') {
    return 0;
  }
  wlen = strlen(word);
  hitdot = (word[0] == '.') ? 1 : 0;
  for (i = 1; i < wlen; i++) {
    if ((word[i] < '0' || word[i] > '9') && word[i] != '.') {
      return 0;
    }
    if (word[i] == '.') {
      if (hitdot == 1) {
	return 0;
      }
      else {
	hitdot = 1;
      }
    }
  }

  return 1;
}

/***=======================================================================***/
/*** WordIsAtomType: function to determine whether a word really looks     ***/
/***                 like an atom type, returning 0 if false or 1 if true. ***/
/***=======================================================================***/
int WordIsAtomType(char* word)
{
  int wlen;

  wlen = strlen(word);
  if (wlen > 2 || wlen == 0) {
    return 0;
  }
  if (!(word[0] >= 'A' && word[0] <= 'Z') &&
      !(word[0] >= 'a' && word[0] <= 'z') &&
      !(word[0] >= '0' && word[0] <= '9')) {
    return 0;
  }
  if (!(word[1] >= 'A' && word[1] <= 'Z') && word[1] != '#' &&
      !(word[1] >= 'a' && word[1] <= 'z') && word[1] != '@' &&
      !(word[1] >= '0' && word[1] <= '9') && word[1] != '*' &&
      word[1] != ' ' && word[1] != '\0') {
    return 0;
  }

  return 1;
}

/***=======================================================================***/
/*** ParseWords: this function stores all words appearing on a line.       ***/
/***=======================================================================***/
cmat ParseWords(char* line)
{
  int i, j, onword, inquotes, nc;
  cmat CW;

  CW = CreateCmat(CountWords(line), MAXNAME);
  i = 0;
  j = 0;
  while (i < CW.row) {

    /*** Advance to the next word ***/
    while (line[j] == ' ') {
      j++;
    }

    /*** Read this word ***/
    onword = 1;
    if (line[j] == 34 || line[j] == 39) {
      inquotes = 1;
      j++;
    }
    else {
      inquotes = 0;
    }
    nc = 0;
    while (onword == 1) {

      if (inquotes == 0 && (line[j] == ' ' || line[j] == '\n')) {
	onword = 0;
	CW.map[i][nc] = '\0';
      }
      else if (inquotes == 1 && (line[j] == 34 || line[j] == 39)) {
	onword = 0;
	inquotes = 0;
	CW.map[i][nc] = '\0';
      }
      else {
	CW.map[i][nc] = line[j];
	nc++;
      }
      j++;
    }
    i++;
  }

  return CW;
}

/***=======================================================================***/
/*** RemoveWhiteSpace: function for removing white space at the beginning  ***/
/***                   of a string.                                        ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   a:     the string                                                   ***/
/***   asize: the length of the string                                     ***/
/***=======================================================================***/
void RemoveWhiteSpace(char* a, int asize)
{
  int i, j, num_rem;

  num_rem = 0;
  while (a[num_rem] == ' ' && num_rem < asize) {
    num_rem++;
  }
  j = 0;
  for (i = num_rem; i < asize; i++) {
    a[j] = a[i];
    j++;
  }
  for (i = asize-num_rem; i < asize; i++) {
    a[i] = ' ';
  }
}

/***=======================================================================***/
/*** EqualSpace: add white space before and after any equal symbols in a   ***/
/***             line UNLESS they occur within ('') or ("") quotes.        ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   line:    the line (from an input command file)                      ***/
/***=======================================================================***/
void EqualSpace(char* line)
{
  int i, j, slen, inquotes;
  slen = strlen(line);

  inquotes = 0;
  for (i = 0; i < slen; i++) {
    if (line[i] == '=' && inquotes == 0) {
      line[slen+2] = '\0';
      for (j = slen-1; j > i; j--) {
        line[j+2] = line[j];
      }
      line[i] = ' ';
      line[i+1] = '=';
      line[i+2] = ' ';
      slen += 2;
      i += 2;
    }

    /*** ASCII standard character set used to identify quotation symbols ***/
    if (line[i] == 34 || line[i] == 39) {
      inquotes = 1 - inquotes;
    }
  }
}

/***=======================================================================***/
/*** RemoveComments: this function removes commented text from a line.     ***/
/***                 Comments are denoted by the appearance of '%', '#',   ***/
/***                 or '$' symbols.                                       ***/
/***=======================================================================***/
void RemoveComments(char* line)
{
  int i;
  const int slen = strlen(line);

  for (i = 0; i < slen; i++) {
    if (line[i] == '$' || line[i] == '%' || line[i] == '#') {
      line[i] = '\0';
      break;
    }
  }
}

/***=======================================================================***/
/*** NixCommaCarriage: replaces all commas and carriage returns on a line  ***/
/***                   with white space.                                   ***/
/***=======================================================================***/
void NixCommaCarriage(char* line)
{
  int i, slen, inquotes;

  slen = strlen(line);
  inquotes = 0;
  for (i = 0; i < slen; i++) {
    if (inquotes == 0 && (line[i] == ',' || line[i] == '\n')) {
      line[i] = ' ';
    }

    /*** ASCII standard character set used to identify quotation symbols ***/
    if (line[i] == 34 || line[i] == 39) {
      inquotes = 1 - inquotes;
    }
  }
}

/***=======================================================================***/
/*** AdvanceToSegment: this function searches a command file from the      ***/
/***                   current file pointer (or, from the beginning, if    ***/
/***                   specified) and stops after finding the heading for  ***/
/***                   the section of interest.                            ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   inp:     the input file                                             ***/
/***   segname: the name of the target segment                             ***/
/***   scan0:   if set to 1, causes the file to be rewound before          ***/
/***            beginning the search; if set to 2, causes the segment name ***/
/***            to be searched directly, with no leading '&' symbol, and   ***/
/***            also prompts a file rewind                                 ***/
/***=======================================================================***/
int AdvanceToSegment(FILE *inp, char* segname, int scan0)
{
  int slen, collect;
  char line[MAXLINE];

  /*** Start by rewinding the file ***/
  slen = strlen(segname);
  if (scan0 == 1 || scan0 == 2) {
    rewind(inp);
  }
  collect = 0;
  while (collect == 0) {
    if (fgets(line, MAXLINE, inp) == NULL) {
      break;
    }
    RemoveWhiteSpace(line, MAXLINE);
    if (scan0 < 2) {
      if (line[0] == '&' && strncmp(&line[1], segname, slen) == 0) {
	collect = 1;
      }
    }
    else if (scan0 == 2) {
      if (strncmp(&line[0], segname, slen) == 0) {
	collect = 1;
      }
    }
  }

  return collect;
}

/***=======================================================================***/
/*** DetectNamelistEnd: detect the end of a namelist, signified by a line  ***/
/***                    containing "&end."                                 ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   line:    the line (from an input command file)                      ***/
/***   errmsg:  the error message to display (the name of the function in  ***/
/***              which the call to DetectNamelistEnd originated)          ***/
/***=======================================================================***/
int DetectNamelistEnd(char* line, char* errmsg)
{
  if (line[0] == '&') {
    if (strncmp(&line[1], "end", 3) == 0) {
      return 0;
    }
    else {
      printf("%s >> Error.  New segment encountered before termination of "
             "%s >> current segment.\n", errmsg, errmsg);
      exit(1);
    }
  }
  else if (line[0] == '/' && (line[1] == ' ' || line[1] == '\0')) {
    return 0;
  }

  return 1;
}

/***=======================================================================***/
/*** SeekString: find a tag for a particular string.                       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   L:       the list of words on this line of the input command file   ***/
/***   val:     the value to set if the string flag is found               ***/
/***   sname:   the string flag being sought                               ***/
/***   salias:  alias for the string flag being sought                     ***/
/***=======================================================================***/
void SeekString(cmat L, char* val, char* sname, char* salias)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        strcpy(val, L.map[i+1]);
      }
      else if (i+2 < L.row) {
        strcpy(val, L.map[i+2]);
      }
      else {
        printf("SeekString >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
      }
    }
  }
}

/***=======================================================================***/
/*** SeekStringTripletInc: find a tag for a particular string triplet, and ***/
/***                       increment a counter variable each time the tag  ***/
/***                       is found.                                       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   L:       the list of words on this line of the input command file   ***/
/***   val1:    the first value to set if the string flag is found         ***/
/***   val2:    the second value to set if the string flag is found        ***/
/***   val3:    the second value to set if the string flag is found        ***/
/***   sname:   the string flag being sought                               ***/
/***   salias:  alias for the string flag being sought                     ***/
/***   counter: the counter variable (gets incremented)                    ***/
/***=======================================================================***/
void SeekStringTripletInc(cmat L, char* val1, char* val2, char* val3,
			  char* sname, char* salias, int *counter)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+3 < L.row && strcmp(L.map[i+1], "=") != 0) {
	strcpy(val1, L.map[i+1]);
	strcpy(val2, L.map[i+2]);
	strcpy(val3, L.map[i+3]);
	*counter += 1;
      }
      else if (i+4 < L.row) {
	strcpy(val1, L.map[i+2]);
	strcpy(val2, L.map[i+3]);
	strcpy(val3, L.map[i+4]);
	*counter += 1;
      }
      else {
        printf("SeekString >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
      }
    }
  }
}

/***=======================================================================***/
/*** SeekTorsionID: detect a torsion identifier for optimization.  This    ***/
/***                routine is specific to GetParamNamelist, but placed    ***/
/***                here to keep Command.c tidy.                           ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   L:       the list of words on this line of the input command file   ***/
/***   xt:      the next adjustable torsion term to specify                ***/
/***   sname:   the string flag being sought                               ***/
/***   salias:  alias for the string flag being sought                     ***/
/***   counter: the counter variable (gets incremented, records the number ***/
/***            of adjustable torsions found thus far)                     ***/
/***=======================================================================***/
void SeekTorsionID(cmat L, prmset *mp, char* sname, char* salias,
		   int *maxhadj)
{
  int i, j, alen, blen, clen, dlen;
  torterm *xt;

  for (i = 0; i < L.row; i++) {
    xt = &mp->hadj[mp->nhadj];
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i + 4 < L.row && strcmp(L.map[i+1], "=") != 0) {
        strcpy(xt->atype, L.map[i+1]);
        strcpy(xt->btype, L.map[i+2]);
        strcpy(xt->ctype, L.map[i+3]);
        strcpy(xt->dtype, L.map[i+4]);
	xt->K = -1.0;
        mp->nhadj += 1;
      }
      else if (i+5 < L.row) {
        strcpy(xt->atype, L.map[i+2]);
        strcpy(xt->btype, L.map[i+3]);
        strcpy(xt->ctype, L.map[i+4]);
        strcpy(xt->dtype, L.map[i+5]);
	xt->K = -1.0;
        mp->nhadj += 1;
      }
      else {
        printf("SeekString >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
      }

      /*** Touch up the atom names ***/
      alen = strlen(xt->atype);
      blen = strlen(xt->btype);
      clen = strlen(xt->ctype);
      dlen = strlen(xt->dtype);
      for (j = 0; j < 4; j++) {
	if (j >= alen) xt->atype[j] = ' ';
	if (j >= blen) xt->btype[j] = ' ';
	if (j >= clen) xt->ctype[j] = ' ';
	if (j >= dlen) xt->dtype[j] = ' ';
      }

      /*** Extend the adjustable torsions array if needed ***/
      if (mp->nhadj == *maxhadj) {
	*maxhadj += 32;
	mp->hadj = (torterm*)realloc(mp->hadj, (*maxhadj)*sizeof(torterm));
      }
    }
  }
}

/***=======================================================================***/
/*** ReadNumericalSuffix: read the numerical suffix to a control variable. ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   term:     the control variable                                      ***/
/***   nstart:   the anticipated start of the numerical part of term       ***/
/***   notanum:  flag to indicate that the numerical part contained        ***/
/***             non-integer characters                                    ***/
/***=======================================================================***/
int ReadNumericalSuffix(char* term, int nstart, int *notanum)
{
  int i, vpos, lflag;

  /*** Determine any numerical suffix to the flag ***/
  lflag = strlen(term);
  *notanum = 0;
  if (lflag == nstart) {
    *notanum = 1;
    return -1;
  }
  for (i = nstart; i < lflag; i++) {
    if (term[i] < '0' || term[i] > '9') {
      *notanum = 1;
      return -1;
    }
  }
  vpos = atoi(&term[nstart])-1;

  return vpos;
}

/***=======================================================================***/
/*** SeekNString: find a tag for a particular string and then check it for ***/
/***              numerical extensions denoting that the string denotes an ***/
/***              input for an alternate molecular system.                 ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   Similar to SeekString above, but val is now a cmat struct with more ***/
/***   than one row.  The integer array fspec denotes whether any rows of  ***/
/***   val had already been specified before this function was called, and ***/
/***   therefore should not be reassigned.                                 ***/
/***=======================================================================***/
cmat SeekNString(cmat L, cmat *val, int* fspec, char* sname, char* salias)
{
  int i, nstart, notanumber, vpos, sfound, afound;

  const int lsname = strlen(sname);
  const int lsalias = strlen(salias);
  for (i = 0; i < L.row; i++) {
    sfound = 0;
    afound = 0;
    if (strncmp(L.map[i], sname, lsname) == 0) {
      sfound = 1;
      nstart = lsname;
    }
    else if (strncmp(L.map[i], salias, lsalias) == 0) {
      afound = 1;
      nstart = lsalias;
    }
    if (sfound == 1 || afound == 1) {
      vpos = ReadNumericalSuffix(L.map[i], nstart, &notanumber);
      if (notanumber == 1) {
	continue;
      }
      if (vpos < 0 || vpos > MAXSYS) {
        printf("SeekNString >> Error.  Invalid numerical value %d specified "
	       "for identifier\nSeekNString >> %s/%s.\n", vpos+1, sname,
	       salias);
	exit(1);
      }
      if (vpos >= val->row) {
	*val = ReallocCmat(val, vpos+1, MAXNAME);
      }
      else if (fspec[vpos] == 1) {
	continue;
      }

      /*** This string value may be set ***/
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        strcpy(val->map[vpos], L.map[i+1]);
      }
      else if (i+2 < L.row) {
        strcpy(val->map[vpos], L.map[i+2]);
      }
      else {
        printf("SeekNString >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
	exit(1);
      }
    }
  }

  return *val;
}

/***=======================================================================***/
/*** SeekReal: find a tag for a particular real value.                     ***/
/***=======================================================================***/
void SeekReal(cmat L, double *val, char* sname, char* salias)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        *val = atof(L.map[i+1]);
      }
      else if (i+2 < L.row) {
        *val = atof(L.map[i+2]);
      }
      else {
        printf("SeekReal >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
	exit(1);
      }
    }
  }
}

/***=======================================================================***/
/*** SeekNReal: find a tag for a particular real value, one of a series.   ***/
/***            Unlike SeekNString, this function does not check for       ***/
/***            the pre-specification of values and does not limit the     ***/
/***            series number to within [ 0, MAXSYS ].                     ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   Same as SeekReal above, except in this case val is an array not a   ***/
/***   pointer.  The val array must be pre-allocated.  maxidx indicates    ***/
/***   the maximum permissible value ID.                                   ***/
/***=======================================================================***/
void SeekNReal(cmat L, double* val, char* sname, char* salias, int maxidx)
{
  int i, lsname, lsalias, vidx, match, nstart, notanumber;

  lsname = strlen(sname);
  lsalias = strlen(salias);
  for (i = 0; i < L.row; i++) {
    match = 0;
    if (strncmp(L.map[i], sname, lsname) == 0) {
      match = 1;
      nstart = lsname;
    }
    else if (strncmp(L.map[i], salias, lsalias) == 0) {
      match = 1;
      nstart = lsalias;
    }
    if (match == 1) {
      vidx = ReadNumericalSuffix(L.map[i], nstart, &notanumber);
      if (notanumber == 1) {
	continue;
      }
      if (vidx < 0 || vidx >= maxidx) {
	printf("SeekNReal >> Error.  Invalid numerical value %d specified "
               "for identifier\nSeekNReal >> %s/%s.\n", vidx+1, sname,
               salias);
	exit(1);
      }
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        val[vidx] = atof(L.map[i+1]);
      }
      else if (i+2 < L.row) {
        val[vidx] = atof(L.map[i+2]);
      }
      else {
        printf("SeekNReal >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
        exit(1);
      }
    }
  }
}

/***=======================================================================***/
/*** SeekInt: find a tag for a particular integer value.                   ***/
/***=======================================================================***/
void SeekInt(cmat L, int *val, char* sname, char* salias)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        *val = atoi(L.map[i+1]);
      }
      else if (i+2 < L.row) {
        *val = atoi(L.map[i+2]);
      }
      else {
        printf("SeekInt >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
	exit(1);
      }
    }
  }
}

/***=======================================================================***/
/*** SeekLLInt: find a tag for a particular long long integer value.       ***/
/***=======================================================================***/
void SeekLLInt(cmat L, long long int *val, char* sname, char* salias)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        sscanf(L.map[i+1], "%lld", val);
      }
      else if (i+2 < L.row) {
        sscanf(L.map[i+2], "%lld", val);
      }
      else {
        printf("SeekLLInt >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
	exit(1);
      }
    }
  }
}

/***=======================================================================***/
/*** ParseAmbmask: parse an ambmask string.  Uses ptrajmask.c              ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   maskstr:   the mask string                                          ***/
/***   tp:        the topology                                             ***/
/***   crd:       the coordinates                                          ***/
/***=======================================================================***/
int* ParseAmbMask(char* maskstr, prmtop *tp, coord *crd)
{
  int i, j;
  int* imask;
  char* cmask;
  NAME* atomNames;
  NAME* residueNames;
  NAME* atomTypes;

  /*** Reshape the atom and residue name arrays ***/
  atomNames = (NAME*)malloc(tp->natom*sizeof(NAME));
  atomTypes = (NAME*)malloc(tp->natom*sizeof(NAME));
  residueNames = (NAME*)malloc(tp->nres*sizeof(NAME));
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < 4; j++) {
      atomNames[i][j] = tp->AtomNames[4*i+j];
      atomTypes[i][j] = tp->AtomTypes[4*i+j];
    }
    atomNames[i][4] = '\0';
    atomTypes[i][4] = '\0';
  }
  for (i = 0; i < tp->nres; i++) {
    for (j = 0; j < 4; j++) {
      residueNames[i][j] = tp->ResNames[4*i+j];
    }
    residueNames[i][4] = '\0';
  }

  /*** Call to Dan Roe's C-implementation parser ***/
  cmask = parseMaskString(maskstr, tp->natom, tp->nres, atomNames,
			  residueNames, tp->ResLims, crd->loc, atomTypes, 0);

  /*** Conversion to mdgx mask format ***/
  imask = (int*)malloc(tp->natom*sizeof(int));
  for (i = 0; i < tp->natom; i++) {
    imask[i] = (cmask[i] == 'T') ? 1 : 0;
  }

  /*** Free allocated memory ***/
  free(cmask);
  free(atomNames);
  free(residueNames);
  free(atomTypes);

  return imask;
}

/***=======================================================================***/
/*** FOpenSafe: open a new file, if and only if it does not already exist  ***/
/***            if file overwriting is not permitted.                      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   fname:  the name of the file to open                                ***/
/***   ovrwrt: flag to authorize overwriting (1 permits, 0 restricts)      ***/
/***=======================================================================***/
FILE* FOpenSafe(char* fname, int ovrwrt)
{
  FILE* outp;

  if (ovrwrt == 1 || (ovrwrt == 0 && (outp = fopen(fname, "r")) == NULL)) {
    outp = fopen(fname, "w");
  }
  else {
    printf("FOpenSafe >> Error.  File %s already exists.\n", fname);
    exit(1);
  }

  return outp;
}

/***=======================================================================***/
/*** ReadNumericalShorthand: converts a string, which may contain terms    ***/
/***                         such as "GB", "MB", "gb", or "KB" into a long ***/
/***                         long int.  "GB" or "G" or "gb" or "g" equals  ***/
/***                         gigabyte, "MB" megabtye, etc.  Numbers parsed ***/
/***                         by this routine cannot be negative.           ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   numstr:      the string that will be converted into a number        ***/
/***=======================================================================***/
long long int ReadNumericalShorthand(char* numstr)
{
  int i, slen;
  long long int product, tensplace, multiplier;

  /*** Parse the string for special characters ***/
  slen = strlen(numstr);
  multiplier = 1;
  for (i = 0; i < slen; i++) {
    if (numstr[i] < '0' || numstr[i] > '9') {
      if (numstr[i] == 'G' || numstr[i] == 'g') {
	multiplier = 1073741824;
      }
      else if (numstr[i] == 'M' || numstr[i] == 'm') {
	multiplier = 1048576;
      }
      else if (numstr[i] == 'K' || numstr[i] == 'k') {
	multiplier = 1024;
      }
      else {
	printf("ReadNummericalShorthand >> Error.  Unable to parse %s into "
	       "digits.\n", numstr);
	exit(1);
      }
      if (i < slen-2 ||
	  (i < slen-1 && !(numstr[i+1] == 'B' || numstr[i+1] == 'b'))) {
	printf("ReadNummericalShorthand >> Error.  Unable to parse %s into "
	       "digits.\n", numstr);
	exit(1);
      }
      slen = i;
    }
  }

  /*** Compose the number ***/
  tensplace = 1;
  product = 0;
  for (i = slen-1; i >= 0; i--) {
    product += (numstr[i] - '0')*tensplace;
    tensplace *= 10;
  }
  product *= multiplier;

  return product;
}
