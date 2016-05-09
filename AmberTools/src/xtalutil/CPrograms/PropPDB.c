#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "pdbRead.h"
#include "crdmanip.h"
#include "matrix.h"
#include "myconstants.h"

int main(int argc, char *argv[])
{
  int i;
  int irep[3];
  double boxd[6];
  char tag[64], outfile[256];
  pdb p, t;

  if (argc == 1) {
    printf("\nPropPDB >> A program for propagating a PDB structure.\n\n"
	   "Options:\n"
           "  -p   : the structure to reassemble (PDB format)\n"
           "  -o   : the output structure (PDB format)\n"
	   "  -X   : the length of the first unit cell vector\n"
	   "  -Y   : the length of the second unit cell vector\n"
	   "  -Z   : the length of the third unit cell vector\n"
	   "  -a   : the box alpha angle (degrees)\n"
	   "  -b   : the box beta angle (degrees)\n"
	   "  -g   : the box gamma angle (degrees)\n" 
	   "  -ix  : number of replicas along _X_ vector\n"
	   "  -iy  : number of replicas along _Y_ vector\n"
	   "  -iz  : number of replicas along _Z_ vector\n\n");
    exit(1);
  }

  /*** Input ***/
  p.source[0] = '\0';
  outfile[0] = '\0';
  boxd[0] = -1.0;
  boxd[1] = -1.0;
  boxd[2] = -1.0;
  boxd[3] = 90.0;
  boxd[4] = 90.0;
  boxd[5] = 90.0;
  irep[0] = 1;
  irep[1] = 1;
  irep[2] = 1;
  for (i = 0; i < argc-2; i += 2) {
    strcpy(tag, *++argv);
    if (strcmp(tag, "-p") == 0) {
      strcpy(p.source, *++argv);
    }
    else if (strcmp(tag, "-o") == 0) {
      strcpy(outfile, *++argv);
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
    else if (strcmp(tag, "-a") == 0) {
      boxd[3] = atof(*++argv);
    }
    else if (strcmp(tag, "-b") == 0) {
      boxd[4] = atof(*++argv);
    }
    else if (strcmp(tag, "-g") == 0) {
      boxd[5] = atof(*++argv);
    }
    else if (strcmp(tag, "-ix") == 0) {
      irep[0] = atoi(*++argv);
    }
    else if (strcmp(tag, "-iy") == 0) {
      irep[1] = atoi(*++argv);
    }
    else if (strcmp(tag, "-iz") == 0) {
      irep[2] = atoi(*++argv);
    }
    else {
      printf("PropPDB >> Error.  Unrecognized tag %s.\n", tag);
      exit(1);
    }
  }

  /*** Check ***/
  if (p.source[0] == '\0') {
    printf("PropPDB >> Error.  Original PDB file not specified!\n");
    exit(1);
  }
  if (outfile[0] == '\0') {
    printf("PropPDB >> Error.  Output PDB file not specified!\n");
    exit(1);
  }
  if (boxd[0] < 0.0 || boxd[1] < 0.0 || boxd[2] < 0.0) {
    printf("PropPDB >> Error.  Box dimensions invalid!\n");
    exit(1);
  }

  /*** Transform the box ***/
  for (i = 0; i < 3; i++) {
    boxd[3+i] *= PI/180.0;
  }

  /*** Get the first PDB ***/
  GetPDB(&p, 1, 0);

  /*** Now, make the new PDB ***/
  t = TilePDB(&p, boxd, irep, 1);
  strcpy(t.source, outfile);

  /*** Print PDB ***/
  ModPdbRA(&t);
  PutPDB(&t, t.source, "STANDARD", "HEADER  ", 1);

  return 0;
}
