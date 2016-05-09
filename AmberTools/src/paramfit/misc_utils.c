/*****************************************************
 * AMBER Bond Angle and Dihedral Parameter Optimiser *
 *                                                   *
 *           Written by: Robin Betz  (2011)          *
 *                       Ross Walker (2004)          *
 *                   UC San Diego                    *
 *           San Diego Supercomputer Center          *
 *            La Jolla, California, 92092            *
 *                       USA                         *
 *****************************************************/

/*misc_utils.c*/

/*This module contains a number of small utilities used in the program*/
/* Proceedures included are:
   void print_close_line_box(int no_spaces) Simple routine to print a set number of spaces and then | - used
                  for drawing boxes etc
   void print_open_line_box(int *i)
   int check_for_valid_filename(const char *data_string, const int length)   Checks if data_string contains any
             illegal characters for filenames   
                 
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"

#include "function_def.h"

void print_close_line_box(int no_spaces)
{
  register int i;
/*This routine prints the closing of a line in a box of style*/
/* ---------------- */
/* |              | */
/* ---------------- */
  for (i=0;i<no_spaces;i++)
  {
    printf(" ");
  }
  printf("|\n");
}

void print_open_line_box(int *i)
{
  /*used with print_close_line_box to print the beginning of the ascii box*/
  *i=printf("        | ");
}

int check_for_valid_filename(const char *data_string, const int length)
{
  /*This won't catch everything but is useful*/
  int temp_count;
  for (temp_count=0;temp_count<length;temp_count++)
  {
    if (
       (*(data_string+temp_count)<32)
       ||
       (((*(data_string+temp_count))>33)&&((*(data_string+temp_count))<43))
       ||
       (((*(data_string+temp_count))>57)&&((*(data_string+temp_count))<65))
       ||
       (*(data_string+temp_count)==92)
       ||
       (*(data_string+temp_count)==124)
       ||
       (*(data_string+temp_count)>126)
       )
       return INVALID_FORMAT; /*error invalid character for filename*/
  }
  return SUCCESS;
}

int s_getline( char *line, int max, FILE *fp )
{
  /*
   *   "safe" version of getline: exits with error if EOF is found
   */

  if( fgets( line, max, fp ) == NULL )
  {
    return FILE_READ_FAIL;
  }
  else
  {
    return strlen( line );    
  }
}

int find_flag( FILE *fptr, char *label )
{
   /*
    *  Used to locate a title in a prmtop - leaves the file pointer
    *  positioned at the line below the selected label
    *  If the label doesn't exist it returns INVALID_LINE
    *  NOTE: rewinds the file pointer before searching
    */

    int retval;
    int i;
    char *buffer;

    buffer = (char *)calloc(BUFFER_SIZE,sizeof(char));
    if (buffer==NULL)
    {
      malloc_failure_char("find_flag", "buffer", BUFFER_SIZE);
      return ALLOC_FAIL;
    }
    
    rewind(fptr);
    while (TRUE) /*loop forever since this will be broken by an internal return statement for either success or an EOF*/
    {
       while (strncmp(buffer, "%FLAG", 5) != 0) /*loop through until we hit a flag*/
       {
         retval = s_getline( buffer, BUFFER_SIZE, fptr );
         /*s_getline returns FILE_READ_FAIL on EOF*/
         if ( retval==FILE_READ_FAIL )
         {
           free(buffer);
          buffer = NULL;
           return INVALID_LINE;
         }
       }

       /*found the flag, now check if it is the correct label*/
       if ( !strncmp( buffer+6, label, strlen(label) ) )
       {
         /*We found the label we were looking for*/
         free(buffer);
	 buffer = NULL;
         return SUCCESS;
       }
       /*if we got to here then it wasn't the correct flag
         clear the beginning of the buffer so the next strncmp doesn't trigger on the current contents
       */
       for (i=0;i<=11;++i)
          buffer[i]=' ';
    }   
}


/*Designed to read 4 characters from a file stream and null terminate
NAME_SIZE is defined in prmtop_params.h*/
int name_copy( FILE *fptr, char *stringp)
{
  int i;
  int ierror;
  for (i = 0; i < NAME_SIZE - 1; ++i)
  {
    do
    {
      ierror = fscanf(fptr, "%c", &stringp[i]);
      if ( ierror == EOF ) /*Hit end of file*/
        return INVALID_LINE;
      if (stringp[i] == ' ' && i!=1) stringp[i] = (char) 0; // no spaces allowed in atom name after 1 char

    } while ( stringp[i] == '\n');

    if ( ierror == 0 ) /*Got a carriage return as part of the atom name*/
        return INVALID_LINE;
  }
  stringp[NAME_SIZE-1] = (char) 0;
  return SUCCESS;
}

/*unObfuscateAtom( int )
   ...inside the topology file, PARM for some unknown reason
   obfuscates each atom number by (when positive), at = ( at + 3 ) / 3
   This routine, when passed an atom number as read from the parm
   file converts it to the true number., e.g number 24 would actually be
   atom number 9.*/
int unObfuscateAtom( int at )
{
  if ( at < 0 )
    at = ( -at + 3 ) / 3;
  else
    at = ( at + 3 ) / 3;
  return( at );
}

/* Re-obfuscate your atom for writing a prmtop */
int ObfuscateAtom( int at )
{
  at = (3*at) - 3;
  return at;
}

/*Calculates the distance between 2 points in 3 dimensions*/
double calc_bond_length(double bond1x, double bond1y, double bond1z, double bond2x, double bond2y, double bond2z)
{

  return(sqrt((bond1x-bond2x)*(bond1x-bond2x)+(bond1y-bond2y)*(bond1y-bond2y)+(bond1z-bond2z)*(bond1z-bond2z)));

}


/*Calculates the angle between 3 points in 3 dimensions*/
double calc_angle_radians(double atom1x, double atom1y, double atom1z, double atom2x, double atom2y, double atom2z,
                          double atom3x, double atom3y, double atom3z)
{
  /*calculates the angle in radians between 3 points in 3 dimensions. Angle is defined as:

    Vector1 = atom1-atom2
    Vector2 = atom3-atom2

                    Vector1 . Vector2
    Angle = cos^-1( ------------------)
                    |Vector1||Vector2|

    Where |Vector1| = sqrt(vector1x^2+vector1y^2+vector1z^2)
  */

  double temp_angle;
  double vector1[3];
  double vector2[3];
  double modvector1;
  double modvector2;
  double vectordot;

  vector1[0]=atom1x-atom2x;
  vector1[1]=atom1y-atom2y;
  vector1[2]=atom1z-atom2z;
  vector2[0]=atom3x-atom2x;
  vector2[1]=atom3y-atom2y;
  vector2[2]=atom3z-atom2z;

  modvector1=sqrt((vector1[0]*vector1[0])+(vector1[1]*vector1[1])+(vector1[2]*vector1[2]));
  modvector2=sqrt((vector2[0]*vector2[0])+(vector2[1]*vector2[1])+(vector2[2]*vector2[2]));

  vectordot=(vector1[0]*vector2[0])+(vector1[1]*vector2[1])+(vector1[2]*vector2[2]);

  temp_angle=acos(vectordot/(modvector1*modvector2));
  return temp_angle;

}

/*calculates the dihedral between 4 points in 3 dimensions*/
double calc_dihedral_radians(double atom1x, double atom1y, double atom1z, double atom2x, double atom2y, double atom2z,
                          double atom3x, double atom3y, double atom3z, double atom4x, double atom4y, double atom4z)
{

 /*Given four atoms A,B,C & D:

 A      D
  \    /
   B--C

 the torsion angle along the torsion axis B--C is the angle between the planes ABC and BCD. The best way to calculate this is
 in terms of 3 vectors:

 a = A->B
 b = B->C
 c = C->D

 In terms of these vectors the dihedral angle is then:

 Dihedral = cos^-1[(a X b) * (b X c)]
                   -----------------
                   (|a X b|*|b X c|)

 The sign of the dihedral is then found from the triple scalar product calculated by evaluating the
 determinant of the matrix:

   | a(x) a(y) a(z) |    | a(x) a(y) a(z) |
   | b(x) b(y) b(z) | -> | b(x) b(y) b(z) |
   | c(x) c(y) c(z) |    | c(x) c(y) c(z) |

   a[x] * (b[y]*c[z]-c[y]*b[z]) - a[y]*(b[x]*c[z]-c[x]*b[z]) + a[z]*(b[x]*c[y]-c[x]*b[y])

 */

 double vector1x, vector1y, vector1z;
 double vector2x, vector2y, vector2z;
 double vector3x, vector3y, vector3z;
 double v1Xv2xv2Xv3;
 double v1Xv2;
 double v2Xv3;
 double cosdihe;
 double dihedral;
 double sign;
 
 vector1x = atom2x - atom1x;
 vector1y = atom2y - atom1y;
 vector1z = atom2z - atom1z;

 vector2x = atom2x - atom3x;
 vector2y = atom2y - atom3y;
 vector2z = atom2z - atom3z;

 vector3x = atom3x - atom4x;
 vector3y = atom3y - atom4y;
 vector3z = atom3z - atom4z;

 v1Xv2xv2Xv3 = (vector1y * vector2z - vector1z * vector2y) * (vector2y * vector3z - vector2z * vector3y) +
               (vector2x * vector1z - vector1x * vector2z) * (vector3x * vector2z - vector2x * vector3z) +
               (vector1x * vector2y - vector2x * vector1y) * (vector2x * vector3y - vector3x * vector2y);
               


 v1Xv2 = (vector1y * vector2z - vector1z * vector2y) * (vector1y * vector2z - vector1z * vector2y) +
         (vector2x * vector1z - vector1x * vector2z) * (vector2x * vector1z - vector1x * vector2z) +
         (vector1x * vector2y - vector2x * vector1y) * (vector1x * vector2y - vector2x * vector1y);

 v2Xv3 = (vector2y * vector3z - vector2z * vector3y) * (vector2y * vector3z - vector2z * vector3y ) +
         (vector3x * vector2z - vector2x * vector3z) * (vector3x * vector2z - vector2x * vector3z ) +
         (vector2x * vector3y - vector3x * vector2y) * (vector2x * vector3y - vector3x * vector2y );
         
 cosdihe = v1Xv2xv2Xv3 / sqrt (v1Xv2 * v2Xv3);
 if (cosdihe > 1.0) // take care of some floating point issues that will sometimes cause an out of bounds to the arccos
   cosdihe = 1.0;
 if (cosdihe < -1.0)
   cosdihe = -1.0;
 dihedral = acos(cosdihe);

 sign = vector1x * (vector2y * vector3z - vector3y * vector2z) - vector1y * (vector2x * vector3z - vector3x * vector2z)
        + vector1z * (vector2x * vector3y - vector3x * vector2y);
        
 if (sign>0)
 {
    if (dihedral>=0) return dihedral;
    if (dihedral<0)  return -dihedral;
 }
 else if (sign<0)
 {
    if (dihedral>0) return -dihedral;
    if (dihedral<=0) return dihedral;
 }

 return dihedral; /*Should never actually gets to this line due to above but this suppresses compiler warnings about
                    hitting the end of a non-void function*/
}

void calc_fit_dimensions(global_options_struct *global_options, parm_struct *parm_data)
{
  int i;
  /*Calculates the number of dimensions of the fit*/
  global_options->NDIMENSIONS=0;
  if (global_options->K_FIT==YES)
    ++global_options->NDIMENSIONS;

  for (i=0;i<parm_data->unique_bonds_found;++i)
  {
    if (parm_data->bond_data[i].DO_FIT_KR==YES)
      ++global_options->NDIMENSIONS;
    if (parm_data->bond_data[i].DO_FIT_REQ==YES)    
      ++global_options->NDIMENSIONS;
  }
  for (i=0;i<parm_data->unique_angles_found;++i)
  {
    if (parm_data->angle_data[i].DO_FIT_KT==YES)
      ++global_options->NDIMENSIONS;
    if (parm_data->angle_data[i].DO_FIT_THEQ==YES)
      ++global_options->NDIMENSIONS;
  }
  for (i=0;i<parm_data->unique_dihedrals_found;++i)
  {
    if (parm_data->dihedral_data[i].DO_FIT_KP==YES)
      ++global_options->NDIMENSIONS;
    if (parm_data->dihedral_data[i].DO_FIT_NP==YES)
      ++global_options->NDIMENSIONS;
    if (parm_data->dihedral_data[i].DO_FIT_PHASE==YES)
      ++global_options->NDIMENSIONS;
  }
}

int modify_params_scratch_data(global_options_struct *global_options, parm_struct *parm_data, double *parameters, readwrite_t MODE)
{
  /*This routine extracts parameters marked as variable and puts them in the linear double array *parameters
   *parameters must have been pre-allocated to NDIMENSIONS long*/

  /*Note, it does this in the very strict order:

    K
    BOND x (KR, KEQ)
    BOND y (KR, KEQ)
    .
    ANGLE x (KT, THEQ)
    ANGLE y (KT, THEQ)
    .
    DIHEDRAL x (KP, PN, PHASE)
    DIHEDRAL y (KP, PN, PHASE)
    .
    To put the data back into the parmtop array you must exactly repeat this procedure

    if MODE=READ - we copy from parm_data to parameters
    if MODE=WRITE - we copy from parameters to parm_data

    Returns the number of variable parameters successfully extracted - in the interests of speed it does not
    check if number_extracted exceeds NDIMENSIONS so will probably seg fault if NDIMENSIONS and thus *parameters
    is too small.

    I am aware that this whole procedure here is clunky and slow but it makes it significantly easier to understand
    and debug. At some point I will replace this whole system with a much more efficient method
    */

  int number_extracted;
  int i;

  number_extracted=0;

  if (MODE==READ)
  {
    if (global_options->K_FIT==YES)
    {
      parameters[number_extracted]=global_options->K;
      ++number_extracted;
    }
    /*Now do all the bonds*/
    for (i=0;i<parm_data->unique_bonds_found;++i)
    {
      /*For each bond in turn check if KR and/or REQ are variable*/
      if (parm_data->bond_data[i].DO_FIT_KR==YES)
      {
         parameters[number_extracted]=parm_data->bond_data[i].rk;
         ++number_extracted;
      }
      if (parm_data->bond_data[i].DO_FIT_REQ==YES)
      {
         parameters[number_extracted]=parm_data->bond_data[i].req;
         ++number_extracted;
      }
    }
    /*Now do the angles*/
    for (i=0;i<parm_data->unique_angles_found;++i)
    {
      /*For each angle in turn check if KT and/or THEQ are variable*/
      if (parm_data->angle_data[i].DO_FIT_KT==YES)
      {
         parameters[number_extracted]=parm_data->angle_data[i].tk;
         ++number_extracted;
      }
      if (parm_data->angle_data[i].DO_FIT_THEQ==YES)
      {
         parameters[number_extracted]=parm_data->angle_data[i].teq;
         ++number_extracted;
      }
    }
    /*now do the dihedrals*/
    for (i=0;i<parm_data->unique_dihedrals_found;++i)
    {
      /*For each diheedral in turn check if KP and/or NP and/or PHASE are variable*/
      if (parm_data->dihedral_data[i].DO_FIT_KP==YES)
      {
         parameters[number_extracted]=parm_data->dihedral_data[i].pk;
         ++number_extracted;
      }
      if (parm_data->dihedral_data[i].DO_FIT_NP==YES)
      {
         parameters[number_extracted]=parm_data->dihedral_data[i].pn;
         ++number_extracted;
      }
      if (parm_data->dihedral_data[i].DO_FIT_PHASE==YES)
      {
         parameters[number_extracted]=parm_data->dihedral_data[i].phase;
         ++number_extracted;
      }
    }
  }
  else if (MODE==WRITE)
  {
    /* Write the variables into the parm structure */
    
    if (global_options->K_FIT==YES)
    {
      global_options->K=parameters[number_extracted];
      ++number_extracted;
    }
    
    /*Now do all the bonds*/
    for (i=0;i<parm_data->unique_bonds_found;++i)
    {
      /* KR stays between 100 and 1000 */
      if (parm_data->bond_data[i].DO_FIT_KR==YES)
      {
        parm_data->bond_data[i].rk=parameters[number_extracted];
        ++number_extracted;
      }
      /* Length between 0 and 3 Angstroms */
      if (parm_data->bond_data[i].DO_FIT_REQ==YES)
      {
        parm_data->bond_data[i].req=parameters[number_extracted];
        ++number_extracted;
      }
    }
    /*Now do the angles*/
    for (i=0;i<parm_data->unique_angles_found;++i)
    {
      /*KT stays between 0 and 170 */
      if (parm_data->angle_data[i].DO_FIT_KT==YES)
      {
        parm_data->angle_data[i].tk=parameters[number_extracted];
        ++number_extracted;
      }
      /* Make sure angle periodicity is okay */
      if (parm_data->angle_data[i].DO_FIT_THEQ==YES)
      {
        parm_data->angle_data[i].teq=parameters[number_extracted];
        ++number_extracted;
      }
    }
    /*now do the dihedrals*/
    for (i=0;i<parm_data->unique_dihedrals_found;++i)
    {
      /*For each diheedral in turn check if KP and/or NP and/or PHASE are variable*/
      // PK stays between -30 and 30 (Glycam ff lists some negative kp, for example)
      if (parm_data->dihedral_data[i].DO_FIT_KP==YES)
      {
        parm_data->dihedral_data[i].pk=parameters[number_extracted];
        ++number_extracted;
      }
      // PN positive or negative integer between 1 and 5, so round it accordingly
      if (parm_data->dihedral_data[i].DO_FIT_NP==YES)
      {
        parm_data->dihedral_data[i].pn = parameters[number_extracted];
        ++number_extracted;
      }
      // Make sure phase stays in a valid angle range and we go around the unit circle if necessary
      if (parm_data->dihedral_data[i].DO_FIT_PHASE==YES)
      {
        parm_data->dihedral_data[i].phase=parameters[number_extracted];
        ++number_extracted;
      }
    }
  }
  else
  {
    printf("   ERROR IN modify_params_scratch_data() - MODE = %d - UNKNOWN MODE\n",MODE);
    printf("          NON-FATAL, ATTEMPTING TO CONTINUE.\n");
    fflush(stdout); /*Flush the printf buffer*/
  }

  return(number_extracted);
}

void print_backtrace(int signal) // this is for debugging and will print a backtrace
{
#ifdef CYGWIN
  fprintf(stderr, "Error: signal %d\n",signal);
#else
  void *array[10];
  size_t size;
  size = backtrace(array, 10);
  
  fprintf(stderr, "Error: signal %d:\n", signal);
  backtrace_symbols_fd(array,size,2);
#endif
  exit(1);
}

void handle_sigint(int param)
{
  /* to be implemented - ideally this will print out the results of the 
  algorithm (esp. the GA) so far if it's taking a really long time and you just
  want to kill it to see how far it's gotten. Requires global_options to be
  actually global though since it can't have any parameters */
  printf("\n!   ERROR: Program terminated.\n");
  exit(ABORT);
}

int calculate_no_fit_params(parm_struct *parm_data, short int MODE)
{
  /*Calculates the number of parameters to be fit corresponding to MODE*/
  int number_params;
  int i;
  number_params=0;

  if (MODE==BONDS)
  {
    for (i=0;i<parm_data->unique_bonds_found;++i)
    {
      /*For each bond in turn check if KR and/or REQ are variable*/
      if (parm_data->bond_data[i].DO_FIT_KR==YES)
         ++number_params;
      if (parm_data->bond_data[i].DO_FIT_REQ==YES)
         ++number_params;
    }
  }
  else if (MODE==ANGLES)
  {
    for (i=0;i<parm_data->unique_angles_found;++i)
    {
      /*For each bond in turn check if KR and/or REQ are variable*/
      if (parm_data->angle_data[i].DO_FIT_KT==YES)
         ++number_params;
      if (parm_data->angle_data[i].DO_FIT_THEQ==YES)
         ++number_params;
    }
  }
  else if (MODE==DIHEDRALS)
  {
    for (i=0;i<parm_data->unique_dihedrals_found;++i)
    {
      /*For each bond in turn check if KR and/or REQ are variable*/
      if (parm_data->dihedral_data[i].DO_FIT_KP==YES)
         ++number_params;
      if (parm_data->dihedral_data[i].DO_FIT_NP==YES)
         ++number_params;
      if (parm_data->dihedral_data[i].DO_FIT_PHASE==YES)
         ++number_params;
    }
  }
  else
  {
    printf("   ERROR IN calculate_no_fit_params() - RW = %d - UNKNOWN MODE\n",MODE);
    printf("             NON-FATAL, ATTEMPTING TO CONTINUE.\n");
    fflush(stdout); /*Flush the printf buffer*/
  }

  return(number_params);
}

int dihedral_types_equal(dihedral_data_struct *first, dihedral_data_struct *second)
{
  if (  !strcmp(first->atom_type1, second->atom_type1) &&
        !strcmp(first->atom_type2, second->atom_type2) &&
        !strcmp(first->atom_type3, second->atom_type3) &&
        !strcmp(first->atom_type4, second->atom_type4) )
    return YES;
  else
    return NO;
}

void print_dihedral(dihedral_data_struct *hi)
{
  printf("%s-%s-%s-%s: KP: %.3f NP: %.3f Phase: %.3f\n", hi->atom_type1, hi->atom_type2, hi->atom_type3, hi->atom_type4, hi->pk, hi->pn, hi->phase*RADIAN_TO_DEGREE);
}

int compare_energy(const void *a, const void *b) // compare function for qsort for coords structures
{
  coords_struct *ia = (coords_struct*)a;
  coords_struct *ib = (coords_struct*)b;
  
  if (ia->energy > ib->energy)
    return 1;
  else if (ia->energy < ib->energy)
    return -1;
  else
    return 0;
}

// This will add more dihedral terms to the existing ones so you can
// (potentially) get a better fit.
int not_enough_dihedrals(parm_struct *parm_data, int n)
{
  // re-alloc the dihedral data
  dihedral_data_struct *more = (dihedral_data_struct*)malloc(sizeof(dihedral_data_struct)*parm_data->unique_dihedrals_found);
  if (more==NULL)
  {
    printf("   ERROR: Failed to reallocate dihedral structures.\n");
    exit(ALLOC_FAIL);
  }
  
  // copy the data into the larger structure
  int i, j, k;
  short int new_copy;
  int counter=0;
  int n_dihedrals = parm_data->unique_dihedrals_found; // since parm_data's will change as more dihedral terms are added
  
  for (i=0; i<n_dihedrals; ++i)
  {
    // only add to dihedrals that will be fit that don't have enough terms already
    if (parm_data->dihedral_data[i].num_terms < n && 
      (parm_data->dihedral_data[i].DO_FIT_KP || parm_data->dihedral_data[i].DO_FIT_NP || parm_data->dihedral_data[i].DO_FIT_PHASE) )
    {
      // Realloc the temporary dihedral data structure for room for n-num more
      parm_data->unique_dihedrals_found += (n-parm_data->dihedral_data[i].num_terms);
      more = (dihedral_data_struct*)realloc(more, sizeof(dihedral_data_struct)*parm_data->unique_dihedrals_found);
      new_copy = YES;
    }
    else // just copy this one over
      new_copy = NO;
    
    for(j=0; j<( new_copy==YES ? n-parm_data->dihedral_data[i].num_terms+1 : 1); ++j)
    {
      // copy basic info about this dihedral
      strcpy(more[counter].atom_type1, parm_data->dihedral_data[i].atom_type1);
      strcpy(more[counter].atom_type2, parm_data->dihedral_data[i].atom_type2);
      strcpy(more[counter].atom_type3, parm_data->dihedral_data[i].atom_type3);
      strcpy(more[counter].atom_type4, parm_data->dihedral_data[i].atom_type4);
      more[counter].improper = parm_data->dihedral_data[i].improper;
      
      // copy info about which parameters will be fit
      more[counter].phase = parm_data->dihedral_data[i].phase;
      more[counter].pk = parm_data->dihedral_data[i].pk;
      more[counter].DO_FIT_PHASE = parm_data->dihedral_data[i].DO_FIT_PHASE;
      more[counter].DO_FIT_KP = parm_data->dihedral_data[i].DO_FIT_KP;
      if (new_copy==YES)
      {
        // give it a random periodicity
        more[counter].pn = (rand()&1) ? (double)(rand()%6)+1.0 : -1.0*((double)(rand()%6)+1.0);
        more[counter].DO_FIT_NP=YES;
      }
      else
      {
        more[counter].pn = parm_data->dihedral_data[i].pn;
        more[counter].DO_FIT_NP = parm_data->dihedral_data[i].DO_FIT_NP;
      }
      
      // copy info about atoms in this type of dihedral
      more[counter].number = parm_data->dihedral_data[i].number;
      for (k=0; k<MAX_DIHEDRALS_PER_TYPE; ++k)
      {
        more[counter].atom1[k] = parm_data->dihedral_data[i].atom1[k];
        more[counter].atom2[k] = parm_data->dihedral_data[i].atom2[k];
        more[counter].atom3[k] = parm_data->dihedral_data[i].atom3[k];
        more[counter].atom4[k] = parm_data->dihedral_data[i].atom4[k];
      }
      ++counter;
    }
  }
  /* remove the old array and replace it with our new one - this needs to be done so   *
   * all our dihedrals stay in order.                                                  */
  free(parm_data->dihedral_data);
  parm_data->dihedral_data = more;
  
  printf("*  Successfully set %d dihedral terms.\n", n);
  return SUCCESS;
}


