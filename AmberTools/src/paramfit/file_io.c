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

/*file_io.c*/
/*
  This module is responsible for all file io with the exception of
  processing the prmtop file

  int read_job_control_file(global_options_struct *global_options)
  int_retval=read_mdcrd(global_options_struct *global_options, coords_struct *coords_data,  int structure);
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "function_def.h"
#include "constants.h"

int read_job_control_file(global_options_struct *global_options)
{
   int  int_retval;
   int num_lines; // number of lines in job control file
   int longest_line;
   int number_settings;
   int line_length;
   register int temp_count;
   char *line=NULL;
   char temp_char;
   FILE *fptr=NULL;
   
   num_lines=0;
   longest_line=0;
   temp_count=0;
   number_settings=0;
   /*We open the file in *ptr_job_control_filename and then set up our settings*/

   /*We read a line at a time from the control file*/
   if (global_options->VERBOSITY>=MEDIUM)
      printf(" Reading job control file: %s\n",global_options->job_control_filename);

   /*Read the entire settings file into memory and then work from there*/
   // Open the file
   if((fptr=fopen(global_options->job_control_filename,"r"))==NULL)
   {
      file_open_failure("read_job_control_file", global_options->job_control_filename);
      return FILE_OPEN_FAIL;
   }

  // Count the number of lines in the file and the length of the longest one
  while (fscanf(fptr,"%c",&temp_char)!=EOF)
  {
  if (temp_char == '\n')
  {
    ++num_lines;
    if (temp_count > longest_line)
    {
      longest_line = temp_count;
      temp_count = 0;
    }
  }
  else
    ++temp_count;
  }
  rewind(fptr);
  
  // Allocate memory - we will read in one line at a time so array should be
  //                   as big as the longest line plus a null character
    line = (char*)calloc((longest_line+1), sizeof(char));
    if (line==NULL)
    {
      malloc_failure_char("read_job_control_file", "line", longest_line+1);
      return ALLOC_FAIL;
    }
    
  // Go through the file and grab and process one line at a time
  for (temp_count=0; temp_count < num_lines; ++temp_count)
  {
    fgets(line, longest_line+1, fptr);
    
    // Get the length of the line
    int_retval=0;
    line_length=0;
    while (int_retval < longest_line+1)
    {
      if (line[int_retval] != '\n' && line[int_retval] != '\0')
	++line_length;
      
      ++int_retval;
    }
    
    // Process the line if it is a valid setting
    if (line[0] != '#' && line_length > 2)
    {
	int_retval=process_job_control_setting(line, line_length, &number_settings, global_options);
      if (int_retval!=SUCCESS)
      {
	if (int_retval==INVALID_FORMAT)
	  printf("!  Error near line %d of job control file, Unknown variable name - ignored.\n",temp_count);
	if (int_retval==INVALID_DATA)
	  printf("!  Error near line %d of job control file, Invalid data for selected variable - ignored.\n",temp_count);
	if (int_retval==INVALID_LINE)
	  printf("!  Error near line %d of settings file, no '=' found - ignored.\n",temp_count);
	if (int_retval==ALLOC_FAIL)
	{
	  printf("!  MALLOC FAILURE in process_setting(), attempting to ignore...\n");
	  printf("!  WARNING - CHECK PARAMETERS - RESULTS MAY BE INVALID\n");
	}
      }
//       else
// 	++number_settings;
    }
    
    // Clean out the line buffer
    for (int_retval=0; int_retval < longest_line+1; ++int_retval)
      line[int_retval] = '\0';
  }
  
  fclose(fptr);

  if (global_options->VERBOSITY>=MEDIUM)
    printf(" Job Control: Read a total of %d lines from job_control file. %d options set.\n\n",temp_count+1, number_settings);
  
  free(line);
  line = NULL;
  
  return(SUCCESS);
}

int read_mdcrd(global_options_struct *global_options, coords_struct *coords_data)
{
  int i;
  FILE *fptr;
  char mdcrd_temp[BUFFER_SIZE];
  int retval;
  int structure;

  /*Check we can open the file*/
  if((fptr=fopen(global_options->mdcrd_filename,"r"))==NULL)
  {
     file_open_failure("read_mdcrd", global_options->mdcrd_filename);
     return FILE_OPEN_FAIL;
  }

  retval=s_getline( mdcrd_temp, BUFFER_SIZE, fptr );  /*blank line*/
  if ( retval==FILE_READ_FAIL )
  {
      printf("*** ERROR IN read_mdcrd SUBROUTINE\n");
      printf("*** HIT EOF WHILE READING INITIAL BLANK LINE OF MDCRD FILE: %s\n",global_options->mdcrd_filename);
      return FILE_READ_FAIL;
  }

  /*Read in all of the structures at once */
  for (structure=0; structure<global_options->NSTRUCTURES; ++structure)
  {
    for (i=0;i<(global_options->NATOMS); ++i)
    {
      /*We should keep reading 3 floats at a time in the order x,y,z
        note, currently there is no allowance for box info in the mdcrd file
      */
      retval=fscanf(fptr,"%lf",&coords_data[structure].x_coord[i]);
      if (retval==EOF)
      {
        printf("*** ERROR IN READ_MDCRD - HIT EOF WHILE READING ELEMENT: %d\n",i);
        return FILE_READ_FAIL;
      }
      else if (retval!=1)
      {
        printf("*** ERROR IN READ_MDCRD - FAILED TO READ X_COORD FOR ELEMENT: %d\n",i);
        return INVALID_FORMAT;
      }
      retval=fscanf(fptr,"%lf",&coords_data[structure].y_coord[i]);
      if (retval==EOF)
      {
        printf("*** ERROR IN READ_MDCRD - HIT EOF WHILE READING ELEMENT: %d\n",i);
        return FILE_READ_FAIL;
      }
      else if (retval!=1)
      {
        printf("*** ERROR IN READ_MDCRD - FAILED TO READ Y_COORD FOR ELEMENT: %d\n",i);
        return INVALID_FORMAT;
      }
      retval=fscanf(fptr,"%lf",&coords_data[structure].z_coord[i]);
      if (retval==EOF)
      {
        printf("*** ERROR IN READ_MDCRD - HIT EOF WHILE READING ELEMENT: %d\n",i);
        return FILE_READ_FAIL;
      }
      else if (retval!=1)
      {
        printf("*** ERROR IN READ_MDCRD - FAILED TO READ Z_COORD FOR ELEMENT: %d\n",i);
        return INVALID_FORMAT;
      }
    }
  }

 fclose(fptr);
 
 return SUCCESS;
}

int read_qm_energy(global_options_struct *global_options, coords_struct *coords_data)
{
   /* This routine will read each of the energies in the QM_data file and puts the relevant value, converted to
      KCal/mol if necessary in the coords_data array*/

   /*Note, the QM data file should consist of NSTRUCTURES worth of doubles with one value per line*/
  int i;
  FILE *fptr;
  int retval;
  double conversion_factor;

  if (global_options->QM_ENERGY_UNITS==HARTREE)
    conversion_factor = HARTREE_TO_KCALMOL;
  else if (global_options->QM_ENERGY_UNITS==KJMOL)
    conversion_factor = KJMOL_TO_KCALMOL;
  else
    conversion_factor = 1.0;

  /*Initially perform a sanity check to see if we are reading at least 1 structure*/
  if ( global_options->NSTRUCTURES < 1 )
  {
     printf("*** ERROR IN read_qm_energy SUBROUTINE\n");
     printf("*** INVALID NUMBER OF STRUCTURES.\n");
     printf("*** NUMBER OF STRUCTURES SPECIFIED WAS %d WHICH IS NOT > 0\n",global_options->NSTRUCTURES);
     return INVALID_DATA;
  }

  /*Check we can open the file*/
  if((fptr=fopen(global_options->energy_filename,"r"))==NULL)
  {
     file_open_failure("read_qm_energy", global_options->energy_filename);
     return FILE_OPEN_FAIL;
  }

  /*Now we read in NSTRUCTURES worth of QM energies and check we don't hit EOF while running*/
  for (i=0;i<global_options->NSTRUCTURES;++i)
  {
     retval=fscanf(fptr,"%lf",&coords_data[i].energy);
     if (retval==EOF)
     {
        printf("*** ERROR IN READ_QM_ENERGY - HIT EOF WHILE READING ELEMENT: %d\n",i);
        return FILE_READ_FAIL;
     }
     else if (retval!=1)
     {
        printf("*** ERROR IN READ_QM_ENERGY - FAILED TO READ ENERGY FOR ELEMENT: %d\n",i);
        return INVALID_FORMAT;
     }
     coords_data[i].energy=coords_data[i].energy*conversion_factor;
  }
  fclose(fptr);
  /* Sort the structures in order of increasing energy */
  qsort(coords_data, global_options->NSTRUCTURES, sizeof(coords_struct), compare_energy);
  
  return SUCCESS;

}

int read_parameter_file(global_options_struct *global_options, parm_struct *parm_data)
{

  FILE *fptr;
  char *line = calloc(128, sizeof(char));
  int line_number = 0;
  short int verify = 0;
  int i = 0;
  short int temp1, temp2, temp3;
  
  // open the file
  if((fptr=fopen(global_options->PARAMETER_FILE_NAME,"r"))==NULL)
  {
    file_open_failure("read_parameter_file", global_options->PARAMETER_FILE_NAME);
    free(line);
    return FILE_OPEN_FAIL;
  }
  
  // locate the section with bond information
  // Skip three lines to get to the section with bonds
  for (i=0; i<4; i++)
  {
    fgets(line,127,fptr);
    ++line_number;
  }
  
  for (i=0; i<parm_data->unique_bonds_found; ++i)
  {
    if (fscanf(fptr, "%hd %hd %hd", &verify, &temp1, &temp2) != 3)
    {
      printf("ERROR reading data from parameter file for bond #%d.\n", i);
      free(line);
      return INVALID_DATA;
    }
    parm_data->bond_data[i].DO_FIT_REQ = temp1;
    parm_data->bond_data[i].DO_FIT_KR = temp2;
  }
  
  // Skip 4 lines to get to the section with angles
  for (i=0; i<4; i++)
  {
    fgets(line,127,fptr);
    ++line_number;
  }
  
  // now put in the angle information
    for (i=0; i<parm_data->unique_angles_found; ++i)
    {
      if(fscanf(fptr, "%hd %hd %hd\n", &verify, &temp1, &temp2) != 3)
      {
        printf("ERROR reading data from parameter file for angle #%d.\n", i);
        free(line);
        return INVALID_DATA;
      }
      parm_data->angle_data[i].DO_FIT_KT = temp1;
      parm_data->angle_data[i].DO_FIT_THEQ = temp2;
    }
    
    // Skip 4 lines to get to the section with angles
    for (i=0; i<3; i++)
    {
      fgets(line,127,fptr);
      ++line_number;
    }
    
    // do the dihedrals
    for (i=0; i<parm_data->unique_dihedrals_found; ++i)
    {
      if (fscanf(fptr, "%hd %hd %hd %hd", &verify, &temp1, &temp2, &temp3) != 4) // number of arguments read
      {
        printf("ERROR reading data from parameter file for dihedral #%d.\n", i);
        free(line);
        return INVALID_DATA;
      }
      parm_data->dihedral_data[i].DO_FIT_KP = temp1;
      parm_data->dihedral_data[i].DO_FIT_NP = temp2;
      if (parm_data->dihedral_data[i].improper == YES) // don't fit improper dihedral phase
        parm_data->dihedral_data[i].DO_FIT_PHASE = NO;
      else
        parm_data->dihedral_data[i].DO_FIT_PHASE = temp3;
    }
  
  free(line);
  fclose(fptr);
  if (global_options->VERBOSITY>=MEDIUM)
    printf("   Prmtop     (info): Successfully read in saved parameter information\n");
  return SUCCESS;
}

/* This will read the gaussian output files for all structures. 
 * The gaussian input files need to have been generated with CREATE_INPUT so the atom numbering is
 * consistent. This also lets us assume that the files are named QMFILEOUTSTART<n>QMFILEOUTEND.out
 */
int read_gaussian_forces(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data)
{ 
  // Determine if a unit conversion is necessary- internal force units are kcal/mol-A
  double units;
  if (global_options->QM_FORCE_UNITS==HARTREE_BOHR)
    units = HARTREE_TO_KCALMOL*BOHR_TO_ANGSTROM;
  else if (global_options->QM_FORCE_UNITS==KCALMOL_ANGSTROM)
    units = 1.0;
  else {
    printf("*** ERROR: Unsupported force units for Gaussian input files.\n");
    return INVALID_DATA;
  }
    
  int i;
  for(i=0; i<global_options->NSTRUCTURES; ++i) {
    // Generate the filename to be opened
    char filename[1000];
    sprintf(filename, "%s%d%s", global_options->QMFILEOUTSTART, i, global_options->QMFILEOUTEND);
    
    // Open the file
    FILE *fptr = fopen(filename, "r");
    if(!fptr)
    {
      file_open_failure("read_gaussian_forces", filename);
      return FILE_OPEN_FAIL;
    }
    
    /* Read and discard lines until axes line is found, discard 4 lines, begin reading.
     * Output will look like this:
     *  ***** Axes restored to original set *****
     * -------------------------------------------------------------------
     * Center     Atomic                   Forces (Hartrees/Bohr)
     * Number     Number              X              Y              Z
     * -------------------------------------------------------------------
     *      1        1           0.011759647    0.003870626    0.010325240
     * ...
     */
    char line[200];
    int check;
    do {
      check = s_getline(line, 200, fptr);
      if ( check==FILE_READ_FAIL ) {
        printf("*** ERROR IN read_gaussian_forces SUBROUTINE\n");
        printf("*** Failed to find forces section in output file: %s\n\n", filename);
        printf("*** Check that this is the correct filename and that\n");
        printf("    the Gaussian job completed successfully!\n");
        return FILE_READ_FAIL;
      }
    } while( !strstr(line, "Axes restored to original set") ); 

    // Discard 3 lines of table header
    int atom;
    for (atom=0; atom<=3; ++atom) {
      s_getline(line, 200, fptr);
    }
    
    // Read in the data for each atom
    for (atom=0; atom<global_options->NATOMS; ++atom) 
    {
      int element = find_atomic_number_from_parm(parm_data, atom+1);
      if ( fscanf(fptr, "%*d %d %15lf %15lf %15lf",  &check, &coords_data[i].force[atom].x,
                  &coords_data[i].force[atom].y, &coords_data[i].force[atom].z) != 4 ) {
        printf("*** ERROR IN read_gaussian_forces SUBROUTINE\n");
        printf("*** Error in reading forces for atom %d\n", atom);
        printf("*** Filename: %s\n", filename);
        return FILE_READ_FAIL;
      }
      if ( check != element ) {
        printf("*** ERROR IN read_gaussian_forces SUBROUTINE\n");
        printf("*** Unexpected element in file: %s\n", filename);
        printf("*** Expected: %d \t Received: %d \n\n", element, check);
        printf("*** Check that the input files used for these output files were\n");
        printf("    generated with paramfit for this same prmtop.\n");
        return FILE_READ_FAIL;
      }
      // Unit conversion
      coords_data[i].force[atom].x *= units;
      coords_data[i].force[atom].y *= units;
      coords_data[i].force[atom].z *= units;
    }
    fclose(fptr);
  }
  return SUCCESS;
}


