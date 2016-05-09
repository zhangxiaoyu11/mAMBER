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

/*write_input.c

  Contains routines for writing the various different formats for input files used to generate
  the energy data
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "function_def.h"
#include "constants.h"

int write_input_gaussian(global_options_struct *global_options, parm_struct *parm_data, coords_struct *current_struct, int num, FILE *fptr)
{

  int i;
  int atomic_number;
  /*The output file should already be open here so we just start writing the data*/
  /* First open the header file and copy it to the top of the input file */
  if (global_options->QMHEADER)
  {
    FILE *header = fopen(global_options->QMHEADER, "r");
    if (!header)
    {
      printf("! ERROR opening QM header file at %s\n", global_options->QMHEADER);
      return INVALID_DATA;
    }

    // Get the size of the header
    fseek(header, 0, SEEK_END);
    long size = ftell(header);
    rewind(header);
    
    // allocate buffer, read in data, and write it to the file
    unsigned char *buffer=(unsigned char *)malloc((size_t)size);
    fread(buffer, size, 1, header);
    fwrite(buffer, size, 1, fptr);
    
    // clean up
    free(buffer);
    fclose(header);
  }
  
  /*Now in gaussian we need to write a blank line, a title and then a blank line, use the current_struct for the title*/
  fprintf(fptr,"\n");
  fprintf(fptr,"Structure: %d\n",num);
  fprintf(fptr,"\n");
  /*Now the charge and multiplicity*/
  fprintf(fptr,"  %d  %d\n",global_options->QM_SYSTEM_CHARGE,global_options->QM_SYSTEM_MULTIPLICITY);

  /*Loop over all atoms*/
  for (i=0;i<global_options->NATOMS;++i)
  {
    /*Stage 1 - get the element*/
    /*we really only have 2 ways to do this, either from the mass or from the first letter of the
      atom name - both are unreliable*/
    /*Lets just use the first letter of the atom name - note, this can give misreads with for example
      N and Na, we can try and correct for these by checking the mass.
    */
    atomic_number=find_atomic_number_from_parm(parm_data, i+1);
    if (atomic_number==UNKNOWN_ELEMENT)
    {
      /*Element not recognised*/
      printf("\n*** ERROR IN WRITE_INPUT FROM CALL TO FIND_ATOMIC_NUMBER_FROM_PARM\n");
      printf("*** UNKNOWN ELEMENT FOR ATOM NAME: %s WITH MASS: %f\n",parm_data->atom[i].igraph, parm_data->atom[i].amass);
      printf("*** CHECK elements.c\n");
      return (UNKNOWN_ELEMENT);
    }
    /*Now write the element as its symbol*/
    fprintf(fptr,"  ");
    print_atomic_number_as_symbol(fptr, atomic_number);
    fprintf(fptr,"   ");
    /*Now the coordinates*/
    fprintf(fptr,"%10.6f %10.6f %10.6f\n",
	    current_struct->x_coord[i], current_struct->y_coord[i], current_struct->z_coord[i]);
  }
    
  
  /*Print a final carriage return*/
  fprintf(fptr,"\n");
  return (SUCCESS);
}

int write_input_parameters(global_options_struct *global_options, parm_struct *parm_data)
{
  FILE *fptr;
  int i;

  /*Now fill the filename_to_write array*/
  printf("   Filename to be written is: %s\n",global_options->PARAMETER_FILE_NAME);
  
  if((fptr=fopen(global_options->PARAMETER_FILE_NAME,"w"))==NULL)
  {
    file_open_failure("create_input", global_options->PARAMETER_FILE_NAME);
    return FILE_OPEN_FAIL;
  }
  // print a header line indicating the data that go with this particular file
  fprintf(fptr,"#AUTO-GENERATED PARAMETER FILE FOR %s - DO NOT MODIFY\n", global_options->prmtop_filename);
  // put information about the desired bond parameters into the file
  fprintf(fptr, "#\n# BOND INFORMATION:\n");
  fprintf(fptr, "#### BOND\tREQ\tKR ####\n");
  for (i=0; i<parm_data->unique_bonds_found; i++)
    fprintf(fptr, "\t%d\t%d\t%d\n", i, parm_data->bond_data[i].DO_FIT_REQ, parm_data->bond_data[i].DO_FIT_KR);
  
  // put information about the desired angle parameters into the file
  fprintf(fptr, "#\n# ANGLE PARAMETERS:\n");
  fprintf(fptr, "#### ANGLE\tKT\tTHEQ ####\n");
  for (i=0; i<parm_data->unique_angles_found; i++)
    fprintf(fptr, "\t%d\t%d\t%d\n", i, parm_data->angle_data[i].DO_FIT_KT, parm_data->angle_data[i].DO_FIT_THEQ);
  
  // put information about the desired dihedral parameters into the file
  fprintf(fptr, "#\n# DIHEDRAL PARAMETERS:\n");
  fprintf(fptr, "#### DIHEDRAL\tKP\tNP\tPHASE ####\n");
  for (i=0; i<parm_data->unique_dihedrals_found; i++)
    fprintf(fptr, "\t%d\t%d\t%d\t%d\n", i, parm_data->dihedral_data[i].DO_FIT_KP, parm_data->dihedral_data[i].DO_FIT_NP, 
	    parm_data->dihedral_data[i].DO_FIT_PHASE);
	    
  fclose(fptr);
  
  return SUCCESS;
}

int write_frcmod(global_options_struct *global_options, parm_struct *parm_data)
{
  FILE *fptr;
  int i;
//   int j;
//   int k;
//   int num_type;
//   int dihedral_names;
//   dihedral_data_struct **dihedrals=NULL;
//   dihedral_data_struct **names=NULL;
  
  /* Make sure the settings are correct */
  if((fptr=fopen(global_options->WRITE_FRCMOD,"w"))==NULL)
  {
    file_open_failure("write_ffrcmod", global_options->WRITE_FRCMOD);
    return FILE_OPEN_FAIL;
  }
  if (global_options->VERBOSITY>=MEDIUM)
    printf(" * Saving ffrcmod file to %s\n", global_options->WRITE_FRCMOD);
  
  // print a header line indicating the data that go with this particular file
  fprintf(fptr, "Generated frcmod with paramfit for %s\n", global_options->prmtop_filename);
  
  // print the bond information, if any
  if (global_options->BOND_PARAMS > 0)
    fprintf(fptr,"BOND\n");
  for (i=0; i<parm_data->unique_bonds_found; ++i)
  {
    if (parm_data->bond_data[i].DO_FIT_KR==YES || parm_data->bond_data[i].DO_FIT_REQ==YES)
    {
      fprintf(fptr, "%s-%s %14.4f %14.4f\n", parm_data->bond_data[i].atom_type1, parm_data->bond_data[i].atom_type2,
                                               parm_data->bond_data[i].rk, parm_data->bond_data[i].req);
    }
  }
  fprintf(fptr, "\n");
  
  // now deal with the angles
  if (global_options->ANGLE_PARAMS > 0)
    fprintf(fptr, "ANGL\n");
  for (i=0; i<parm_data->unique_angles_found; ++i)
  {
    if (parm_data->angle_data[i].DO_FIT_KT==YES || parm_data->angle_data[i].DO_FIT_THEQ==YES)
    {
      fprintf(fptr, "%s-%s-%s %14.4f %14.4f\n", parm_data->angle_data[i].atom_type1, parm_data->angle_data[i].atom_type2,
                                                      parm_data->angle_data[i].atom_type3, parm_data->angle_data[i].tk,
                                                      parm_data->angle_data[i].teq*RADIAN_TO_DEGREE);
    }
  }
  fprintf(fptr, "\n");
  
  /* finally dihedrals - these are a little more complicated */
  
  // count how many unique dihedral names we have and allocate accordingly
//  names = (dihedral_data_struct**)malloc(sizeof(dihedral_data_struct*));                                            // 1d array of pointers to unique dihedral types
//  if (names==NULL)
//  {
//    printf("!  ERROR in counting dihedrals.\n");
//    return ALLOC_FAIL;
//  }
//  dihedral_names = 1;                                                                                                                      // number of unique dihedrals found
//  
//  names[0] = &parm_data->dihedral_data[0];
//  for (i=1; i<parm_data->unique_dihedrals_found; ++i)
//  {
//    j=0;
//    while (j<dihedral_names && dihedral_types_equal(&parm_data->dihedral_data[i], names[j])==NO) // compare to all unique names found so far
//      ++j;
//    if (j==dihedral_names)       // dihedral name has not yet been found
//    {
//      ++dihedral_names;
//      names = (dihedral_data_struct**)realloc(names, dihedral_names*sizeof(struct _dihedral_data_struct*));
//      if (names==NULL)
//      {
//        printf("!  ERROR in counting dihedrals.\n");
//        return ALLOC_FAIL;
//      }
//      names[dihedral_names-1] = &parm_data->dihedral_data[i];
//    }
//  }

// print DIHE
/* 
 // Now count each dihedral name, sort it in order of decreasing periodicity, and write it to the file 
  for (i=0; i<dihedral_names; ++i)
  {
    num_type=0;
    // Count how many of this type and put it in an array
    for (j=0; j<parm_data->unique_dihedrals_found; ++j)
    {
      if (dihedral_types_equal(names[i], &parm_data->dihedral_data[j])==YES)
      {
        ++num_type;
        dihedrals = (dihedral_data_struct**)realloc(dihedrals, num_type*sizeof(struct _dihedral_data_struct*));
        if (dihedrals==NULL)
        {
          printf("!  ERROR in sorting dihedrals for writing the ffrcmod.\n");
          return ALLOC_FAIL;
        }
        dihedrals[num_type-1] = &parm_data->dihedral_data[j];
      }
    }

    // Sort this dihedral type in order of decreasing periodicity
    for (j=1; j<num_type; ++j)
    {
      dihedral_data_struct *temp = dihedrals[j];
   
      k = j-1;
      while (k>=0 && dihedrals[k]->pn < temp->pn)
      {
        dihedrals[k+1] = dihedrals[k];
        --k;
      }
      dihedrals[k+1] = temp;
      temp = NULL;
    }
    // Write this dihedral type to the ffrcmod file
    for (j=0; j<num_type; ++j)
    {
       if (dihedrals[j]->DO_FIT_KP || dihedrals[j]->DO_FIT_NP || dihedrals[j]->DO_FIT_PHASE)
       {
        if (j!=num_type-1) // if it's not the last one, pn should be negative to show existence of other terms
          dihedrals[j]->pn*=-1.0;
        fprintf(fptr, "%s-%s-%s-%s    1 %14.4f %14.3f %14.3f\n", dihedrals[j]->atom_type1, dihedrals[j]->atom_type2,
                dihedrals[j]->atom_type3, dihedrals[j]->atom_type4,
                dihedrals[j]->pk, dihedrals[j]->phase*RADIAN_TO_DEGREE, // frcmod files have degrees
                dihedrals[j]->pn);
      } 
    }
   }
*/

 // Now the dihedrals
 if (global_options->DIHEDRAL_PARAMS > 0)
   fprintf(fptr, "DIHE\n");
 int impropers = 0;
 for (i=0; i<parm_data->unique_dihedrals_found; ++i)
 {
   if (parm_data->dihedral_data[i].improper==NO &&
       (parm_data->dihedral_data[i].DO_FIT_KP==YES || parm_data->dihedral_data[i].DO_FIT_NP==YES ||
       parm_data->dihedral_data[i].DO_FIT_PHASE==YES) )
   {
     fprintf(fptr, "%s-%s-%s-%s    1 %14.4f %14.3f %14.0f\n", parm_data->dihedral_data[i].atom_type1,
             parm_data->dihedral_data[i].atom_type2, parm_data->dihedral_data[i].atom_type3,
             parm_data->dihedral_data[i].atom_type4, parm_data->dihedral_data[i].pk,
             parm_data->dihedral_data[i].phase*RADIAN_TO_DEGREE, parm_data->dihedral_data[i].pn);
   }
   else
     ++impropers;
 }
 fprintf(fptr, "\n");
 
 if (impropers > 0)
   fprintf(fptr, "IMPR\n");
 for (i=0; i<parm_data->unique_dihedrals_found; ++i)
 {
   if (parm_data->dihedral_data[i].improper==YES &&
       (parm_data->dihedral_data[i].DO_FIT_KP==YES || parm_data->dihedral_data[i].DO_FIT_NP==YES ||
       parm_data->dihedral_data[i].DO_FIT_PHASE==YES) )
   {
     fprintf(fptr, "%s-%s-%s-%s     %14.4f %14.3f %14.3f\n", parm_data->dihedral_data[i].atom_type1,
             parm_data->dihedral_data[i].atom_type2, parm_data->dihedral_data[i].atom_type3,
             parm_data->dihedral_data[i].atom_type4, parm_data->dihedral_data[i].pk,
             parm_data->dihedral_data[i].phase*RADIAN_TO_DEGREE, parm_data->dihedral_data[i].pn);
   }
 }
 fprintf(fptr, "\n");
 fclose(fptr);
 
//  free(dihedrals);
//  free(names);
//  
//  dihedrals = NULL;
//  names = NULL;
 
 return SUCCESS;
}

int write_input_adf(global_options_struct *global_options, parm_struct *parm_data, coords_struct *current_struct, int num, FILE *fptr)
{
  // write some basic header information
  fprintf(fptr, "TITLE structure %d generated with paramfit\n", num);
  if (global_options->QMHEADER)
    fprintf(fptr, "INLINE %s\n!\n", global_options->QMHEADER);
  fprintf(fptr, "UNITS\n length Angstrom\n angle Radian\nEnd\n!\n");
  
  // write the charge and multiplicity
  fprintf(fptr, "CHARGE %d %d\n!\n", global_options->QM_SYSTEM_CHARGE, global_options->QM_SYSTEM_MULTIPLICITY);
  fprintf(fptr, "TOTALENERGY\n!\n");
  
  // write the cartesian matrix
  fprintf(fptr, "ATOMS Cartesian\n");
  int i, atomic_number;
  for (i=0;i<global_options->NATOMS;++i)
  {
    atomic_number=find_atomic_number_from_parm(parm_data, i+1);
    if (atomic_number==UNKNOWN_ELEMENT)
    {
      /*Element not recognised*/
      printf("\n*** ERROR IN WRITE_INPUT FROM CALL TO FIND_ATOMIC_NUMBER_FROM_PARM\n");
      printf("*** UNKNOWN ELEMENT FOR ATOM NAME: %s WITH MASS: %f\n",parm_data->atom[i].igraph, parm_data->atom[i].amass);
      printf("*** CHECK elements.c\n");
      return (UNKNOWN_ELEMENT);
    }
    fprintf(fptr, " %i ", i);
    print_atomic_number_as_symbol(fptr, atomic_number);
    fprintf(fptr, " %.6f %.6f %.6f\n",
	    current_struct->x_coord[i], current_struct->y_coord[i], current_struct->z_coord[i]);
  }
  fprintf(fptr, "End\n");
  return SUCCESS;
}

int write_input_gamess(global_options_struct *global_options, parm_struct *parm_data, coords_struct *current_struct, int num, FILE *fptr)
{
  // include the user's header, which will probably be the $SYSTEM section and the $BASIS section
  time_t now = time(NULL);
  fprintf(fptr, "! Gamess input file for structure %d automatically generated with paramfit\n", num);
  fprintf(fptr, "! Generated on %s for prmtop at %s", asctime(localtime(&now)), global_options->prmtop_filename);
  if (global_options->QMHEADER)
  {
    FILE *header = fopen(global_options->QMHEADER, "r");
    if (!header)
    {
      printf("! ERROR opening QM header file at %s\n", global_options->QMHEADER);
      return INVALID_DATA;
    }
    // Get the size of the header
    fseek(header, 0, SEEK_END);
    long size = ftell(header);
    rewind(header);
    
    // allocate buffer, read in data, and write it to the file
    unsigned char *buffer=(unsigned char *)malloc((size_t)size);
    fread(buffer, size, 1, header);
    fwrite(buffer, size, 1, fptr);
    
    // clean up
    free(buffer);
    fclose(header);
  }
  
  // write the basic $CONTRL section
  fprintf(fptr, " $CONTRL\n");
  fprintf(fptr, "  RUNTYP=ENERGY\n");
  fprintf(fptr, "  COORD=UNIQUE\n");
  fprintf(fptr, "  ICHARG=%d\n  MULT=%d\n", global_options->QM_SYSTEM_CHARGE, global_options->QM_SYSTEM_MULTIPLICITY);
  fprintf(fptr, " $END\n!\n");
  
  // write the $DATA section
  fprintf(fptr, " $DATA\n");
  fprintf(fptr, "  Structure %d generated by paramfit\n", num);
  fprintf(fptr, "  C1\n"); // don't specify symmetry since molecule could be anything
  
  // write all the atoms
  int i, atomic_number;
  for (i=0;i<global_options->NATOMS;++i)
  {
    atomic_number=find_atomic_number_from_parm(parm_data, i+1);
    if (atomic_number==UNKNOWN_ELEMENT)
    {
      /*Element not recognised*/
      printf("\n*** ERROR IN WRITE_INPUT FROM CALL TO FIND_ATOMIC_NUMBER_FROM_PARM\n");
      printf("*** UNKNOWN ELEMENT FOR ATOM NAME: %s WITH MASS: %f\n",parm_data->atom[i].igraph, parm_data->atom[i].amass);
      printf("*** CHECK elements.c\n");
      return (UNKNOWN_ELEMENT);
    }
    fprintf(fptr, "  ");
    print_atomic_number_as_symbol(fptr, atomic_number);
    fprintf(fptr, " %d %10.6f %10.6f %10.6f\n", atomic_number,
	    current_struct->x_coord[i], current_struct->y_coord[i], current_struct->z_coord[i]);
  }
  fprintf(fptr, " $END\n");
  return SUCCESS;
}

int write_energy(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data, int generation)
{
  FILE *fptr;
  int i;
  double temp_energy;
  char filename[1000];
  
  // create the file
  if (generation >= 0)
    sprintf(filename, "%d_%s", generation, global_options->WRITE_ENERGY);
  else
    sprintf(filename, "%s", global_options->WRITE_ENERGY);
  
  if((fptr=fopen(filename,"w"))==NULL)
  {
    file_open_failure("write_energies", global_options->WRITE_ENERGY);
    return FILE_OPEN_FAIL;
  }
  if (global_options->VERBOSITY>=MEDIUM)
    printf(" * Saving energy file with %d structures to %s\n", global_options->NSTRUCTURES, filename);

  // print a header line indicating the data that go with this particular file
  fprintf(fptr, "# Final energies generated by paramfit for %s\n", global_options->prmtop_filename);
  fprintf(fptr, "# Num\tAmber+K\t\tQuantum\n");
  
  // evaluate the energy for all of the input structures and print it to the file
  // here we put K back in since it represents the difference between quantum and classical so everything lines up
  for (i=0; i<global_options->NSTRUCTURES; ++i)
  {
    temp_energy = eval_amber_std_for_single_struct(global_options, parm_data, &coords_data[i]);
    fprintf(fptr, "%d\t%10.5f\t%10.5f\n", i, (temp_energy+global_options->K), coords_data[i].energy);
  }
  fclose(fptr);
  
  return SUCCESS;
}

int write_prmtop(global_options_struct *global_options, parm_struct *parm_data)
{
  FILE *fptr;
  int i, j, index;
  
  if ((fptr=fopen(global_options->WRITE_PRMTOP, "w"))==NULL)
  {
    file_open_failure("write_prmtop", global_options->WRITE_PRMTOP);
    return FILE_OPEN_FAIL;
  }
  
  printf(" * Saving prmtop to %s\n", global_options->WRITE_PRMTOP);
  
  // Print header information
  fprintf(fptr, "\n%%FLAG TITLE\n");
  fprintf(fptr, "%%FORMAT(20a4)\n");
  fprintf(fptr, "Generated prmtop with paramfit for %s\n", global_options->prmtop_filename);
  
  // Print the flag pointers
  fprintf(fptr, "%%FLAG POINTERS\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  fprintf(fptr, "%8i%8i%8i%8i%8i%8i%8i\n",//%8i%8i%8i\n", 
          parm_data->NTOTAT, parm_data->NTYPES, parm_data->NBONH, parm_data->NBONA, parm_data->NTHETH,
          parm_data->NTHETA, parm_data->NPHIH);//, parm_data->NPHIA, parm_data->JHPARM, parm_data->JPARM);
  fprintf(fptr, "%8i%8i%8i%8i%8i%8i%8i%8i%8i%8i\n", 
          parm_data->NEXT, parm_data->NTOTRS, parm_data->MBONA, parm_data->MTHETS, parm_data->MPHIA,
          parm_data->MUMBND, parm_data->MUMANG, parm_data->MPTRA, parm_data->NATYP, parm_data->NHB);
  fprintf(fptr, "%8i%8i%8i%8i%8i%8i%8i\n",//%8i%8i%8i\n", 
          parm_data->IFPERT, parm_data->NBPER, parm_data->NGPER, parm_data->NDPER, parm_data->MBPER,
          parm_data->MGPER, parm_data->MDPER);// DOUBLE CHECK parm_data->IFBOX, parm_data->NMXRS, parm_data->IFCAP);
  fprintf(fptr, "%8i\n", parm_data->NUMEXTRA);
  
  // Print atom names
  fprintf(fptr, "%%FLAG ATOM_NAME\n");
  fprintf(fptr, "%%FORMAT(20a4)\n");
  for (i=0; i<parm_data->NTOTAT; ++i)
  {
    fprintf(fptr, "%4s", parm_data->atom[i].igraph);
  }
  
  // Print atom charges
  fprintf(fptr, "\n%%FLAG CHARGE\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<parm_data->NTOTAT; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->atom[i].chrg);
    if (i%5==1&&i<parm_data->NTOTAT-1)
      fprintf(fptr, "\n");
  }
  
  // Print atom masses
  fprintf(fptr, "\n%%FLAG MASS\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<parm_data->NTOTAT; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->atom[i].amass);
    if (i%5==1&&i<parm_data->NTOTAT-1)
      fprintf(fptr, "\n");
  }
  
  // Print atom type indices
  fprintf(fptr, "\n%%FLAG ATOM_TYPE_INDEX\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  for (i=0; i<parm_data->NTOTAT; ++i)
  {
    fprintf(fptr, "%8i", parm_data->atom[i].iac);
    if (i%10==1&&i<parm_data->NTOTAT-1)
      fprintf(fptr, "\n");
  }
  
  // Print number of excluded atoms
  fprintf(fptr, "\n%%FLAG NUMBER_EXCLUDED_ATOMS\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  for (i=0; i<parm_data->NTOTAT;++i)
  {
    fprintf(fptr, "%8i", parm_data->atom[i].numex);
    if (i%10==1&&i<parm_data->NTOTAT-1)
      fprintf(fptr, "\n");
  }
  
  // Print nonbonded parameter index
  fprintf(fptr, "\n%%FLAG NONBONDED_PARM_INDEX\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  for (i=0; i<parm_data->NTYPES*parm_data->NTYPES; ++i)
  {
    fprintf(fptr, "%8i", parm_data->nno[i]);
    if (i%10==1&&i<parm_data->NTYPES*parm_data->NTYPES-1)
      fprintf(fptr, "\n");
  }
  
  // Print residue labels
  fprintf(fptr, "\n%%FLAG RESIDUE_LABEL\n");
  fprintf(fptr, "%%FORMAT(20A4)\n");
  for (i=0; i<parm_data->NTOTRS; ++i)
  {
    fprintf(fptr, "%4s", parm_data->residue[i].labres);
    if (i%20==1&&i<parm_data->NTOTRS-1)
      fprintf(fptr, "\n");
  }
  
  // Now the residue pointers
  fprintf(fptr, "\n%%FLAG RESIDUE_POINTER\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  for (i=0; i<parm_data->NTOTRS; ++i)
  {
    fprintf(fptr, "%8i", parm_data->residue[i].ipres);
    if (i%10==1&&i<parm_data->NTOTRS-1)
      fprintf(fptr, "\n");
  }
  
  /* Put the data that were fitted back into the arrays for each     *
   * individual bond, angle, and dihedral instead of just the unique *
   * ones. Loop through each one in the prmtop and find which one it *
   * belongs to in the fitted data structures.                       */
  // Start with bonds- look through ones with H then ones without H
  for (i=0; i<parm_data->unique_bonds_found; ++i)
  {
    // find number of each bond 
    for (j=0; j<parm_data->bond_data[i].number; ++j)
    {
      index = 0;
      int atom1 = ObfuscateAtom(parm_data->bond_data[i].atom1[j]);
      int atom2 = ObfuscateAtom(parm_data->bond_data[i].atom2[j]);
      
      while ( index < parm_data->NBONH && ( parm_data->pbondH[index].ib != atom1 || parm_data->pbondH[index].jb != atom2 ))        
      {
        ++index;
      }
      if (index == parm_data->NBONH) // not found in bonds with H
      {
        index = 0;
        while ( (parm_data->pbond[index].ib != atom1 || parm_data->pbond[index].jb != atom2))
        {
          ++index;
        }
        
        index = parm_data->pbond[index].icb-1;
      }
      else
        index = parm_data->pbondH[index].icb-1;
      
      // now index has the bond we're looking for. Copy its parameters into rk and req
      parm_data->rk[index] = parm_data->bond_data[i].rk;
      parm_data->req[index] = parm_data->bond_data[i].req;
    }
  }
  for (i=0; i<parm_data->unique_angles_found; ++i) // now for the angles
  {
    for (j=0; j<parm_data->angle_data[i].number; ++j)
    {
      int atom1 = ObfuscateAtom(parm_data->angle_data[i].atom1[j]);
      int atom2 = ObfuscateAtom(parm_data->angle_data[i].atom2[j]);
      int atom3 = ObfuscateAtom(parm_data->angle_data[i].atom3[j]);
      index = 0;
      while (index < parm_data->NTHETH &&
        (parm_data->pangleH[index].it != atom1 || parm_data->pangleH[index].jt != atom2 || parm_data->pangleH[index].kt != atom3))
        ++index;
      if (index == parm_data->NTHETH)
      {
        index=0;
        while (parm_data->pangle[index].it != atom1 || parm_data->pangle[index].jt != atom2 || parm_data->pangle[index].kt != atom3)
          ++index;
        index = parm_data->pangle[index].ict-1;
      }
      else
        index = parm_data->pangleH[index].ict-1;
      // copy the parameters
        parm_data->tk[index] = parm_data->angle_data[i].tk;
        parm_data->teq[index] = parm_data->angle_data[i].teq;
    }
  }
  
  for (i=0; i<parm_data->unique_dihedrals_found; ++i)
  {
    for (j=0; j<parm_data->dihedral_data[i].number; ++j)
    {
      int atom1 = ObfuscateAtom(parm_data->dihedral_data[i].atom1[j]);
      int atom2 = ObfuscateAtom(parm_data->dihedral_data[i].atom2[j]);
      int atom3 = ObfuscateAtom(parm_data->dihedral_data[i].atom3[j]);
      int atom4 = ObfuscateAtom(parm_data->dihedral_data[i].atom4[j]);
      index = 0;
      while (index < parm_data->NPHIH && (fabs(parm_data->pdihedralH[index].ip)!=atom1 || 
                                          fabs(parm_data->pdihedralH[index].jp)!=atom2 ||
                                          fabs(parm_data->pdihedralH[index].kp)!=atom3 || 
                                          fabs(parm_data->pdihedralH[index].lp)!=atom4))
      {
        ++index;
      }
      if (index == parm_data->NPHIH)
      {
        index=0;
        while ((fabs(parm_data->pdihedral[index].ip)!=atom1 || 
                fabs(parm_data->pdihedral[index].jp)!=atom2 ||
                fabs(parm_data->pdihedral[index].kp)!=atom3 || 
                fabs(parm_data->pdihedral[index].lp)!=atom4))
        {
          ++index;
        }
        index = parm_data->pdihedral[index].icp-1;
      }
      else
        index = parm_data->pdihedralH[index].icp-1;
      parm_data->pk[index] = parm_data->dihedral_data[i].pk;
      parm_data->pn[index] = parm_data->dihedral_data[i].pn;
      parm_data->phase[index] = parm_data->dihedral_data[i].phase;
    }
  }
  
  // The bond force constants
  fprintf(fptr, "\n%%FLAG BOND_FORCE_CONSTANT\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<parm_data->MUMBND; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->bond_data[i].rk);
    if (i%5==1&&i<parm_data->MUMBND-1)
      fprintf(fptr, "\n");
  }
  
  // Bond equilibrium values
  fprintf(fptr, "\n%%FLAG BOND_EQUIL_VALUE\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<parm_data->MUMBND; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->bond_data[i].req);
    if (i%5==1&&i<parm_data->MUMBND-1)
      fprintf(fptr, "\n");
  }
  
  // Angle force constants
  fprintf(fptr, "\n%%FLAG ANGLE_FORCE_CONSTANT\n");
  fprintf(fptr, "%%FORMAT(5E16.8\n");
  for (i=0; i<parm_data->MUMANG; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->angle_data[i].tk);
    if (i%5==1&&i<parm_data->MUMANG-1)
      fprintf(fptr, "\n");
  }
  
  // Angle equilibrium angle
  fprintf(fptr, "\n%%FLAG ANGLE_EQUIL_VALUE\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<parm_data->MUMANG; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->angle_data[i].teq);
    if (i%5==1&&i<parm_data->MUMANG-1)
      fprintf(fptr, "\n");
  }
  
  // Dihedral force constant
  fprintf(fptr, "\n%%FLAG DIHEDRAL_FORCE_CONSTANT\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<parm_data->MPTRA; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->dihedral_data[i].pk);
    if (i%5==1&&i<parm_data->MPTRA-1)
      fprintf(fptr, "\n");
  }
  
  // Dihedral periodicity
  fprintf(fptr, "\n%%FLAG DIHEDRAL_PERIODICITY\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<parm_data->MPTRA; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->dihedral_data[i].pn);
    if (i%5==1&&i<parm_data->MPTRA-1)
      fprintf(fptr, "\n");
  }
  
  // Dihedral phase
  fprintf(fptr, "\n%%FLAG DIHEDRAL_PHASE\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<parm_data->MPTRA; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->dihedral_data[i].phase);
    if (i%5==1&&i<parm_data->MPTRA-1)
      fprintf(fptr, "\n");
  }
  
  // Electrostatic scaling constant
//   fprintf(fptr, "%%FLAG SCEE_SCALE_FACTOR\n");
//   fprintf(fptr, "%%FORMAT(5E16.8)\n");
//   for (i=0; i<parm_data->MPTRA; ++i)
//   {
//   }

  // Solty coefficients, even though they're unused
  fprintf(fptr, "\n%%FLAG SOLTY\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<parm_data->NATYP; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->solty[i]);
    if (i%5==1&&i<parm_data->NATYP-1)
      fprintf(fptr, "\n");
  }
  
  // Lennard Jones coefficients
  index = parm_data->NTYPES*(parm_data->NTYPES+1)/2;
  fprintf(fptr, "\n%%FLAG LENNARD_JONES_ACOEF\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<index; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->cn1[i]);
    if (i%5==1&&i<index-1)
      fprintf(fptr, "\n");
  }
  fprintf(fptr, "\n%%FLAG LENNARD_JONES_BCOEF\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<index; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->cn2[i]);
    if (i%5==1&&i<index-1)
      fprintf(fptr, "\n");
  }
  
  // Bonds with hydrogen
  fprintf(fptr, "\n%%FLAG BONDS_INC_HYDROGEN\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  for (i=0; i<parm_data->NBONH; ++i)
  {
    fprintf(fptr, "%8i%8i%8i",
            parm_data->pbondH[i].ib, parm_data->pbondH[i].jb, parm_data->pbondH[i].icb);
    if ((i*3)%10==1&&i<parm_data->NBONH-1)
      fprintf(fptr, "\n");
  }
  
  // Bonds without hydrogen
  fprintf(fptr, "\n%%FLAG BONDS_WITHOUT_HYDROGEN\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  for (i=0; i<parm_data->NBONA; ++i)
  {
    fprintf(fptr, "%8i%8i%8i",
            parm_data->pbond[i].ib, parm_data->pbond[i].jb, parm_data->pbond[i].icb);
    if ((i*3)%10==1&&i<parm_data->NBONA-1)
      fprintf(fptr, "\n");
  }
  
  // Angles with hydrogen
  fprintf(fptr, "\n%%FLAG ANGLES_INC_HYDROGEN\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  for (i=0; i<parm_data->NTHETH; ++i)
  {
    fprintf(fptr, "%8i%8i%8i%8i",
            parm_data->pangleH[i].it, parm_data->pangleH[i].jt, parm_data->pangleH[i].kt, parm_data->pangleH[i].ict);
    if ((i*4)%10==1&&i<parm_data->NTHETH-1)
      fprintf(fptr, "\n");
  }
  
  // Angles without hydrogen
  fprintf(fptr, "\n%%FLAG ANGLES_WITHOUT_HYDROGEN\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  for (i=0; i<parm_data->MTHETS; ++i)
  {
    fprintf(fptr, "%8i%8i%8i%8i",
            parm_data->pangle[i].it, parm_data->pangle[i].jt, parm_data->pangle[i].kt, parm_data->pangle[i].ict);
    if ((i*4)%10==1&&i<parm_data->MTHETS-1)
      fprintf(fptr, "\n");
  }
  
  // Dihedrals with hydrogen
  fprintf(fptr, "\n%%FLAG DIHEDRALS_INC_HYDROGEN\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  for (i=0; i<parm_data->NPHIH; ++i)
  {
    fprintf(fptr, "%8i%8i%8i%8i%8i",
            parm_data->pdihedralH[i].ip, parm_data->pdihedralH[i].jp, parm_data->pdihedralH[i].kp,
            parm_data->pdihedralH[i].lp, parm_data->pdihedralH[i].icp);
    if ((i*5)%10==1&&i<parm_data->NPHIH-1)
      fprintf(fptr, "\n");
  }
  
  // Dihedrals without hydrogen
  fprintf(fptr, "\n%%FLAG DIHEDRALS_WITHOUT_HYDROGEN\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  for (i=0; i<parm_data->MPHIA; ++i)
  {
    fprintf(fptr, "%8i%8i%8i%8i%8i",
            parm_data->pdihedral[i].ip, parm_data->pdihedral[i].jp, parm_data->pdihedral[i].kp,
            parm_data->pdihedral[i].lp, parm_data->pdihedral[i].icp);
    if ((i*5)%10==1&&i<parm_data->MPHIA-1)
      fprintf(fptr, "\n");
  }
  
  // Excluded atom list
  fprintf(fptr, "\n%%FLAG EXCLUDED_ATOMS_LIST\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  for (i=0; i<parm_data->NEXT; ++i)
  {
    fprintf(fptr, "%8i", parm_data->natex[i]);
    if (i%10==1&&i<parm_data->NEXT-1)
      fprintf(fptr, "\n");
  }
  
  // Hydrogen bonding
  fprintf(fptr, "\n%%FLAG HBOND_ACOEF\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<parm_data->NHB; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->ag[i]);
    if (i%5==1&&i<parm_data->NHB-1)
      fprintf(fptr, "\n");
  }
  fprintf(fptr, "%%HBOND_BCOEF\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<parm_data->NHB; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->bg[i]);
    if (i%5==1&&i<parm_data->NHB-1)
      fprintf(fptr, "\n");
  }
  
  // Hbcut even though not in use :P
  fprintf(fptr, "\n%%FLAG HBCUT\n");
  fprintf(fptr, "%%FORMAT(5E16.8)\n");
  for (i=0; i<parm_data->NHB; ++i)
  {
    fprintf(fptr, "%16.8E", parm_data->hbcut[i]);
    if (i%5==1&&i<parm_data->NHB-1)
      fprintf(fptr, "\n");
  }
  
  // AMBER atom types
  fprintf(fptr, "\n%%FLAG AMBER_ATOM_TYPE\n");
  fprintf(fptr, "%%FORMAT(20A4)\n");
  for (i=0; i<parm_data->NTOTAT; ++i)
  {
    fprintf(fptr, "%4s", parm_data->atom[i].isymbl);
    if (i%20==1&&i<parm_data->NTOTAT-1)
      fprintf(fptr, "\n");
  }
  
  // Tree chain
  fprintf(fptr, "\n%%FLAG TREE_CHAIN_CLASSIFICATION\n");
  fprintf(fptr, "%%FORMAT(20a4)\n");
  for (i=0; i<parm_data->NTOTAT; ++i)
  {
    fprintf(fptr, "%4s", parm_data->atom[i].itree);
    if (i%20==1&&i<parm_data->NTOTAT-1)
      fprintf(fptr, "\n");
  }
  
  // Join array
  fprintf(fptr, "\n%%FLAG JOIN_ARRAY\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  for (i=0; i<parm_data->NTOTAT; ++i)
  {
    fprintf(fptr, "%8i", parm_data->atom[i].join);
    if (i%10==1&&i<parm_data->NTOTAT-1)
      fprintf(fptr, "\n");
  }
  
  // Rotation (not useful for sander)
  fprintf(fptr, "\n%%FLAG IROTAT\n");
  fprintf(fptr, "%%FORMAT(10I8)\n");
  for (i=0; i<parm_data->NTOTAT; ++i)
  {
    fprintf(fptr, "%8i", parm_data->atom[i].irotat);
    if (i%10==1&&i<parm_data->NTOTAT-1)
      fprintf(fptr, "\n");
  }
  
  // That's all the important stuff from the prmtop.
  return SUCCESS;
}

int write_mdcrd(global_options_struct *global_options, coords_struct *coords_data)
{
  FILE *fptr;
  int structure;
  int i;
  
  if ((fptr=fopen(global_options->mdcrd_filename, "w"))==NULL)
  {
    file_open_failure("write_mdcrd", global_options->mdcrd_filename);
    return FILE_OPEN_FAIL;
  }
  fprintf(fptr, "\n");
  
  for (structure=0; structure<global_options->NSTRUCTURES; ++structure)
  {
    if (coords_data[structure].energy != 0.0)
    {
      for (i=0; i<global_options->NATOMS; ++i)
      {
        fprintf(fptr, "%f %f %f ", coords_data[structure].x_coord[i], coords_data[structure].y_coord[i], coords_data[structure].z_coord[i]);
      }
      fprintf(fptr, "\n");
    }
  }
  fclose(fptr);
  
  if ((fptr=fopen(global_options->energy_filename, "w"))==NULL)
  {
    file_open_failure("write_mdcrd", global_options->energy_filename);
    return FILE_OPEN_FAIL;
  }
  
  for (structure=0; structure<global_options->NSTRUCTURES; ++structure)
  {
    if (coords_data[structure].energy != 0.0)
      fprintf(fptr, "%f\n", coords_data[structure].energy);
  }
  return SUCCESS;
}


