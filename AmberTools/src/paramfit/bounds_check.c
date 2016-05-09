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

/*bounds_check.c*/

/*

Contains functions that ensure the minimisation algorithms
have not wandered off into areas where there are too little
data in the input structures.

int calculate_structure_diversity
int check_bounds
int check_bonds
int check_angles
int check_dihedrals
  
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"

#include "function_def.h"

int calculate_structure_diversity(global_options_struct *global_options, bounds_struct *bounds_data, parm_struct *parm_data, coords_struct *coords_data)
{
  /* Allocate necessary memory for the arrays */
  if (global_options->VERBOSITY>=HIGH)
    printf("Allocating %d bytes for *bond_lengths.\n",(int)((global_options->NSTRUCTURES)*(parm_data->unique_bonds_found)*sizeof(double)));
  bounds_data->bond_lengths=alloc_2D_double(parm_data->unique_bonds_found, global_options->NSTRUCTURES);
  if (bounds_data->bond_lengths==NULL)
  {
    malloc_failure_double("calculate_structure_diversity", "bond_lengths", (global_options->NSTRUCTURES));
    return ALLOC_FAIL;
  }
  if (global_options->VERBOSITY>=HIGH)
    printf("Allocating %d bytes for *angle_thetas.\n",(int)((global_options->NSTRUCTURES)*(parm_data->unique_angles_found)*sizeof(double)));
  bounds_data->angle_thetas=alloc_2D_double(parm_data->unique_angles_found, global_options->NSTRUCTURES);
  if (bounds_data->angle_thetas==NULL)
  {
    malloc_failure_double("calculate_structure_diversity", "angle_thetas", (global_options->NSTRUCTURES));
    return ALLOC_FAIL;
  }
  if (global_options->VERBOSITY>=HIGH)    
    printf("Allocating %d bytes for *dihedral_thetas.\n",(int)((global_options->NSTRUCTURES)*(parm_data->unique_dihedrals_found)*sizeof(double)));
  bounds_data->dihedral_thetas=alloc_2D_double(parm_data->unique_dihedrals_found, global_options->NSTRUCTURES);
  if (bounds_data->dihedral_thetas==NULL)
  {
    malloc_failure_double("calculate_structure_diversity", "dihedral_thetas", (global_options->NSTRUCTURES));
    return ALLOC_FAIL;
  }
  bounds_data->mem_allocated+=(global_options->NSTRUCTURES*3*sizeof(double));
  
  /* Calculate the actual available parameters and put them all in the arrays */
  int structure=0;
  int i, j;
  double tempx1, tempy1, tempz1, tempx2, tempy2, tempz2;
  double tempx3, tempy3, tempz3, tempx4, tempy4, tempz4;  
  
  if (global_options->VERBOSITY>=HIGH)
  {
    printf("Gathering structure information\n");
    printf("There are %d unique angles and %d structs\n", parm_data->unique_angles_found, global_options->NSTRUCTURES);
  }
  
    for (structure=0; structure < global_options->NSTRUCTURES; ++structure)
    {
      /* Calculate available bond lengths and put them in the arrays */
      for (i=0;i<parm_data->unique_bonds_found;++i)
      {
        for (j=0;j<parm_data->bond_data[i].number;++j)
        {
          tempx1=coords_data[structure].x_coord[parm_data->bond_data[i].atom1[j]-1];
          tempy1=coords_data[structure].y_coord[parm_data->bond_data[i].atom1[j]-1];
          tempz1=coords_data[structure].z_coord[parm_data->bond_data[i].atom1[j]-1];
          tempx2=coords_data[structure].x_coord[parm_data->bond_data[i].atom2[j]-1];
          tempy2=coords_data[structure].y_coord[parm_data->bond_data[i].atom2[j]-1];
          tempz2=coords_data[structure].z_coord[parm_data->bond_data[i].atom2[j]-1];
          
          bounds_data->bond_lengths[i][structure] += calc_bond_length(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2);
        }
        // average bond length for this type of bond since we only optimize types
        bounds_data->bond_lengths[i][structure] /= parm_data->bond_data[i].number;
      }
    
      /* Calculate available bond angles and put them in array */
      for (i=0;i<parm_data->unique_angles_found;++i)
      {
        for (j=0;j<parm_data->angle_data[i].number;++j)
        {
          tempx1=coords_data[structure].x_coord[parm_data->angle_data[i].atom1[j]-1];
          tempy1=coords_data[structure].y_coord[parm_data->angle_data[i].atom1[j]-1];
          tempz1=coords_data[structure].z_coord[parm_data->angle_data[i].atom1[j]-1];
          tempx2=coords_data[structure].x_coord[parm_data->angle_data[i].atom2[j]-1];
          tempy2=coords_data[structure].y_coord[parm_data->angle_data[i].atom2[j]-1];
          tempz2=coords_data[structure].z_coord[parm_data->angle_data[i].atom2[j]-1];
          tempx3=coords_data[structure].x_coord[parm_data->angle_data[i].atom3[j]-1];
          tempy3=coords_data[structure].y_coord[parm_data->angle_data[i].atom3[j]-1];
          tempz3=coords_data[structure].z_coord[parm_data->angle_data[i].atom3[j]-1];

          /*find angle*/
          double rad = calc_angle_radians(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2,tempx3,tempy3,tempz3);
          if (rad != rad) // check for NaN error
          {
            printf("ERROR: NaN error:\n");
            printf("Unique angle %i, structure %i, mult %i\n", i, structure, j);
            printf("(%f,%f,%f) (%f,%f,%f) (%f,%f,%f) = %f\n", tempx1,tempy1,tempz1,tempx2,tempy2,tempz2,tempx3,tempy3,tempz3, rad);
          }
          bounds_data->angle_thetas[i][structure] += rad;
          // Fix cases where angle is out of range 0 to PI
          if (bounds_data->angle_thetas[i][structure] >= PI)
            bounds_data->angle_thetas[i][structure] -= 2*PI;
          if (bounds_data->angle_thetas[i][structure] < 0.0)
            bounds_data->angle_thetas[i][structure] += 2*PI;
        }
        bounds_data->angle_thetas[i][structure] /= (double)parm_data->angle_data[i].number;
      }
      
      /* Calculate available dihedral angles and put them in their array */
      /*Dihedrals*/
      for (i=0;i<parm_data->unique_dihedrals_found;++i)
      {
        for (j=0;j<parm_data->dihedral_data[i].number;++j)
        {
          /*First off we need to find the dihedral between the atoms in the angle*/
          /*Get the coordinates*/
          tempx1=coords_data[structure].x_coord[parm_data->dihedral_data[i].atom1[j]-1];
          tempy1=coords_data[structure].y_coord[parm_data->dihedral_data[i].atom1[j]-1];
          tempz1=coords_data[structure].z_coord[parm_data->dihedral_data[i].atom1[j]-1];
          tempx2=coords_data[structure].x_coord[parm_data->dihedral_data[i].atom2[j]-1];
          tempy2=coords_data[structure].y_coord[parm_data->dihedral_data[i].atom2[j]-1];
          tempz2=coords_data[structure].z_coord[parm_data->dihedral_data[i].atom2[j]-1];
          tempx3=coords_data[structure].x_coord[parm_data->dihedral_data[i].atom3[j]-1];
          tempy3=coords_data[structure].y_coord[parm_data->dihedral_data[i].atom3[j]-1];
          tempz3=coords_data[structure].z_coord[parm_data->dihedral_data[i].atom3[j]-1];
          tempx4=coords_data[structure].x_coord[parm_data->dihedral_data[i].atom4[j]-1];
          tempy4=coords_data[structure].y_coord[parm_data->dihedral_data[i].atom4[j]-1];
          tempz4=coords_data[structure].z_coord[parm_data->dihedral_data[i].atom4[j]-1];
          
          /*find dihedral*/
          bounds_data->dihedral_thetas[i][structure] += (PI-calc_dihedral_radians(tempx1,tempy1,tempz1,tempx2,tempy2,tempz2,tempx3,tempy3,tempz3,tempx4,tempy4,tempz4));

          // fix out of bounds
          while (bounds_data->dihedral_thetas[i][structure] >= PI)
            bounds_data->dihedral_thetas[i][structure] -= 2*PI;
          while (bounds_data->dihedral_thetas[i][structure] < -1*PI)
            bounds_data->dihedral_thetas[i][structure] += 2*PI;
        }
        bounds_data->dihedral_thetas[i][structure] /= (double)parm_data->dihedral_data[i].number;
      }
    }

  /* Now sort each of the arrays in ascending order */
  if (global_options->VERBOSITY>=HIGH)
    printf("Sorting structure information\n");
  /* Sort each type of bond */
  for (i=0; i<parm_data->unique_bonds_found; ++i)
  {
    for (structure=1; structure<global_options->NSTRUCTURES; ++structure)
    {
      double temp = bounds_data->bond_lengths[i][structure];
      
      j = structure-1;
      while(j >=0 && bounds_data->bond_lengths[i][j] > temp)
      {
        bounds_data->bond_lengths[i][j+1] = bounds_data->bond_lengths[i][j];
        --j;
      }
      bounds_data->bond_lengths[i][j+1] = temp;
    }
  }
  /* Now sort angles */
  for (i=0; i<parm_data->unique_angles_found; ++i)
  {
    for (structure=1; structure<global_options->NSTRUCTURES; ++structure)
    {
      double temp = bounds_data->angle_thetas[i][structure];
      
      j = structure-1;
      while(j >=0 && bounds_data->angle_thetas[i][j] > temp)
      {
        bounds_data->angle_thetas[i][j+1] = bounds_data->angle_thetas[i][j];
        --j;
      }
      bounds_data->angle_thetas[i][j+1] = temp;
    }
  }
  /* Sort the dihedrals */
  for (i=0; i<parm_data->unique_dihedrals_found; ++i)
  {
    for (structure=1; structure<global_options->NSTRUCTURES; ++structure)
    {
      double temp = bounds_data->dihedral_thetas[i][structure];
      
      j = structure-1;
      while(j >=0 && bounds_data->dihedral_thetas[i][j] > temp)
      {
        bounds_data->dihedral_thetas[i][j+1] = bounds_data->dihedral_thetas[i][j];
        --j;
      }
      bounds_data->dihedral_thetas[i][j+1] = temp;
    }
  }
  
 return SUCCESS;
}

void clean_up_bounds(bounds_struct *bounds_data)
{
  double_2D_array_free(bounds_data->bond_lengths);
  double_2D_array_free(bounds_data->angle_thetas);
  double_2D_array_free(bounds_data->dihedral_thetas);
}

int check_dihedrals(global_options_struct *global_options, parm_struct *parm_data, bounds_struct *bounds_data)
{
  double increment = PI/global_options->DIHEDRAL_SPAN;     // there should be a dihedral value every this many rads
  double rad=0;
  double invalid_min=-1.0, invalid_max=-1.0;               // to give information about where data is missing.
  int dihedral=0;
  int num_missing=0;
  short valid=YES;
  int theta = 0;        // keeps track of unique dihedrals
  
  for (theta=0; theta<parm_data->unique_dihedrals_found; ++theta) // check each dihedral is well represented
  {
    if (parm_data->dihedral_data[theta].DO_FIT_PHASE==YES || parm_data->dihedral_data[theta].DO_FIT_KP==YES ) // only check ones that are being fitted
    {
      // Adjust phase to be between 0 and 2PI if necessary
      if (parm_data->dihedral_data[theta].phase >= 2.0*PI)
        parm_data->dihedral_data[theta].phase = fmod(parm_data->dihedral_data[theta].phase, 2.0*PI);
      else if (parm_data->dihedral_data[theta].phase < 0.0)
        parm_data->dihedral_data[theta].phase = 2.0*PI - fmod(parm_data->dihedral_data[theta].phase, 2.0*PI);
      
      // Create scatter plot file if desired
      if (global_options->SCATTERPLOTS==TRUE)
  {
        FILE *plot;
        char filename[10];
        sprintf(filename, "%d.diheq", theta);
        printf("  Writing dihedral scatter plot %s\n", filename);
        plot = fopen(filename, "w");
        fprintf(plot, "#%s-%s-%s-%s\n", parm_data->dihedral_data[theta].atom_type1, parm_data->dihedral_data[theta].atom_type2,
                parm_data->dihedral_data[theta].atom_type3, parm_data->dihedral_data[theta].atom_type4);
        for (dihedral=0; dihedral<global_options->NSTRUCTURES; ++dihedral) {
          fprintf(plot, "%f  1\n", bounds_data->dihedral_thetas[theta][dihedral]);
        }
        fclose(plot);
      }
      rad = 0;
      dihedral = 0;
      valid = YES;
      num_missing = 0;
      
      while (rad+increment <= PI)
      {
        dihedral = 0;
        valid = YES;
        int done = NO;
        // attempt to find a dihedral with value in this range
        while ( dihedral < global_options->NSTRUCTURES && done==NO )
        {
          if ( (bounds_data->dihedral_thetas[theta][dihedral]> rad) && (bounds_data->dihedral_thetas[theta][dihedral] < rad+increment) )
            done = YES;
          else
            ++dihedral;
        }
        
        // check if a structure could be found
        if (dihedral==global_options->NSTRUCTURES)
        {
          if (invalid_min==-1.0)
            invalid_min = rad;
          invalid_max = rad+increment;
          ++num_missing;
          valid = NO;
          
          if (global_options->VERBOSITY>=HIGH)
            printf("%s-%s-%s-%s dihedral FAILED span check.  No data in range %.4f - %.4f rads.\n",  parm_data->dihedral_data[theta].atom_type1,
                                                                                                        parm_data->dihedral_data[theta].atom_type2,
                                                                                                        parm_data->dihedral_data[theta].atom_type3,
                                                                                                        parm_data->dihedral_data[theta].atom_type4,
                                                                                                        rad, rad+increment);
        }
        rad+=increment;
      }
      
      // print out results if it passed
      if (valid==YES && global_options->VERBOSITY>=HIGH)
        printf("%s-%s-%s-%s dihedral PASSED span check for range %.4f - %.4f rads.\n",  parm_data->dihedral_data[theta].atom_type1,
               parm_data->dihedral_data[theta].atom_type2,
               parm_data->dihedral_data[theta].atom_type3,
               parm_data->dihedral_data[theta].atom_type4,
               rad, rad+increment);
               
      if (valid==NO && global_options->VERBOSITY>=MEDIUM)
      {
        if (global_options->CHECK_BOUNDS==YES)
          printf("   ERROR: ");
        else
          printf("   WARNING: ");
        printf("%s-%s-%s-%s dihedral is missing %.i data points in the range %.4f to %.4f radians.\n",  parm_data->dihedral_data[theta].atom_type1,
                                                                                                           parm_data->dihedral_data[theta].atom_type2,
                                                                                                           parm_data->dihedral_data[theta].atom_type3,
                                                                                                           parm_data->dihedral_data[theta].atom_type4,
                                                                                                           num_missing, invalid_min, invalid_max);
      }
    }
  }
  
  if (valid==YES)
  {
    if (global_options->VERBOSITY >= MEDIUM)
      printf(" * Input structures passed dihedral span check.\n");
    return SUCCESS;
  }
  else
  {
    printf("\n\n");
    if (global_options->CHECK_BOUNDS==YES)
      printf("   ERROR: ");
    else
      printf("   WARNING: ");
    printf("Insufficient dihedral information in sample structures.\n");
    printf("            Your settings require at least %i samples with data    \n", global_options->DIHEDRAL_SPAN);
    printf("            at least every %.3f radians (%.2f degrees).            \n", increment, increment*RADIAN_TO_DEGREE);
    printf("            Either 1) Add the missing input data or                \n");
    printf("                   2) Set DIHEDRAL_SPAN to a smaller value or      \n");
    printf("                   3) Set BOUNDS_CHECK to warn (not recommended).  \n");
    printf("            Please read the help and/or documentation.             \n");
    printf("\n");
    fflush(stdout);
    if (global_options->CHECK_BOUNDS==YES)
      return FAILURE;
    else
      return SUCCESS;
  }
}

int check_angles(global_options_struct *global_options, parm_struct *parm_data, bounds_struct *bounds_data)
{
  // ensure that there are angles around the existing minimum
  int param, pos;
  double diff_next, diff_prev;
  short passed=YES;
  
  for (param=0; param<parm_data->unique_angles_found; ++param)
  {
    // only check angles that are being optimised
    if (parm_data->angle_data[param].DO_FIT_THEQ==YES)
    {
      // Create scatter plot file if desired
      FILE *plot;
      if (global_options->SCATTERPLOTS==TRUE) {
        char filename[10];
        sprintf(filename, "%d.angleq", param);
        printf("  Writing angle scatter plot %s\n", filename);
        plot = fopen(filename, "w");
        fprintf(plot, "#%s-%s-%s\n", parm_data->angle_data[param].atom_type1, parm_data->angle_data[param].atom_type2, parm_data->angle_data[param].atom_type3);
        for (pos=0; pos<global_options->NSTRUCTURES; ++pos) {
          fprintf(plot, "%f  1\n", bounds_data->angle_thetas[param][pos]);
        }
        fclose(plot);
      }
      
      // adjust phase to be between 0 and PI if necessary
      if (parm_data->angle_data[param].teq > PI)
        parm_data->angle_data[param].teq = fmod(parm_data->angle_data[param].teq, PI);
      else if (parm_data->angle_data[param].teq < 0.0)
        parm_data->angle_data[param].teq = PI - fmod(parm_data->angle_data[param].teq, PI);
      
      // find position in parameter array
      pos=0;
      while (pos < global_options->NSTRUCTURES && parm_data->angle_data[param].teq > bounds_data->angle_thetas[param][pos])
        ++pos;
      
      // find distance to neighbors in array
      if (pos==0)
      {
        diff_next = bounds_data->angle_thetas[param][pos] - parm_data->angle_data[param].teq;
        diff_prev = -1.0;
      }
      else if (pos==global_options->NSTRUCTURES)
      {
        diff_next = -1.0;
        diff_prev = parm_data->angle_data[param].teq - bounds_data->angle_thetas[param][pos-1];
      }
      else
      {
        diff_next = bounds_data->angle_thetas[param][pos] - parm_data->angle_data[param].teq;
        diff_prev = parm_data->angle_data[param].teq - bounds_data->angle_thetas[param][pos-1];
      }
      
      if (diff_next > global_options->ANGLE_LIMIT || diff_prev > global_options->ANGLE_LIMIT || diff_next < 0.0 || diff_prev < 0.0)
      {
        passed = NO;
        if (global_options->VERBOSITY>=MEDIUM)
        {
          if (global_options->CHECK_BOUNDS==YES)
            printf("   ERROR: ");
          else
            printf("   WARNING: ");
          if (diff_next == -1.0)
            printf("%s-%s-%s angle has no sample structures after it, nearest sample is %.4f radians before.\n",  
                   parm_data->angle_data[param].atom_type1,
                   parm_data->angle_data[param].atom_type2,
                   parm_data->angle_data[param].atom_type3,
                   diff_prev);
          else if (diff_prev == -1.0)
            printf("%s-%s-%s angle has no sample structures before it, nearest sample is %.4f radians after.\n", 
                   parm_data->angle_data[param].atom_type1,
                   parm_data->angle_data[param].atom_type2,
                   parm_data->angle_data[param].atom_type3,
                   diff_next);
          else
            printf("%s-%s-%s angle has no sample structures within %.4f radians, nearest are %.4f and %.4f radians away.\n",  
                   parm_data->angle_data[param].atom_type1,
                   parm_data->angle_data[param].atom_type2,
                   parm_data->angle_data[param].atom_type3,
                   global_options->ANGLE_LIMIT, diff_prev, diff_next);
        }
      }
      else
        if (global_options->VERBOSITY>=HIGH)
          printf("%s-%s-%s angle PASSED: data exist within %.4f and %.4f radians, passes limit of %.4f.\n",  parm_data->angle_data[param].atom_type1,
                 parm_data->angle_data[param].atom_type2,
                 parm_data->angle_data[param].atom_type3,
                 diff_prev, diff_next, global_options->ANGLE_LIMIT);
    }
  }
  
  if (passed==YES)
  {
    if (global_options->VERBOSITY >= MEDIUM)
      printf(" * Result passed angle validity check.\n");
    return SUCCESS;
  }
  else
  {
    printf("\n\n");
    if (global_options->CHECK_BOUNDS==YES)
      printf("   ERROR: ");
    else
      printf("   WARNING: ");
    printf("Insufficient angle information in sample structures.   \n");
    printf("            Your settings require that sample data exist within    \n");
    printf("            %.3f radians of the algorithm's result.                \n", global_options->ANGLE_LIMIT);
    printf("            Either 1) Add the missing input data                 or\n");
    printf("                   2) Set ANGLE_LIMIT to a larger value          or\n");
    printf("                   3) Set BOUNDS_CHECK to warn (not recommended).  \n");
    printf("            Please read the help and/or documentation.             \n");
    printf("\n");
    fflush(stdout);
    if (global_options->CHECK_BOUNDS==YES)
      return FAILURE;
    else
      return SUCCESS;
  }
  
}

int check_bonds(global_options_struct *global_options, parm_struct *parm_data, bounds_struct *bounds_data)
{
  // ensure that there are bonds around the existing minimum
  int param, pos;
  double diff_next, diff_prev;
  short passed=YES;
  FILE *plot;
  
  for (param=0; param<parm_data->unique_bonds_found; ++param)
  {
    // only check bonds that are being optimised
    if (parm_data->bond_data[param].DO_FIT_REQ==YES)
    {
      // Create scatter plot file if desired
      if (global_options->SCATTERPLOTS==TRUE) {
        char filename[10];
        sprintf(filename, "%d.bondeq\n", param);
        plot = fopen(filename, "w");
        fprintf(plot, filename, "%s-%s", parm_data->bond_data[param].atom_type1, parm_data->bond_data[param].atom_type2);
        for (pos=0; pos<global_options->NSTRUCTURES; ++pos) {
          fprintf(plot, "%f\n", bounds_data->bond_lengths[param][pos]);
        }
        fclose(plot);
      }
      
      // find position in parameter array
      pos=0;
      while (pos < global_options->NSTRUCTURES && parm_data->bond_data[param].req > bounds_data->bond_lengths[param][pos])
        ++pos;
      // find distance to neighbors in array
      if (pos==0)
      {
        diff_next = bounds_data->bond_lengths[param][pos] - parm_data->bond_data[param].req;
        diff_prev = -1.0;
      }
      else if (pos==global_options->NSTRUCTURES)
      {
        diff_next = -1.0;
        diff_prev = parm_data->bond_data[param].req - bounds_data->bond_lengths[param][pos-1];
      }
      else
      {
        diff_next = bounds_data->bond_lengths[param][pos] - parm_data->bond_data[param].req;
        diff_prev = parm_data->bond_data[param].req - bounds_data->bond_lengths[param][pos-1];
      }
      
      if (diff_next > global_options->BOND_LIMIT || diff_prev > global_options->BOND_LIMIT || diff_next < 0.0 || diff_prev < 0.0)
      {
        passed = NO;
        if (global_options->VERBOSITY>=MEDIUM)
        {
          if (global_options->CHECK_BOUNDS==YES)
            printf("   ERROR: ");
          else
            printf("   WARNING: ");
          
          if (diff_next == -1.0)
            printf("%s-%s bond has no sample structures larger than it, nearest sample is %.4f A smaller.\n", parm_data->bond_data[param].atom_type1,
                                                                                                                          parm_data->bond_data[param].atom_type2,
                                                                                                                          diff_prev);
          else if (diff_prev == -1.0)
            printf("%s-%s bond has no sample structures smaller than it, nearest sample is %.4f A larger.\n", parm_data->bond_data[param].atom_type1,
                                                                                                                          parm_data->bond_data[param].atom_type2, 
                                                                                                                          diff_next);
          else
            printf("%s-%s bond has no sample structures within %f A, nearest are %f and %f A\n", parm_data->bond_data[param].atom_type1,
                                                                                                              parm_data->bond_data[param].atom_type2,
                                                                                                              global_options->BOND_LIMIT, diff_prev, diff_next);
        }
        
      }
      else
      {
        if (global_options->VERBOSITY>=HIGH)
          printf("PASSED %s-%s bond: Data exist within %f and %f A, passes limit of %f.\n", parm_data->bond_data[param].atom_type1,
                                                                                            parm_data->bond_data[param].atom_type2,
                                                                                            diff_prev, diff_next, global_options->BOND_LIMIT);
      }
    }
  }
  
  if (passed==YES)
  {
    if (global_options->VERBOSITY >= MEDIUM)
      printf(" * Result passed bond validity check.\n");
    return SUCCESS;
  }
  else
  {
    printf("\n\n");
    if (global_options->CHECK_BOUNDS==YES)
      printf("   mERROR: ");
    else
      printf("   WARNING: ");
    printf("Insufficient bond information in sample structures.    \n");
    printf("            Your settings require that sample data exist within    \n");
    printf("            %.3f A of the algorithm's result.                      \n", global_options->BOND_LIMIT);
    printf("            Either 1) Add the missing input data                 or\n");
    printf("                   2) Set BOND_LIMIT to a larger value           or\n");
    printf("                   3) Set BOUNDS_CHECK to warn (not recommended).  \n");
    printf("            Please read the help and/or documentation.             \n");
    printf("\n");
    fflush(stdout);
    if (global_options->CHECK_BOUNDS==YES)
      return FAILURE;
    else
      return SUCCESS;
  }
} 

/* Checks that the bonds, angles, and dihedrals are in a valid range before
 * conducting a function evaluation. If they are not, it will change any invalid
 * parameter to a random value within the valid range. This prevents the algorithms
 * (especially the simplex algorithm) from crawling into corners of the valid solution
 * space and getting stuck there because there is gradient that points into an 
 * invalid area.
 */
void check_range(global_options_struct *global_options, parm_struct *parm_data)
{
  // Do not bounds check the simplex, or it will mess up.
  if (global_options->ALGORITHM==SIMPLEX) return;
   
  int i;
  for (i=0; i<parm_data->unique_bonds_found; ++i) {
    // Bond force constant between 100 and 1000
    if (parm_data->bond_data[i].rk < 100.0 || parm_data->bond_data[i].rk > 1000.0) {
      parm_data->bond_data[i].rk = 900.0*((double)rand()/(double)RAND_MAX) + 100.0;
    }
    
    // Bond length between 0 and 3 Angstroms
    if (parm_data->bond_data[i].req <= 0.0 || parm_data->bond_data[i].req > 3.0)
      parm_data->bond_data[i].req = 3.0*((double)rand()/((double)RAND_MAX+1));
  }
  for (i=0; i<parm_data->unique_angles_found; ++i) {
    // Angle force constant between 0 and 200
    if (parm_data->angle_data[i].tk < 0.0 || parm_data->angle_data[i].tk > 200.0)
      parm_data->angle_data[i].tk = 200.0*((double)rand()/(double)RAND_MAX);
    
    // Angle should stay between 0 and PI, wrap it around if it's wrong
    if (parm_data->angle_data[i].teq > PI)
        parm_data->angle_data[i].teq = fmod(parm_data->angle_data[i].teq, PI);
    else if (parm_data->angle_data[i].teq < 0.0)
        parm_data->angle_data[i].teq = PI - fmod(parm_data->angle_data[i].teq, PI);
  }
  for (i=0; i<parm_data->unique_dihedrals_found; ++i) {
    // Dihedral force constant between -30 and 30 (Glycam has some negative Kp for example)
    if (parm_data->dihedral_data[i].pk < -30.0 || parm_data->dihedral_data[i].pk > 30.0)
      parm_data->dihedral_data[i].pk = 60*((double)rand()/(double)RAND_MAX) + 30.0;
    
    // Dihedral periodicity between -5 and 5
    if (parm_data->dihedral_data[i].pn < -5.0 || parm_data->dihedral_data[i].pn > 5.0)
      parm_data->dihedral_data[i].pn = 10.0*((double)rand()/(double)RAND_MAX) + 5.0;
    
    // Dihedral phase between 0 and 2*PI, wrap around if wrong
    if (parm_data->dihedral_data[i].phase >= 2.0*PI)
      parm_data->dihedral_data[i].phase = fmod(parm_data->dihedral_data[i].phase, 2.0*PI);
    else if (parm_data->dihedral_data[i].phase < 0.0)
      parm_data->dihedral_data[i].phase = 2.0*PI - fmod(parm_data->dihedral_data[i].phase, 2.0*PI);
  }
}

