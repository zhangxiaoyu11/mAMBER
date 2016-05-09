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

/*param_summary.c*/

/*This routine prints out a summary of all of the parameters currently stored in parm_data*/

#include <stdio.h>
#include "function_def.h"
#include "constants.h"

void print_parameter_summary(global_options_struct *global_options, parm_struct *parm_data)
{
  int param;
  
  if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD || 
      global_options->FUNC_TO_FIT==AMBER_FORCES ||
      global_options->FUNC_TO_FIT==DIHEDRAL_LEAST_SQUARES )
  {
     printf("parameters for force field equation: AMBER_STANDARD:\n");
     if (global_options->VERBOSITY >= HIGH)
      printf("   (* means parameter is NOT constant during fit)\n");
     else
       printf("   (* means parameter is NOT constant during fit)\n");
     

     /*First off, if K was fitted print the value of K*/
     if (global_options->K_FIT==YES)
        printf("                         *K = %8.6f kcal/mol\n",global_options->K);
     else
       printf("                         K = %8.6f kcal/mol\n", global_options->K);

     /*next, print the bond parameters*/
     for (param=0;param<parm_data->unique_bonds_found;++param)
     {
        printf("             (%s-%s)",parm_data->bond_data[param].atom_type1,parm_data->bond_data[param].atom_type2);
        if (parm_data->bond_data[param].DO_FIT_KR==YES)
        {
          if (global_options->VERBOSITY>=HIGH)
            printf("*");
          else
            printf("*");
        }
        else
           printf(" ");
        printf("Kr = %8.4f kcal/(mol A)^2,",parm_data->bond_data[param].rk);
        if (parm_data->bond_data[param].DO_FIT_REQ==YES)
        {
          if (global_options->VERBOSITY>=HIGH)
            printf("*");
          else
            printf("*");
        }
        else
           printf(" ");
        printf("r_eq = %8.4f A \n",parm_data->bond_data[param].req);
     }
     /*next, print the angle parameters*/
     for (param=0;param<parm_data->unique_angles_found;++param)
     {
        printf("        (%s-%s-%s)",parm_data->angle_data[param].atom_type1,parm_data->angle_data[param].atom_type2,parm_data->angle_data[param].atom_type3);
        if (parm_data->angle_data[param].DO_FIT_KT==YES)
        {
          if (global_options->VERBOSITY>=HIGH)
            printf("*");
          else
            printf("*");
        }
        else
           printf(" ");
        printf("Kt = %8.4f kcal/(mol rad)^2, ",parm_data->angle_data[param].tk);
        if (parm_data->angle_data[param].DO_FIT_THEQ==YES)
        {
          if (global_options->VERBOSITY>=HIGH)
            printf("*");
          else
            printf("*");
        }
        else
           printf(" ");
        printf("th_eq = %8.4f deg \n", parm_data->angle_data[param].teq*RADIAN_TO_DEGREE);
      }
     /*finally, print the dihedral parameters*/
     for (param=0;param<parm_data->unique_dihedrals_found;++param)
     {
				if (parm_data->dihedral_data[param].improper==YES)
					printf("   IMP ");
				else
					printf("       ");
        printf("(%s-%s-%s-%s)",parm_data->dihedral_data[param].atom_type1,parm_data->dihedral_data[param].atom_type2,parm_data->dihedral_data[param].atom_type3,parm_data->dihedral_data[param].atom_type4);
        if (parm_data->dihedral_data[param].DO_FIT_KP==YES)
        {
          if (global_options->VERBOSITY>=HIGH)
            printf("*");
          else
            printf("*");
        }
        else
          printf(" ");
        printf("Kp = %8.4f kcal/mol, ",parm_data->dihedral_data[param].pk);
        if (parm_data->dihedral_data[param].DO_FIT_NP==YES)
        {
          if (global_options->VERBOSITY>=HIGH)
            printf("*");
          else
            printf("*");
        }
        else
           printf(" ");
        printf("Np = %6.4f, ",parm_data->dihedral_data[param].pn);
        if (parm_data->dihedral_data[param].DO_FIT_PHASE==YES)
        {
          if (global_options->VERBOSITY>=HIGH)
            printf("*");
          else
            printf("*");
        }
        else
           printf(" ");
        printf("Phase = %8.4f Deg \n",parm_data->dihedral_data[param].phase*RADIAN_TO_DEGREE);           
     }
   }
  else
  {
     printf("  FORCE FIELD EQUATION IS NOT CURRENTLY IMPLEMENTED.\n");
  }

  
}

