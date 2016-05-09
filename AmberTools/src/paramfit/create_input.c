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

/*
  create_input.c

  This file contains the routines for making the job files from
  a prmtop and mdcrd file that can then be used to obtain the energy
  data to be fitted against
*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "function_def.h"

int create_input(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data)
{
  FILE *fptr;
  int current_struct=0;
  int retval;
  char filename_to_write[1024];
  int temp_len;
  temp_len=0;
  /*Here we will be creating input files based on the structures in the mdcrd file*/

  if (global_options->VERBOSITY>=HIGH)
  {
     printf(" Creating input files. Will write a total of %d sequentially numbered files.\n",global_options->NSTRUCTURES);
     printf("           %6d > .",0);
  }

  /* Read the mdcrd if necessary */
  if (coords_data == NULL)
  {
    retval = read_mdcrd(global_options,coords_data);
    if (retval!=SUCCESS)
    {
      printf("*** ERROR IN create_input SUBROUTINE\n");
      printf("*** FAILED TO READ STRUCTURE %d FROM MDCRD FILE.\n",current_struct);
      return FILE_READ_FAIL;
    }
  }
  
  /* Print out a warning if fitting forces */
  if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
    if (global_options->QMFILEFORMAT==GAUSSIAN)
    printf("!  Will fit forces- ensure the route line in your header has the force keyword!\n");
    else {
      printf("*** ERROR: Force calculations are only supported in Gaussian.\n");
      return INVALID_DATA;
    }
  }
  
  /*loop over all the structures*/
  for (current_struct=0; current_struct<global_options->NSTRUCTURES;++current_struct)
  {
     /*Open the output file with the relevant name*/
     /*filename to write is made up of QMFILEOUTSTART, current_struct and QMFILEOUTEND
     the char variable filename_to_write has been allocated as 1024 bytes long but we should
     still really have a sanity check here to ensure we don't exceed it - not that many people
     would ever want a filename over 1023 characters long*/
     temp_len += strlen(global_options->QMFILEOUTSTART);
     temp_len += strlen(global_options->QMFILEOUTEND);
     ++temp_len;

     /*Note, a potential overflow could occur here - currently NOT checked for*/
     
     /*Now fill the filename_to_write array*/
     sprintf(filename_to_write,"%s%d%s",global_options->QMFILEOUTSTART,current_struct,global_options->QMFILEOUTEND);
     if (global_options->VERBOSITY>=HIGH)
           printf("   Filename to be written is: %s\n",filename_to_write);
     
     if((fptr=fopen(filename_to_write,"w"))==NULL)
     {
        file_open_failure("create_input", filename_to_write);
        return FILE_OPEN_FAIL;
     }

     /*print info for each file we write*/
     if (global_options->VERBOSITY>=HIGH)
     {
       /*check if current_struct is divisible by 50, if it is print a new line*/
       if(current_struct%50==0)
         printf("\n           %6d > .",current_struct);
       else
         printf(".");
       fflush(stdout); /*Flush the printf buffer*/
     }

     /*Now we have the file open call the relevant routines for writing the data*/
     if (global_options->QMFILEFORMAT==GAUSSIAN)
     {
        retval=write_input_gaussian(global_options, parm_data, &coords_data[current_struct], current_struct, fptr);
     }
     else if (global_options->QMFILEFORMAT==ADF)
     {
       retval = write_input_adf(global_options, parm_data, &coords_data[current_struct], current_struct, fptr);
     }
     else if (global_options->QMFILEFORMAT==GAMESS)
     {
       retval = write_input_gamess(global_options, parm_data, &coords_data[current_struct], current_struct, fptr);
     }
     else
     {
       /*Unknown option*/
       printf("   ERROR IN create_input SUBROUTINE\n");
       printf("   UNKNOWN FILE FORMAT TO WRITE: %d.\n",global_options->QMFILEFORMAT);
      return UNKNOWN_OPT;
     }
     if (retval!=SUCCESS)
       return retval;
                                           
     fclose(fptr);
  }

  printf("*   Successfully wrote %i files\n", current_struct-1);
  return SUCCESS;
}

