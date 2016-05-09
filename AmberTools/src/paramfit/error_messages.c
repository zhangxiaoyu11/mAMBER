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

/*error_messages.c*/
/*This module contains a number of routines for printing out error messages*/

#include <stdio.h>
#include <stdlib.h>

#include "function_def.h"

void process_retval(int err_code, verbosity_t VERBOSITY)
{
  if (err_code != SUCCESS)
  {
    /*Failure or abort flag due to printing of program history / help etc..*/
    /*Lower (less -ve) than ABORT = not a failure, an abort signal*/
    if (VERBOSITY>=HIGH && err_code <= ABORT)
      printf("*** PROGRAM ABORT CODE: %d\n",err_code);
    else if (VERBOSITY>=HIGH)
      printf("*** PROGRAM ERROR CODE: %d\n",err_code);
    
    exit(err_code);
  }  
}

void malloc_failure_char(char *routine, char *var_name, int chars_requested)
{
   /*This routine prints out a standard malloc failure message with the routine
     name, variable name and the number of bytes requested for char data types*/
   printf("\n*** FATAL ERROR IN %s FUNCTION\n",routine);
   printf("m*** COMMAND WAS: %s = malloc(%d * sizeof(char));\n",var_name,chars_requested);
   printf("*** STATUS - FAILED TO ALLOCATE %d bytes\n",(int)(chars_requested*sizeof(char)));
  exit(ALLOC_FAIL);
}

void malloc_failure_int(char *routine, char *var_name, int ints_requested)
{
   /*This routine prints out a standard malloc failure message with the routine
     name, variable name and the number of bytes requested for int data types*/
   printf("\n*** FATAL ERROR IN %s FUNCTION\n",routine);
   printf("*** COMMAND WAS: %s = malloc(%d * sizeof(int));\n",var_name,ints_requested);
   printf("*** STATUS - FAILED TO ALLOCATE %d bytes\n",(int)(ints_requested*sizeof(int)));
   exit(ALLOC_FAIL);
}

void malloc_failure_short_int(char *routine, char *var_name, int short_ints_requested)
{
   /*This routine prints out a standard malloc failure message with the routine
     name, variable name and the number of bytes requested for short int data types*/
   printf("\n*** FATAL ERROR IN %s FUNCTION\n",routine);
   printf("*** COMMAND WAS: %s = malloc(%d * sizeof(short int));\n",var_name,short_ints_requested);
   printf("*** STATUS - FAILED TO ALLOCATE %d bytes\n",(int)(short_ints_requested*sizeof(short int)));
   exit(ALLOC_FAIL);
}


void malloc_failure_double(char *routine, char *var_name, int doubles_requested)
{
   /*This routine prints out a standard malloc failure message with the routine
     name, variable name and the number of bytes requested for char data types*/
   printf("\n*** FATAL ERROR IN %s FUNCTION\n",routine);
   printf("*** COMMAND WAS: %s = malloc(%d * sizeof(double));\n",var_name,doubles_requested);
   printf("*** STATUS - FAILED TO ALLOCATE %d bytes\n",(int)(doubles_requested*sizeof(double)));
   exit(ALLOC_FAIL);
}

void file_open_failure(char *routine, char *var_name)
{
   /*This routine prints out a standard file open failure message with the routine
     name and file name*/
   printf("\n*** FATAL ERROR IN %s FUNCTION\n",routine);
   printf("*** STATUS - FAILED TO OPEN FILE: %s \n",var_name);
   exit(FILE_OPEN_FAIL);
}
