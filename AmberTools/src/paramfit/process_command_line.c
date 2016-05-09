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

/*process_command_line.c*/

/* int process_command_line(int argc, char *cmd_line_options[], struct global_options_struct *global_options);

   This module takes the arguments int argc and char *argv[] from main containing
   the number of command line options and a pointer to an array containing each of
   the options with argb[0] containing the name of the program.

   The function returns the following codes:
      SUCCESS
      UNKNOWN_OPT
      TOO_MANY_OPT
      ALLOC_FAIL 
      CMD_HELP_REQ
      HIST_REQ
      HELP_REQ  

command line options are:

./prog.x -i Job_Control.in -p prmtop -c mdcrd -q QM_data.dat -d ON/OFF/DEBUG

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "function_def.h"

#define MAX_CMDLINE_OPTIONS 13  /* Prog name + 12 options*/

int process_command_line(int argc, char **argv, global_options_struct *global_options)
{
  int i;
  /*Step 1, check to see if too many options have been specified */
  if (argc > MAX_CMDLINE_OPTIONS)
  {
     printf("*** ERROR - Too many command line options specified.\n");
     command_line_help(argv);
     return TOO_MANY_OPT;
  }

/*If the user specified no command line options I will switch on diagnostics
and inform the user that we are using default options*/
  if (argc<2)
  {
      printf("\n!  No command line options given - LOADING DEFAULTS\n");
      printf("!  Setting verbosity to medium.\n\n");
      global_options->VERBOSITY=MEDIUM;
      return SUCCESS;
  }

  if (!(strcmp(argv[1],"/?")) || !(strcmp(argv[1],"-?")) || !(strcmp(argv[1],"/help")) || !(strcmp(argv[1],"-help")) || !(strcmp(argv[1],"--help")))
  {
     /*User has requested help */
     command_line_help(argv);
     return CMD_HELP_REQ;
  }
  if (!(strcmp(argv[1],"/history")) || !(strcmp(argv[1],"/HISTORY")) || !(strcmp(argv[1],"/History")))
  {
     /*User has requested program history */
     print_program_history();
     return HIST_REQ;
  }
  /*step 2, check each command line option to see if it is valid*/
  for (i=1; i<argc; i++)
  {
     /*Command line option 1 = -i Job_Control.in*/
     if(!(strcmp(argv[i],"-i")))
     {
        /*option is for job_control_filename*/
        ++i; /*Increment i to get next command line option */

        if (i >= argc) /*User didn't specify anything after the -i */
        {
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help(argv);
          return UNKNOWN_OPT;
        }
        /* global_options->job_control_filename has not been allocated by default so we alloc it.
           so we have to use realloc */
        global_options->job_control_filename = (char *)malloc( ( strlen(argv[i])+1 ) * sizeof(char) );
        if (global_options->job_control_filename == NULL)
        {
           malloc_failure_char("process_command_line", "global_options->job_control_filename", (strlen(argv[i])+1));
           return ALLOC_FAIL;
        }
        global_options->mem_allocated+=(((strlen(argv[i])+1)*sizeof(char)));
        strcpy(global_options->job_control_filename,argv[i]);
     }
     else if(!(strcmp(argv[i],"-p")))
     {        /*Command line option 2 = -p prmtop*/
        /*option is for prmtop filename*/
         ++i; /*Increment i to get next command line option */

         if (i >= argc) /*User didn't specify anything after the -p */
         {
            printf("*** ERROR - Unknown command line option.\n");
            command_line_help(argv);
            return UNKNOWN_OPT;
         }

         /* global_options->prmtop_filename was already allocated by set_default_options()
            so we have to use realloc */
         global_options->prmtop_filename = (char *)realloc(global_options->prmtop_filename,(strlen(argv[i])+1));
         if (global_options->prmtop_filename == NULL)
         {
            malloc_failure_char("process_command_line", "global_options->prmtop_filename", ((strlen(argv[i])+1)*sizeof(char)));
            return ALLOC_FAIL;
         }
         global_options->mem_allocated+=(((strlen(argv[i])+1)*sizeof(char))-(7*sizeof(char)));
         strcpy(global_options->prmtop_filename,argv[i]);
     }
     else if(!(strcmp(argv[i],"-c")))
     {   /* Option 3 = -c mdcrd */
         /*option is for mdcrd filename*/
         ++i; /*Increment i to get next command line option */

         if (i >= argc) /*User didn't specify anything after the -c */
         {
            printf("*** ERROR - Unknown command line option.\n");
            command_line_help(argv);
            return UNKNOWN_OPT;
         }
         /* global_options->mdcrd_filename was already allocated by set_default_options()
            so we have to use realloc */
         global_options->mdcrd_filename = (char *)realloc(global_options->mdcrd_filename,((strlen(argv[i])+1)*sizeof(char)));
         if (global_options->mdcrd_filename == NULL)
         {
            malloc_failure_char("process_command_line", "global_options->mdcrd_filename", (strlen(argv[i])+1));
            return ALLOC_FAIL;
         }
         global_options->mem_allocated+=(((strlen(argv[i])+1)*sizeof(char))-(6*sizeof(char)));         
         strcpy(global_options->mdcrd_filename,argv[i]);
     }
     else if(!(strcmp(argv[i],"-q")))
     {   /* Option 4 = -q QM_data.dat */
         /*option is for energy filename*/
         ++i; /*Increment i to get next command line option */

         if (i >= argc) /*User didn't specify anything after the -q */
         {
            printf("*** ERROR - Unknown command line option.\n");
            command_line_help(argv);
            return UNKNOWN_OPT;
         }
         /* global_options->energy_filename was already allocated by set_default_options()
            so we have to use realloc */
         global_options->energy_filename = (char *)realloc(global_options->energy_filename,((strlen(argv[i])+1)*sizeof(char)));
         if (global_options->energy_filename == NULL)
         {
            malloc_failure_char("process_command_line", "global_options->energy_filename", (strlen(argv[i])+1));
            return ALLOC_FAIL;
         }
         global_options->mem_allocated+=(((strlen(argv[i])+1)*sizeof(char))-(12*sizeof(char)));         
         strcpy(global_options->energy_filename,argv[i]);
     }
     else if(!(strcmp(argv[i],"-v")))
     { /* Option 5 = -d Verbosity setting */
        ++i;
        if (i >= argc) /*User didn't specify anything after the -v */
        {
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help(argv);
          return UNKNOWN_OPT;
        }             
        /*Valid options are HIGH, MEDIUM, LOW, or lower case versions */
        if ((strcmp(argv[i],"HIGH")) && (strcmp(argv[i],"MEDIUM")) && (strcmp(argv[i],"LOW")) && 
            (strcmp(argv[i],"high")) && (strcmp(argv[i],"medium")) && (strcmp(argv[i],"low")))
        {
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help(argv);
          return UNKNOWN_OPT; /*User specified unknown debug level*/
        }
        /*determine what level the user has specified*/
        if (!(strcmp(argv[i],"LOW"))||!(strcmp(argv[i],"low")))
        {
          global_options->VERBOSITY=LOW;
        }
        else if (!(strcmp(argv[i],"MEDIUM"))||!(strcmp(argv[i],"medium")))
        {
          global_options->VERBOSITY=MEDIUM;
        }
        else if (!(strcmp(argv[i],"HIGH"))||!(strcmp(argv[i],"high")))
        {
          global_options->VERBOSITY=HIGH;
        }
        else
        {
          /*Unknown option*/
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help(argv);
          return UNKNOWN_OPT;
        }
      }
      else if (!(strcmp(argv[i], "--random-seed")))
      {
        ++i; 
        if (i>=argc)
        {
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help(argv);
          return UNKNOWN_OPT;
        }
        if (sscanf(argv[i], "%d", &global_options->RANDOM_SEED) != 1) 
        {
          // Invalid data for the random seet
          printf("*** ERROR - invalid random seed.\n");
          command_line_help(argv);
          return UNKNOWN_OPT;
        }
      }
      else if (!(strcmp(argv[i], "-d")))
      {
        // support for deprecated -d option
        ++i;
        if (i >= argc)
        {
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help(argv);
          return UNKNOWN_OPT;
        }
        if(!(strcmp(argv[i],"ON")))
          global_options->VERBOSITY=MEDIUM;
        else if (!(strcmp(argv[i],"OFF")))
          global_options->VERBOSITY=LOW;
        else if (!(strcmp(argv[i],"DEBUG")))
          global_options->VERBOSITY=HIGH;
        else
        {
          printf("*** ERROR - Unknown command line option.\n");
          command_line_help(argv);
          return UNKNOWN_OPT;
        }
        printf("!  WARNING - The -d command line parameter is deprecated. Use the -v option instead.\n");
      }
      else
      {
        printf("*** ERROR - Unknown command line option.\n");
        command_line_help(argv);
        return UNKNOWN_OPT; /*Unknown command line option*/
      }
  }
  return SUCCESS;
}




