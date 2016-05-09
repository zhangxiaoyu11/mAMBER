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
/*This module contains a number of Help routines for printing out help messages*/

#include <stdio.h>

#include "function_def.h"

/*This module prints out the command line help when the user either requests
  help with /? -? /help or -help or if the user provides unknown command line
  arguments */

void command_line_help(char *cmd_line_options[])
{
  printf("Usage is:\n");
  printf("  %s -i Job_Control.in -p prmtop -c mdcrd -q QM_data.dat -v [LOW/MEDIUM/HIGH]\n",cmd_line_options[0]);
  printf("\n");
  printf("All switches are optional. If not specified the following defaults are loaded:\n");
  printf("     -i Job_Control.in\n");
  printf("     -p prmtop\n");
  printf("     -c mdcrd\n");
  printf("     -q QM_data.dat\n");
  printf("     -v MEDIUM\n");
  printf("     --random-seed 0 (for debugging only, no default value)\n\n");
  printf("     /history prints program development history\n\n");
  printf("     For HELP please see the documentation\n");
  printf("\n");
}
