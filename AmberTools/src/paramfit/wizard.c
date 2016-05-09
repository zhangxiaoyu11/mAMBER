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

/* Defines a series of prompts that can be used to create a job control file
 * Hopefully this will simplify the process of writing a job control file for the
 * user and make paramfit overall easier to use.
 */

#include <stdlib.h>
#include <string.h>
#include "function_def.h"

int get_option(int min, int max);
double get_float();
int job_control_wizard(global_options_struct *global_options);
void genetic_wizard(global_options_struct *global_options, FILE *file);
void simplex_wizard(global_options_struct *global_options, FILE *file);

int job_control_wizard(global_options_struct *global_options)
{
  int input;
  char save[100];
  short write = 0;
  FILE *file;
  
  printf("PARAMFIT SETTINGS WIZARD\n");
  printf("------------------------\n\n");
  printf("Enter a job control filename to save settings as, or hit enter\n> ");
  fgets(save, 100, stdin);
  sscanf(save, "%s", &save);
  if (save[0] != '\n') {
    write = 1;
    printf("Saving job control settings as \"%s\"\n", save);
    file = fopen(save, "w");
  }
  else {
    printf("Not saving job control settings\n");
    file = NULL;
  }
 
  // Number of structures
  printf("\nHow many input structures do you have?\n");
  input = get_option(1, 100000);
  global_options->NSTRUCTURES=input;
  if (write) fprintf(file, "NSTRUCTURES=%d\n", input);
     
  // Runtype
  printf("\nWhat type of job do you want to run?\n\t1) Create input files for quantum program\n\t2) Set which parameters to fit\n\t3) Conduct a fit\n");
  input=get_option(1,3);
  
  if (input == 1) { // CREATE_INPUT
    global_options->RUNTYPE=CREATE_INPUT;
    if (write) fprintf(file, "RUNTYPE=CREATE_INPUT\n");
    
    // Format of quantum input/output files
    printf("\nWhat format should the quantum input files be written in?\n\t1) Gaussian\n\t2) GAMESS\n\t3)ADF\n> ");
    input = get_option(1,3);
    if (input==1) {
      global_options->QMFILEFORMAT=GAUSSIAN;
      if(write) fprintf(file, "QMFILEFORMAT=GAUSSIAN\n");
    }
    else if (input==2) {
      global_options->QMFILEFORMAT=GAMESS;
      if (write) fprintf(file, "QMFILEFORMAT=GAMESS\n");
    }
    else if (input==3) {
      global_options->QMFILEFORMAT=ADF;
      if (write) fprintf(file, "QMFILEFORMAT=ADF\n");
    }
    
    // Header for quantum input/output files
    printf("\nType the location of a header file to be prepended to all quantum output files, or press enter\n> ");
    fgets(save, 100, stdin);
    sscanf(save, "%s", &save);
    while (save[0] != '\n' && check_for_valid_filename(save, strlen(save))!=SUCCESS) {
      printf("\tInvalid response! Enter a valid file name:\n> ");
      fgets(save, 100, stdin);
      sscanf(save, "%s", &save);
    }
    if (save[0] != '\n') {
      global_options->QMHEADER = (char*)malloc( (strlen(save)+1)*sizeof(char) );
      global_options->mem_allocated += ( strlen(save)+1)*sizeof(char);
      strcpy(global_options->QMHEADER, save);
      if (write) fprintf(file, "QMHEADER=%s\n", save);
    }
    
    // Names for quantum input/output files
    printf("\nWhat should the prefix for all quantum input files be? Files will be named PREFIX000SUFFIX\n> ");
    fgets(save, 100, stdin);
    sscanf(save, "%s", &save);
    while (save[0] == '\n' || check_for_valid_filename(save, strlen(save))!=SUCCESS) {
      printf("\tInvalid response! Enter a valid filename prefix:\n> ");
      fgets(save, 100,stdin);
      sscanf(save, "%s", &save);
    }
    global_options->QMFILEOUTSTART=(char*)malloc( (strlen(save)+1)*sizeof(char) );
    global_options->mem_allocated += ( strlen(save)+1 )*sizeof(char);
    strcpy(global_options->QMFILEOUTSTART, save);
    if (write) fprintf(file, "QMFILEOUTSTART=%s\n", save);
    
    printf("\nWhat should the suffix for all quantum input files be? Files will be named %s000SUFFIX\n> ", global_options->QMFILEOUTSTART);
    fgets(save, 100, stdin);
    sscanf(save, "%s", &save);
    while (save[0] == '\n' || check_for_valid_filename(save, strlen(save))!=SUCCESS) {
      printf("\tInvalid response! Enter a valid filename suffix:\n> ");
      fgets(save, 100, stdin);
      sscanf(save, "%s", &save);
    }
    global_options->QMFILEOUTEND=(char*)malloc( (strlen(save)+1)*sizeof(char) );
    global_options->mem_allocated += (strlen(save)+1)*sizeof(char);
    strcpy(global_options->QMFILEOUTEND, save);
    if (write) fprintf(file, "QMFILEOUTEND=%s\n", save);
    
  }
  else if (input == 2) { // SET_PARAMS
    global_options->RUNTYPE=SET_PARAMS;
    if (write) fprintf(file, "RUNTYPE=SET_PARAMS\n");
   
    // Parameter file to save into
    printf("\nEnter a file name to define the parameter list in:\n> ");
    fgets(save, 100, stdin);
    sscanf(save, "%s", &save);
    while (save[0] == '\n' || check_for_valid_filename(save, strlen(save))!=SUCCESS) {
      printf("\tInvalid response! Enter a file name:\n> ");
      fgets(save, 100, stdin);
      sscanf(save, "%s", &save);
    }
    global_options->PARAMETER_FILE_NAME=(char*)malloc( (strlen(save)+1)*sizeof(char) );
    global_options->mem_allocated += ( strlen(save)+1 )*sizeof(char);
    strcpy(global_options->PARAMETER_FILE_NAME, save);
    if (write) fprintf(file, "PARAMETER_FILE_NAME=%s\n", save);
  }
  else if (input == 3) { // FIT
    global_options->RUNTYPE=FIT;
    if (write) fprintf(file, "RUNTYPE=FIT\n");
    
    // Parameters to fit
    printf("\nWhich parameters would you like to fit?\n\t1) Defaults (see manual)\n\t2) Load list from file\n\t3) K only\n> ");
    input=get_option(1,3);
    if (input==1) { // default
      global_options->PARAMETERS_TO_FIT=DEFAULT;
      if (write) fprintf(file, "PARAMETERS_TO_FIT=DEFAULT\n");
    }
    else if (input==2) { // load parameters to fit
      global_options->PARAMETERS_TO_FIT=LOAD;
      if (write) fprintf(file, "PARAMETERS_TO_FIT=LOAD\n");
      // prompt for file name
      printf("Enter file name:\n> ");
      fgets(save, 100, stdin);
      sscanf(save, "%s", &save);
      while (save[0] == '\n' || check_for_valid_filename(save, strlen(save))!=SUCCESS) {
        printf("\tInvalid response! Enter a file name:\n> ");
        fgets(save, 100, stdin);
        sscanf(save, "%s", &save);
      }
      global_options->PARAMETER_FILE_NAME=save;
      if (write) fprintf(file, "PARAMETER_FILE_NAME=%s\n", save);
      // prompt for K
      printf("Enter the initial value for K:\n> "); 
      global_options->K=get_float();
      if (write) fprintf(file, "K=%f\n", global_options->K);
    }
    else if (input==3) { // k only
      global_options->PARAMETERS_TO_FIT=K_ONLY;
      if (write) fprintf(file, "PARAMETERS_TO_FIT=K_ONLY\n");
    }
    
    // Fitting function
   printf("\nWhich function should the fitting be done to?\n\t1) AMBER standard energies\n\t2) AMBER forces on each atom\n\t3) Dihedral least squares fitting\n> ");
   input=get_option(1,3);
   if (input==1) {
     global_options->FUNC_TO_FIT=SUM_SQUARES_AMBER_STANDARD;
     if (write) fprintf(file, "FUNC_TO_FIT=SUM_SQUARES_AMBER_STANDARD\n");
     // QM energy units
     printf("\nWhich units are the energies from the quantum program in?\n\t1) Hartree\n\t2) KCal/mol\n\t3) KJ/mol\n> ");
     input=get_option(1,3);
     if (input==1) {
       global_options->QM_ENERGY_UNITS=HARTREE;
       if (write) fprintf(file, "QM_ENERGY_UNITS=HARTREE\n");
     }
     else if (input==2) {
       global_options->QM_ENERGY_UNITS=KCALMOL;
       if(write) fprintf(file, "QM_ENERGY_UNITS=KCALMOL\n");
     }
     else if (input==3) {
       global_options->QM_ENERGY_UNITS=KJMOL;
       if (write) fprintf(file, "QM_ENERGY_UNITS=KJMOL\n");
     }
   }
   else if (input==2) {
     global_options->FUNC_TO_FIT=AMBER_FORCES;
     if (write) fprintf(file, "FUNC_TO_FIT=AMBER_FORCES\n");
     // QM force units
     printf("\nWhich units are the forces from the quantum program in?\n\t1) Hartree*Bohr\n\t2) KCal/mol*Angstrom\n> ");
     input=get_option(1,2);
     if (input==1) {
      global_options->QM_FORCE_UNITS=HARTREE_BOHR;
      if (write) fprintf(file, "QM_FORCE_UNITS=HARTREE_BOHR\n");
     }
     else if (input==2) {
       global_options->QM_FORCE_UNITS=KCALMOL_ANGSTROM;
       if (write) fprintf(file, "QM_FORCE_UNITS=KCALMOL_ANGSTROM\n");
     }
     printf("\nWhat is the prefix for all quantum output files? Files are expected to be named PREFIX###SUFFIX\n> ");
     fgets(save, 100, stdin);
     sscanf(save, "%s", &save);
     while (save[0] == '\n' || check_for_valid_filename(save, strlen(save))!=SUCCESS) {
       printf("\tInvalid response! Enter a valid filename prefix:\n> ");
     fgets(save, 100,stdin);
     sscanf(save, "%s", &save);
     }
     global_options->QMFILEOUTSTART=(char*)malloc( (strlen(save)+1)*sizeof(char) );
     global_options->mem_allocated += ( strlen(save)+1 )*sizeof(char);
     strcpy(global_options->QMFILEOUTSTART, save);
     if (write) fprintf(file, "QMFILEOUTSTART=%s\n", save);
     
     printf("\nWhat is the suffix for all quantum output files? Files are expected to be named %s###SUFFIX\n> ", global_options->QMFILEOUTSTART);
     fgets(save, 100, stdin);
     sscanf(save, "%s", &save);
     while (save[0] == '\n' || check_for_valid_filename(save, strlen(save))!=SUCCESS) {
       printf("\tInvalid response! Enter a valid filename suffix:\n> ");
       fgets(save, 100, stdin);
       sscanf(save, "%s", &save);
     }
     global_options->QMFILEOUTEND=(char*)malloc( (strlen(save)+1)*sizeof(char) );
     global_options->mem_allocated += (strlen(save)+1)*sizeof(char);
     strcpy(global_options->QMFILEOUTEND, save);
     if (write) fprintf(file, "QMFILEOUTEND=%s\n", save);
     
   }
   else if (input==3) {
     global_options->FUNC_TO_FIT=DIHEDRAL_LEAST_SQUARES;
     if (write) fprintf(file, "FUNC_TO_FIT=DIHEDRAL_LEAST_SQUARES\n");
      // fit phase
     printf("\nWhat do you want to fit?\n\t1) Dihedral force constants only\n\t2) Dihedral force constants and phases\n> ");
     input=get_option(1,2);
     if (input==1) {
       global_options->FIT_PHASE=NO;
       if (write) fprintf(file, "FIT_PHASE=NO\n");
     } else if (input==2) {
       global_options->FIT_PHASE=YES;
       if (write) fprintf(file, "FIT_PHASE=YES\n");
     }
   
   // Algorithm
   if (global_options->FUNC_TO_FIT!=DIHEDRAL_LEAST_SQUARES) {
    printf("\nWhich algorithm should the fitting be done with?\n\t1) Simplex\n\t2) Genetic algorithm\n\t3) Simplex, then genetic algorithm\n> ");
    input=get_option(1,3);
    if (input==1) {
      global_options->ALGORITHM=SIMPLEX;
      if (write) fprintf(file, "ALGORITHM=SIMPLEX\n");
      // Simplex options
      simplex_wizard(global_options, file);
    }
    else if (input==2) {
      global_options->ALGORITHM=GENETIC;
      if (write) fprintf(file, "ALGORITHM=GENETIC\n");
      // Genetic algorithm options
      genetic_wizard(global_options, file);
    }
    else if (input==3) {
      global_options->ALGORITHM=BOTH;
      if (write) fprintf(file, "ALGORITHM=BOTH\n");
      // Simplex and GA options
      simplex_wizard(global_options, file);
      genetic_wizard(global_options, file);
    }
   }
   
   // Output files
   printf("\nEnter a filename for saving calculated vs. original energies in to aid in gauging result quality, or hit enter\n> ");
   fgets(save, 100, stdin);
   sscanf(save, "%s", &save);
   while (save[0] != '\n' && check_for_valid_filename(save, strlen(save))!=SUCCESS) {
     printf("\tInvalid response! Enter a valid file name:\n> ");
     fgets(save, 100, stdin);
     sscanf(save, "%s", &save);
   }
   if (save[0] != '\n') {
     global_options->WRITE_ENERGY=(char*) malloc( (strlen(save)+1)*sizeof(char) );
     global_options->mem_allocated += ( strlen(save)+1 )*sizeof(char);
     strcpy(global_options->WRITE_ENERGY, save);
     if (write) fprintf(file, "WRITE_ENERGY=%s\n", save);
   }
   else printf("Not saving.\n");
   
   printf("\nEnter a filename to save a frcmod with the new parameters in, or hit enter\n> ");
   fgets(save, 100, stdin);
   sscanf(save, "%s", &save);
   while (save[0] != '\n' && check_for_valid_filename(save, strlen(save))!=SUCCESS) {
     printf("\tInvalid response! Enter a valid file name:\n> ");
     fgets(save, 100, stdin);
     sscanf(save, "%s", &save);
   }
   if (save[0] != '\n') {
     global_options->WRITE_FRCMOD=(char*) malloc( (strlen(save)+1)*sizeof(char) );
     global_options->mem_allocated += ( strlen(save)+1 )*sizeof(char);
     strcpy(global_options->WRITE_FRCMOD, save);
     if (write) fprintf(file, "WRITE_ENERGY=%s\n", save);
   }
   else printf("Not saving.\n");
  }
 
   printf("\nDo you want to save scatterplots of bonds, angles, and dihedrals in the input conformations?\n\t1) Yes\n\t2) No\n> ");
   input=get_option(1,2);
   if (input==1) {
     global_options->SCATTERPLOTS=YES;
     if (write) fprintf(file, "SCATTERPLOTS=YES\n");
   } 
   else if (input==2)
     global_options->SCATTERPLOTS=NO;
   }
     
  if (write) fclose(file);
  printf("Successfully set options!\n");
  return SUCCESS;
}

// Gets an int between min and max from stdin
int get_option(int min, int max)
{
  char raw[10];
  int input = 0;
  fgets(raw, 10, stdin);
  
  int num = sscanf(raw, "%d", &input); 
  
  while (num != 1 || input < min || input > max) {
    num = 1;
    printf("\tInvalid response! Enter a choice between %d and %d\n> ", min, max);
    fgets(raw, 10, stdin);
    num = sscanf(raw, "%d", &input);
  }
  return input;
}

double get_float()
{
  char raw[100];
  double input = 0;
  fgets(raw, 100, stdin);
  int num = sscanf(raw, "%lf", &input);
  
  while (num != 1) {
    num = 1;
    printf("\tInvalid response! Enter a valid floating point number\n> ");
    fgets(raw, 100, stdin);
    num = sscanf(raw, "%lf", &input); 
  }
  return input;
}

// Reads in options specific to the genetic algorithm
void genetic_wizard(global_options_struct *global_options, FILE *file)
{
  int input;
  printf("\n\nGENETIC ALGORITHM SETTINGS\n");
  printf("--------------------------\n\n");
  
  printf("Enter number of optimizations: (recommended 20-100)\n> ");
  global_options->NOPTIMIZATIONS=get_option(10,1000);
  if (file) fprintf(file, "OPTIMIZATIONS=%d\n", global_options->NOPTIMIZATIONS);
  
  printf("\nEnter maximum number of generations: (recommended >1000)\n> ");
  global_options->MAX_GENERATIONS=get_option(0, 10000000);
  if (file) fprintf(file, "MAX_GENERATIONS=%d\n", global_options->NOPTIMIZATIONS);
  
  printf("\nEnter number of generations that must pass without improvement for convergence to be reached: (recommend 10-100)\n> ");
  global_options->GENERATIONS_TO_CONVERGE=get_option(1, 1000);
  if (file) fprintf(file, "GENERATIONS_TO_CONV=%d\n", global_options->GENERATIONS_TO_CONVERGE);
  
  printf("\nEnter %% of original value to search around, or -1 to search entire valid solution space\n> ");
  global_options->SEARCH_SPACE=get_float(0,10000);
  while (global_options->SEARCH_SPACE != -1. && global_options->SEARCH_SPACE < 0.) {
    printf("\tInvalid response! Enter a choice greater than 0 or -1 to search entire solution space\n> ");
    global_options->SEARCH_SPACE=get_float(0,10000);
  }
  if (file) fprintf(file, "SEARCH_SPACE=%f\n", global_options->SEARCH_SPACE);
}

// Reads in options specific to the simplex algorithm
void simplex_wizard(global_options_struct *global_options, FILE *file)
{
  printf("\n\nSIMPLEX SETTINGS\n");
  printf("--------------\n\n");
  
  if (global_options->PARAMETERS_TO_FIT==K_ONLY) {
    printf("\nEnter the step size of K: (recommended 10.0)\n> ");
    global_options->K_dx = get_float(0, 1000);
    if (file) fprintf(file, "K_dx=%f\n", global_options->K_dx);
  }
  else { 
    printf("\nEnter bond force constant step size: (recommended 5.0)\n> ");
    global_options->BONDFC_dx = get_float(0,100);
    if (file) fprintf(file, "BONDFC_dx=%f\n", global_options->BONDFC_dx);
    
    printf("\nEnter bond equilibrium distance step size: (recommended 0.02)\n> ");
    global_options->BONDEQ_dx = get_float(0,1);
    if (file) fprintf(file, "BONDEQ_dx=%f\n", global_options->BONDEQ_dx);
    
    printf("\nEnter angle force constant step size: (recommended 1.0)\n> ");
    global_options->ANGLEFC_dx= get_float(0,100);
    if (file) fprintf(file, "ANGLEFC_dx=%f\n", global_options->ANGLEFC_dx);
    
    printf("\nEnter angle equilibrium value step size: (recommended 0.05)\n> ");
    global_options->ANGLEEQ_dx= get_float(0,1);
    if (file) fprintf(file, "ANGLEEQ_dx=%f\n", global_options->ANGLEEQ_dx);
    
    printf("\nEnter dihedral force constant step size: (recommended 0.2)\n> ");
    global_options->DIHEDRALBH_dx= get_float(0,5);
    if (file) fprintf(file, "DIHEDRALBH_dx=%f\n", global_options->DIHEDRALBH_dx);
    
    printf("\nEnter dihedral periodicity step size: (recommended 0.01)\n> ");
    global_options->DIHEDRALN_dx= get_float(0,1);
    if (file) fprintf(file, "DIHEDRALN_dx=%f\n", global_options->DIHEDRALN_dx);
    
    printf("\nEnter dihedral phase step size: (recommended 0.05)\n> ");
    global_options->DIHEDRALG_dx= get_float(0,1);
    if (file) fprintf(file, "DIHEDRALG_dx=%f\n", global_options->DIHEDRALG_dx);
    
    printf("\nEnter dihedral phase step size: (recommended 0.05)\n> ");
    global_options->DIHEDRALG_dx= get_float(0,1);
    if (file) fprintf(file, "DIHEDRALG_dx=%f\n", global_options->DIHEDRALG_dx);
  }
  
  printf("\nEnter a convergence limit, smaller is more stringent: (recommended 1.0E-8)\n> ");
  global_options->CONV_LIMIT = get_float(0,1);
  if (file) fprintf(file, "CONV_LIMIT=%E\n", global_options->CONV_LIMIT);
}
