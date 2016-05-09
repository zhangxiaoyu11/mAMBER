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

/* fitting_control.c

   contains routines to control the fitting process calling
   the relevant energy minimisation routines etc.
*/

#include <stdio.h>
#include <stdlib.h>

#include "function_def.h"

int do_fit(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data)
{	
  double initial_function_value;
  double initial_r_squared_value;
  double initial_energy_struc_1;
  double final_r_squared_value;
  double final_function_value;
  double final_energy_struc_1;
  int retval;
  bounds_struct bounds_data;
  
  /* Get initial information about the structures and ensure that if any dihedral phases are being *
  * parameterised that the phases are well represented in the sample structures.                  */
  if (global_options->CHECK_BOUNDS != NO)
  {
    retval = calculate_structure_diversity(global_options, &bounds_data, parm_data, coords_data);
    if (retval != SUCCESS)
      return retval;
    retval = check_dihedrals(global_options, parm_data, &bounds_data);
    if (retval != SUCCESS)
      return retval;
  }

  /*The initial parameters and details on the fitting etc will have been listed by the options_summary
    routine. So the first thing we do here is do an initial evaluation of our function*/

  if (global_options->VERBOSITY >= MEDIUM || global_options->RANDOM_SEED != 0) // don't print for test cases
  {
    printf("      --------------------------------- INITIAL PARAMETERS ------------------------------\n");
    printf("   Initial ");
    print_parameter_summary(global_options,parm_data);
    printf("      -----------------------------------------------------------------------------------\n");
  }
  if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
  {   
      /* Calculate how many bonds, angles, and dihedrals will be fit. */
      global_options->BOND_PARAMS = calculate_no_fit_params(parm_data, BONDS);
      global_options->ANGLE_PARAMS = calculate_no_fit_params(parm_data, ANGLES);
      global_options->DIHEDRAL_PARAMS = calculate_no_fit_params(parm_data, DIHEDRALS);
      

      // Calculate initial function and R^2 value if desired
      if (global_options->VERBOSITY >= MEDIUM)
      {
        initial_function_value = eval_sum_squares_amber_std(global_options, parm_data, coords_data);
        initial_r_squared_value = calc_r_squared(global_options, parm_data, coords_data);
        initial_energy_struc_1=eval_amber_std_for_single_struct(global_options, parm_data, &coords_data[0]);
        printf("   Sum of squares for initial parameters = %15.10f kcal^2/mol^2\n",initial_function_value);
        printf("   R^2 value for initial parameters      = %10.6f\n",initial_r_squared_value);
        printf("   Calculated energy with initial parameters for structure 1 = %10.6f KCal/mol\n\n",initial_energy_struc_1);
      }
   }
  else if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
     /* Calculate how many bonds, angles, and dihedrals will be fit. */
     global_options->BOND_PARAMS = calculate_no_fit_params(parm_data, BONDS);
     global_options->ANGLE_PARAMS = calculate_no_fit_params(parm_data, ANGLES);
     global_options->DIHEDRAL_PARAMS = calculate_no_fit_params(parm_data, DIHEDRALS);
      
    // Mark the atoms that are to be fit
    parm_data->fit_atom = (int *)malloc(global_options->NATOMS*sizeof(int));
    if (parm_data->fit_atom == NULL)
      malloc_failure_int("do_fit", "parm_data->fit_atom", global_options->NATOMS);
    mark_relevant_atoms(global_options, parm_data);
     
    // Calculate initial function and R^2 value if desired
    if (global_options->VERBOSITY >= MEDIUM)
    {
      initial_function_value = eval_sum_amber_forces(global_options, parm_data, coords_data);
      initial_r_squared_value = calc_r_squared(global_options, parm_data, coords_data);
      initial_energy_struc_1=eval_amber_std_for_single_struct(global_options, parm_data, &coords_data[0]);
      
      force_struct *eval = (force_struct*)malloc(global_options->NATOMS*sizeof(force_struct));
      if (eval==NULL) return ALLOC_FAIL;
      eval_amber_forces_single_struct(global_options, parm_data, coords_data, eval, 0);
    
      printf("   Sum of force difference magnitudes for initial parameters = %15.10f kcal/mol-A\n", initial_function_value);
      printf("   R^2 value for initial parameters                          = %10.6f\n", initial_r_squared_value);
      printf("   Calculated energy with initial parameters for structure 1 = %10.6f KCal/mol\n",initial_energy_struc_1);
      printf("   Calculated forces with initial parameters atom 1, structure 1 = %3.4f %3.4f %3.4f kcal/mol-A\n\n", eval[0].x, eval[0].y, eval[0].z);
      free(eval);
    }
  }
  else if (global_options->FUNC_TO_FIT==DIHEDRAL_LEAST_SQUARES) {
    /* Calculate how many bonds, angles, and dihedrals will be fit. */
    printf("!  Warning- setting all dihedral terms to fit pk\n");
    if (global_options->FIT_PHASE) printf("!  Warning- setting all dihedral terms to fit phase\n");
    else printf("!  Warning- fitting no dihedral phases\n");
    int i;
    for (i=0; i<parm_data->unique_dihedrals_found; ++i) {
      parm_data->dihedral_data[i].DO_FIT_KP=YES;
      if(global_options->FIT_PHASE)
        parm_data->dihedral_data[i].DO_FIT_PHASE=YES;
      else
        parm_data->dihedral_data[i].DO_FIT_PHASE=NO;
    }
    
    global_options->BOND_PARAMS = calculate_no_fit_params(parm_data, BONDS);
    global_options->ANGLE_PARAMS = calculate_no_fit_params(parm_data, ANGLES);
    global_options->DIHEDRAL_PARAMS = calculate_no_fit_params(parm_data, DIHEDRALS);
    printf("   Using dihedral minimization algorithm described by Chad Hopkins and Adrian Roitberg.\n");
  }
  else
  {
    printf("   ERROR IN DO_FIT - FUNCTION %d IS NOT YET IMPLEMENTED\n",global_options->FUNC_TO_FIT);
    return NOT_IMPLEMENTED;
  }
   
  // If K is being fit along with other parameters, give an error
  if (global_options->K_FIT==YES && global_options->PARAMETERS_TO_FIT!=K_ONLY)
  {
    printf("   ERROR: K is being fit along with other parameters. This will produce inaccurate\n");
    printf("         results and is disallowed. Specify a value for K or set PARAMETERS_TO_FIT=K_ONLY\n");
    return INVALID_DATA;
  }
  
  if (global_options->FUNC_TO_FIT==DIHEDRAL_LEAST_SQUARES) { // fitting function and algorithm combined
    dihedral_least_squares(global_options, parm_data, coords_data, global_options->FIT_PHASE);
  } 
  else if (global_options->ALGORITHM==SIMPLEX) {
     retval=minimise_function_simplex(global_options, parm_data, coords_data);
  } 
  else if (global_options->ALGORITHM==GENETIC)
  {
    if (global_options->VERBOSITY>=HIGH)
      print_parameter_summary(global_options,parm_data);
    retval = minimise_function_genetic(global_options, parm_data, coords_data);
  }
  else if (global_options->ALGORITHM==BOTH)
  {
    retval = minimise_function_genetic(global_options, parm_data, coords_data);
    retval = minimise_function_simplex(global_options, parm_data, coords_data);
  }
  else if (global_options->ALGORITHM==NONE)
  {
    // Filler for printing information or something
  }
  else
   {
     printf("   ERROR IN do_fit() - UNKNOWN FITTING FUNCTION: %d\n",global_options->ALGORITHM);
     return UNKNOWN_OPT;
   }
   if (retval!=SUCCESS )
   {
      /*We can trap one error here - that is if we quit due to exceeding maximum iterations*/
      /*Other errors are fatal*/
      if (retval==EXCEEDEDMAXITERATIONS || retval==MINSTATIC)
        printf("!  Warning - final parameters do NOT represent a converged fit.\n");
      else
        return(retval);
   }

    /* Check that the results are reasonable based on the given data in the input structures.       *
    * Dihedrals are already checked at the beginning, but do bonds and angles here to make sure    *
    * the result is in a well-represented area in the input data.                                  */
    if (global_options->CHECK_BOUNDS)
    {
      retval = check_angles(global_options, parm_data, &bounds_data);
      retval = check_bonds(global_options, parm_data, &bounds_data);
      clean_up_bounds(&bounds_data);
      if (retval != SUCCESS) return FAILURE;
    }
   
    // Print out a final summary
    if (global_options->VERBOSITY>=MEDIUM)
    {
      printf("   ---------------------------------- FINAL PARAMETERS -------------------------------\n");
      printf("   Fitted ");
      print_parameter_summary(global_options,parm_data);
      printf("   -----------------------------------------------------------------------------------\n");
      
      if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD) {
        final_function_value = eval_sum_squares_amber_std(global_options, parm_data, coords_data);
        final_r_squared_value=calc_r_squared(global_options, parm_data, coords_data);
        final_energy_struc_1=eval_amber_std_for_single_struct(global_options, parm_data, &coords_data[0]);
        
        printf("   Function value with fitted parameters  =  %12.4f, R^2 = %12.4f\n",final_function_value,final_r_squared_value);
        printf("   Calculated energy with fitted parameters for structure 1 = %11.4f KCal/mol\n",final_energy_struc_1);
        printf("\n");
      }
      else if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
        final_function_value = eval_sum_amber_forces(global_options, parm_data, coords_data);
        final_r_squared_value=calc_r_squared(global_options, parm_data, coords_data);
        final_energy_struc_1=eval_amber_std_for_single_struct(global_options, parm_data, &coords_data[0]);
        
        force_struct *eval = (force_struct*)malloc(global_options->NATOMS*sizeof(force_struct));
        if (eval==NULL) return ALLOC_FAIL;
        eval_amber_forces_single_struct(global_options, parm_data, coords_data, eval, 0);
        
        printf("   Function value with fitted parameters  =  %12.4f, R^2 = %12.4f\n",final_function_value,final_r_squared_value);
        printf("   Calculated energy with fitted parameters for structure 1 = %11.4f KCal/mol\n",final_energy_struc_1);
        printf("   Calculated forces with fitted parameters atom 1, structure 1 = %3.4f %3.4f %3.4f kcal/mol-A\n\n", eval[0].x, eval[0].y, eval[0].z);
        printf("\n");
        free(eval);
      
      // TEST TODO print out difference files
      printf("   Experimental: Writing theta.txt and mag.txt with difference for atom 0\n");
      print_forces(global_options, parm_data, coords_data, 0);
      }
    }
    
    // if desired, save the frcmod file
    if (global_options->WRITE_FRCMOD)
      write_frcmod(global_options, parm_data);
    
    // if desired, write the list of final energies for comparison
    if (global_options->WRITE_ENERGY)
      write_energy(global_options, parm_data, coords_data, -1);
    
    // if desired, save a new toplogy and parameter file
    if (global_options->WRITE_PRMTOP)
      write_prmtop(global_options, parm_data);

   return SUCCESS;
}


