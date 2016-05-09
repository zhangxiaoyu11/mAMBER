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

/*genetic_algorithm.c*/

/*Contains the genetic algorithm fitting routine:*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "function_def.h"
#include "constants.h"

int minimise_function_genetic(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data)
{
  
  /* variables to use go here */
  double **data_matrix;           // will hold all of the possible optimizations
  double *function_results;       // will hold the amber function results for each optimization
  double previous_best;           // used for calculating the convergence ratio
  double initial_result;          // for comparison with the converged results
  double mean=0.0;                // useful statistics
  double mutation_rate = 0.35;    // percentage of values randomly changed each generation
  double parent_percent = 0.25;   // this constitutes the top percent of solutions
  double weighting;               // used to change values by a certain amount
  
  double conv_f;                  // keeps track of difference in best fitness from this generation to next
  double init_f;                  // keeps track of previous generation's fitness

  int generation;                 // how many generations have passed
  int converged;                  // how generations in a row have passed with a low first derivative
  int function_calls=0;           // keeps track of how many times the amber function has been called
  int col, row;
  int retval;
  int k_offset = (global_options->K_FIT==YES) ? 1 : 0;
  int i;                          // counting variable for loops
  
  /* print out a header with some information */
  if (global_options->VERBOSITY>=MEDIUM)
  {
    printf("   ---------------------------GENETIC ALGORITHM MINIMISATION ---------------------------\n");
      if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
        printf("   Minimising function SUM_SQUARES_AMBER_STANDARD, using the GENETIC ALGORITHM\n");
      else if (global_options->FUNC_TO_FIT==AMBER_FORCES)
        printf("   Minimising function AMBER_FORCES using the GENETIC ALGORITHM\n");
      else
        printf("   Minimising function UNKNOWN, using the GENETIC ALGORITHM\n");      
      if (global_options->VERBOSITY>=MEDIUM)
        printf("   ------------------------------------ CONVERGENCE ------------------------------------\n");
      else
        printf("   -------------------------------------------------------------------------------------\n");
      fflush(stdout); /*Flush the printf buffer*/
  }
 
 /* Allocate the 1D matrix to hold the function results for each optimization */
 if (global_options->VERBOSITY>=HIGH)
   printf("Allocating %d bytes for *function_results.\n",(int)((global_options->NOPTIMIZATIONS)*sizeof(double)));
 function_results=(double *)calloc(global_options->NOPTIMIZATIONS, sizeof(double));
 if (function_results==NULL)
 {
   malloc_failure_double("minimise_function_genetic", "function_results", (global_options->NOPTIMIZATIONS));
   return ALLOC_FAIL;
 }
 
  /* Allocate space for the optimizations */
  /* The data matrix is 2 dimensional, with num_optimizations rows and NDIMENSIONS cols */
  data_matrix=alloc_data_matrix(global_options->NOPTIMIZATIONS, global_options->NDIMENSIONS);
  if (data_matrix==NULL)
  {
    printf("   ERROR - MALLOC FAILED FOR data_matrix\n");
    printf("   COMMAND WAS data_matrix=alloc_2D_double(global_options->NOPTIMIZATIONS,global_options->NDIMENSIONS)\n");
    return ALLOC_FAIL;
  }

  if ( (int)(global_options->NOPTIMIZATIONS*parent_percent) <= 2)
  {
    printf("   ERROR: NOPTIMIZATIONS=%d is too small to conduct the algorithm.\n", global_options->NOPTIMIZATIONS);
    free_data_matrix(data_matrix, global_options->NOPTIMIZATIONS);
    free(function_results);
    return INVALID_DATA;
  }
  
// end memory allocation
  
  /* Initially fill the data matrix with the initial parameters. The first *
    * potential solution will be the starting parameters and will be        *
    * changed only at the end of the initial generation creation.           */
  retval = modify_params_scratch_data(global_options, parm_data, data_matrix[0], READ);
  if (retval != global_options->NDIMENSIONS)
  {
    printf("   ERROR IN minimise_function_genetic()\n");
    printf("             modify_params_scratch_data() DOES NOT MATCH\n");
    printf("             NDIMENSIONS OF %d. ABORTING RUN.\n",global_options->NDIMENSIONS);
    fflush(stdout); /*Flush the printf buffer*/
    free_data_matrix(data_matrix, global_options->NOPTIMIZATIONS);
    free(function_results);
    /*This is most probably a bug so return an unknown error*/
    return(FAILURE);
  }
  for (row=1; row<global_options->NOPTIMIZATIONS; ++row)
    for (col=0; col<global_options->NDIMENSIONS; ++col)
      data_matrix[row][col] = data_matrix[0][col];
 
  /* Evaluate the function for the initial parameters */
  initial_result = eval_sum_squares_amber_std(global_options, parm_data, coords_data);
  int params_found = 0;
  int count = 0;
  double delta = 0;

  if (global_options->VERBOSITY>=HIGH)
  {
    printf("ORIGINAL:\n");
    for (col=0; col<global_options->NDIMENSIONS; ++col)
      printf("%8.3f", data_matrix[0][col]);
    printf("\tFUNCTION: %f\n---\n", initial_result);
  }
  
  /* Create the initial generation randomly. Some values, like K, can be any number to start.      *
  * Others will be random within an allowed domain, such as angles and dihedral phases, using the *
  * initial parameters as a starting point.                                                       */
  for (row=0; row < global_options->NOPTIMIZATIONS; ++row)
  {
    params_found = 0;
    for(col=0; col<global_options->NDIMENSIONS; ++col)
    {
      // search space value determines how far to look or to look everywhere
      if (global_options->SEARCH_SPACE <= 0.0) // search everything
      {
        delta = 0.0;
        data_matrix[row][col] = 0.0;
        if (global_options->K_FIT==YES)
          global_options->K = (rand() & 1) ? rand() : -1.0*rand();
      }
      else
      {
        // deviation from initial parameter is random, can be positive or negative
        delta = global_options->SEARCH_SPACE * ( (double)rand()/(double)RAND_MAX );
        if (rand() & 1) delta*= -1.0;
      }

      // K can be any value, positive or negative, is not based on initial parameters unless
      // an initial guess is provided
      if (col==0 && global_options->K_FIT==YES && global_options->SEARCH_SPACE > 0.0)
      {
          data_matrix[row][col] += data_matrix[row][col] * delta;
      }
      else if (col-k_offset<global_options->BOND_PARAMS)
      {
        /*Current column represents a bonding parameter. Loop through the bonds until we find it*/
        params_found=0;
        for (count=0;count<parm_data->unique_bonds_found;++count)
        {
          // Base Kr off the initial value or make it between 100 and 1000
          if (parm_data->bond_data[count].DO_FIT_KR==YES)
          {
            /*At this point we check to see if our params_found to date match the column we are on*/
            if (params_found==col-k_offset)
            {
              if (data_matrix[row][col] == 0.0)
          data_matrix[row][col] = 900*((double)rand()/(double)RAND_MAX) + 100.0;
              else
          data_matrix[row][col] += data_matrix[row][col] * delta;
              /*no point continuing for this column*/
              break;
            }
            ++params_found;
          }
          // Req is based on the initial value or between 0 and 3
          if (parm_data->bond_data[count].DO_FIT_REQ==YES)
          {
            if (params_found==col-k_offset)
            {
              if (data_matrix[row][col] == 0.0)
                data_matrix[row][col] = 3.0*((double)rand()/(double)RAND_MAX);
              else
                data_matrix[row][col] += data_matrix[row][col] * delta;
              break;
            }
            ++params_found;
          }
        }
      }
      else if (col-k_offset<global_options->ANGLE_PARAMS+global_options->BOND_PARAMS)
      {
        /*Current column represents an angle parameter. Loop through the angles until we find it*/
        params_found=0;
        for (count=0;count<parm_data->unique_angles_found;++count)
        {
          // Kt based off initial parameters or between 0 and 170
          if (parm_data->angle_data[count].DO_FIT_KT==YES)
          {
            if (params_found+global_options->BOND_PARAMS==col-k_offset)
            {
              if (data_matrix[row][col] == 0.0)
                data_matrix[row][col] = 170.0*((double)rand()/(double)RAND_MAX);
              else
                data_matrix[row][col] += data_matrix[row][col] * delta;
              break;
            }
            ++params_found;
          }
          // Theta needs to stay between zero and PI
          if (parm_data->angle_data[count].DO_FIT_THEQ==YES)
          {
            if (params_found+global_options->BOND_PARAMS==col-k_offset)
            {
              data_matrix[row][col] += data_matrix[row][col]*delta;
              if (data_matrix[row][col] <= 0.0 || data_matrix[row][col] > PI)
                data_matrix[row][col] = ((double)rand()/(double)RAND_MAX) * PI;
              break;
            }
            ++params_found;
          }
        }
      }
      else if (col-k_offset<global_options->DIHEDRAL_PARAMS+global_options->ANGLE_PARAMS+global_options->BOND_PARAMS)
      {
        /*Current column represents a dihedral parameter. Loop through the dihedrals until we find it*/
        params_found=0;
        for (count=0;count<parm_data->unique_dihedrals_found;++count)
        {
          // Kp based on initial value or between -30 and 30
          if (parm_data->dihedral_data[count].DO_FIT_KP==YES)
          {
            if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset)
            {
              if (data_matrix[row][col] <= 0.0)
                data_matrix[row][col] = -30.0 + 60.0*((double)rand()/(double)RAND_MAX);
              else
                data_matrix[row][col] += data_matrix[row][col]*delta;
              break;
            }
            ++params_found;
          }
          // Np based off initial value, is an integer, or between 1 and 5
          if (parm_data->dihedral_data[count].DO_FIT_NP==YES)
          {
            if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset)
            {
              if (fabs(data_matrix[row][col])>5.0 || data_matrix[row][col]==0.0)
              {
                data_matrix[row][col] = (double)(rand() % 6) + 1.0;
                if (rand()&1) data_matrix[row][col] *= -1.0;
              }
              else
                data_matrix[row][col] += (double)( (data_matrix[row][col]*delta) );
              break;
            }
            ++params_found;
          }
          // Phase between 0 and PI, based off initial value or random
          if (parm_data->dihedral_data[count].DO_FIT_PHASE==YES)
          {
            if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset)
            {
              data_matrix[row][col] += data_matrix[row][col]*delta;
              if (data_matrix[row][col] <= 0.0 || data_matrix[row][col] > PI)
                data_matrix[row][col] =  PI*((double)rand()/(double)RAND_MAX);
              break;
            }
            ++params_found;
          }
        }
      }
    }
  }
  
/* Evaluate the function for each of the initial optimisations and sort them by this fitness so    *
 * that the top row contains the best optimisation so far.                                         */
for(row=0; row<global_options->NOPTIMIZATIONS; ++row)
{
  // put the current parameters into parm data for the function evaluation
  retval = modify_params_scratch_data(global_options, parm_data, data_matrix[row], WRITE);
  if (retval != global_options->NDIMENSIONS) // check result is as expected
  {
    printf("   ERROR IN minimise_function_genetic() - RETURN VALUE OF %d FROM\n",retval);
    printf("             modify_params_scratch_data() DOES NOT MATCH\n");
    printf("             NDIMENSIONS OF %d. ABORTING RUN.\n",global_options->NDIMENSIONS);
    fflush(stdout);
    free_data_matrix(data_matrix, global_options->NOPTIMIZATIONS);
    free(function_results);
    return(FAILURE);
  }
  
  // evaluate function for the row- evaluation done in parallel within the function
  if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
  {
    function_results[row] = eval_sum_squares_amber_std(global_options, parm_data, coords_data);
    ++function_calls;
  }
  else if (global_options->FUNC_TO_FIT==AMBER_FORCES)
  {
    function_results[row] = eval_sum_amber_forces(global_options, parm_data, coords_data);
    ++function_calls;
  }
  else // function hasn't been implemented
  {
    printf("   ERROR: FUNCTION %d IS NOT YET IMPLEMENTED\n",global_options->FUNC_TO_FIT);
    free_data_matrix(data_matrix, global_options->NOPTIMIZATIONS);
    free(function_results);
    return(NOT_IMPLEMENTED);
  }
} // end fitness function

  // sort the results by fitness (insertion sort)
  // data[0] will point to the row with smallest function evaluation
  for (row=1; row<global_options->NOPTIMIZATIONS; ++row)
  {
    double value = function_results[row];
    double *temp = data_matrix[row];
    
    col = row-1;
    while(col >=0 && function_results[col] > value)
    {
      data_matrix[col+1] = data_matrix[col];
      function_results[col+1] = function_results[col];
      col--;
    }
    data_matrix[col+1] = temp;
    function_results[col+1] = value;
  }
  
  init_f = function_results[0];
  
  generation = 0;
  converged  = 0;

/* This is the main loop, which continues through the process of recombination, function          *
 * evaluation, and bounds checking until there is neglible improvement for                        *
 * global_options->GENERATIONS_TO_CONVERGE or the maximum number of generations is hit.           */
  while (generation < global_options->MAX_GENERATIONS && (converged < global_options->GENERATIONS_TO_CONVERGE || generation < 100))
  {
    // hold on to stats from the previous generation
    previous_best = function_results[0];
      
  /* Do recombination and mutation row-by row in this large loop */
    int parent1, parent2;
    for (count=0; count<global_options->NOPTIMIZATIONS-(int)(global_options->NOPTIMIZATIONS*parent_percent); ++count) // make NOPTIMIZATIONS children each generation
    {
      params_found=0;
      row = (rand() % (int)((1.0-parent_percent)*global_options->NOPTIMIZATIONS + 0.5)) + (int)(global_options->NOPTIMIZATIONS*parent_percent);
      
      // Do recombination, check for clones and mutate on a column-by-column basis- treat variables as independent of one another.
      for (col=0; col<global_options->NDIMENSIONS; ++col)
      {
        parent2 = rand() % global_options->NOPTIMIZATIONS*parent_percent;
        /* Recombination */
        do {
          parent1 = rand() % (int)(global_options->NOPTIMIZATIONS*parent_percent);
        } while (parent1==parent2);

        // Child value is somewhere in the range between the two parents
        weighting = ( (double)rand() / (double)RAND_MAX );
        if (data_matrix[parent1][col] > data_matrix[parent2][col])
          data_matrix[row][col] = data_matrix[parent2][col] + weighting*(data_matrix[parent1][col]-data_matrix[parent2][col]);
        else
          data_matrix[row][col] = data_matrix[parent1][col] + weighting*(data_matrix[parent2][col]-data_matrix[parent1][col]);

        /* Mutation */
        if (rand() % (int)100*mutation_rate == 0)
          do_mutation(global_options, parm_data, data_matrix[row], col, YES);
        
        // Check for clones of this value, if it exists, mutate this value away
          for (i=0; i<global_options->NOPTIMIZATIONS; ++i)
          {
            if (i==row) continue;    // don't check it against itself
            if (fabs(data_matrix[row][col] - data_matrix[i][col]) < 1.0e-10)
              do_mutation(global_options, parm_data, data_matrix[row], col, YES);
          }
      } // end col loop
	
    if (modify_params_scratch_data(global_options, parm_data, data_matrix[row], WRITE) != global_options->NDIMENSIONS)
      {
        printf("   ERROR IN minimise_function_genetic() - RETURN VALUE OF FROM\n");
        printf("             modify_params_scratch_data() DOES NOT MATCH\n");
        printf("             NDIMENSIONS OF %d. ABORTING RUN.\n",global_options->NDIMENSIONS);
        fflush(stdout);
        free_data_matrix(data_matrix, global_options->NOPTIMIZATIONS);
        free(function_results);
        return(FAILURE);
      } 
      
    // Get the fitness of the new child optimization
    if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
    {
      function_results[row] = eval_sum_squares_amber_std(global_options, parm_data, coords_data);
      ++function_calls;
    }
    else if (global_options->FUNC_TO_FIT==AMBER_FORCES)
    {
      function_results[row] = eval_sum_amber_forces(global_options, parm_data, coords_data);
      ++function_calls;
    }
    else
    {
      printf("   ERROR: FUNCTION %d IS NOT YET IMPLEMENTED\n",global_options->FUNC_TO_FIT);
      free_data_matrix(data_matrix, global_options->NOPTIMIZATIONS);
      free(function_results);
      return(NOT_IMPLEMENTED);
    }
    
    // sort the array by fitness
    for (row=1; row<global_options->NOPTIMIZATIONS; ++row)
    {
      double value = function_results[row];
      double *temp = data_matrix[row];
      col = row-1;
      while(col >=0 && function_results[col] > value)
      {
        data_matrix[col+1] = data_matrix[col];
        function_results[col+1] = function_results[col];
        col--;
      }
      data_matrix[col+1] = temp;
      function_results[col+1] = value;
    }
    /* Force the whole bottom 20% to be completely mutated- this can pop it out of local minima. *
     * The fitness of these won't be evaluated until after the next child is created, but that's *
     * okay since they're too low in the hierarchy to get to be parents most of the time.        */
    for (row=global_options->NOPTIMIZATIONS-(0.20*global_options->NOPTIMIZATIONS); row<global_options->NOPTIMIZATIONS; ++row)
      for (col=0; col<global_options->NDIMENSIONS; ++col)
        do_mutation(global_options, parm_data, data_matrix[row], col, YES);
    
  } // end recombination counter
  
  
  /* Check for convergence- after a while the algorithm will begin to stagnate, improving the best  *
   * optimization less and less often. When there is no improvement in the best solution for many   *
   * generations in a row (when the first derivative is zero), the algorithm is considered to have  *
   * converged. The number of generations that must happen with this condition is held in           *
   * global_options->GENERATIONS_TO_CONVERGE.                                                       */
  
    /* Print the top five -- for debugging */
    // for(row=0; row<global_options->NOPTIMIZATIONS; row+=5)
    // {
    //       for (col=0; col<global_options->NDIMENSIONS; ++col)
    //         printf("%7.2f", data_matrix[row][col]);
    //       printf("\t\t%.4E\tR%d\n", function_results[row], row);
    // }
    
    if ( fabs(previous_best - function_results[0]) < 1.0e-6 )
      ++converged;
    else
      converged = 0;
      
    // calculate some useful statistics
    mean = 0.0;
    for (row=0; row<global_options->NOPTIMIZATIONS; ++row)
      mean += function_results[row];
    mean /= (double)global_options->NOPTIMIZATIONS;
    
    conv_f = previous_best-function_results[0];

    if (generation < 100)
    printf("   Gen %4i:  Best= %10.5f  \tMean= %10.5g Elapsed=%10d/%d\n",
           generation, function_results[0], mean, generation, 100);
    else
      printf("   Gen %4i:  Best= %10.5f  \tMean= %10.5g Conv=%10d/%d\n",
      generation, function_results[0], mean, converged, global_options->GENERATIONS_TO_CONVERGE);
      

    fflush(stdout);
    
    ++generation;
  } // end main generation loop
  
  // Once convergence has been reached, copy the answer into the parm_data struct
  printf("| Took    %d generations to converge.\n", generation);
  retval = modify_params_scratch_data(global_options, parm_data, data_matrix[0], WRITE);
  if (retval != global_options->NDIMENSIONS) // check result is as expected
  {
    printf("    ERROR IN minimise_function_genetic() - RETURN VALUE OF %d FROM\n",retval);
    printf("             modify_params_scratch_data() DOES NOT MATCH\n");
    printf("             NDIMENSIONS OF %d. ABORTING RUN.\n",global_options->NDIMENSIONS);
    fflush(stdout);
    free_data_matrix(data_matrix, global_options->NOPTIMIZATIONS);
    free(function_results);
    return(FAILURE);
  }
  if (global_options->VERBOSITY>=MEDIUM)
    printf("   Called the amber function %d times.\n", function_calls);

  // Do some cleanup
  free_data_matrix(data_matrix, global_options->NOPTIMIZATIONS);
  free(function_results);
  function_results = NULL;
    
  // Warn if it's exited because of having passed the max number of generations
  if (generation >= global_options->MAX_GENERATIONS)
  {
    printf("\n");
    printf("   WARNING: Generation count of %d exceeded maximum of %d.\n", generation, global_options->MAX_GENERATIONS);
    printf("            Either 1) Fit fewer dimensions or             \n");
    printf("                   2) Increase MAX_GENERATIONS.         \n\n");
    fflush(stdout);
    return EXCEEDEDMAXITERATIONS;
  }

  return SUCCESS;
 }

// does mutation or a validity check for the given element in the given row
int do_mutation(global_options_struct *global_options, parm_struct *parm_data, double *row, int col, bool_t do_mutate)
{
  double mutation_amount = 0.10;      // how much mutation can change each value             
  double delta = 0.0;                 // how much value will change by (Check bounds = 0)
  int k_offset = (global_options->K_FIT==YES) ? 1 : 0;
  int params_found = 0;
  int count = 0;

  if (do_mutate==YES)
  {
    delta = (double)rand() / ((double)RAND_MAX+1) * mutation_amount; // a value > 0, < rec_deviation
    if (rand() & 1) delta*=-1.0;
    // If the value is zero, change it a bit so mutation will work and it will stay nonzero
    if (row[col] == 0.0) row[col] = (rand() & 1) ? 1.0 : -1.0;
  }

  /* This part either mutates or does a bounds check. If it is set to mutate, delta is nonzero
  * and so the data point will be changed. Otherwise, delta is zero so nothing is added to the
  * data and a bounds check happens.
  */
  if (col==0 && global_options->K_FIT==YES)
  {
    // Don't mutate K too much since it is probably a big number and will just kill your fitness
    // since K doesn't have a local minimum and is fit pretty much instantly. Instead let it
    // change by +- 10.0 for small adjustments when other parameters change.
    row[col] += ( (double)rand() / ((double)RAND_MAX+1) )*10.0;
  }
  else if (col-k_offset<global_options->BOND_PARAMS)
  {
    params_found=0;
    for (count=0; count<parm_data->unique_bonds_found; ++count)
    {
      if (parm_data->bond_data[count].DO_FIT_KR==YES)
      {
        if (params_found==col-k_offset)
        {
          row[col] += row[col] * delta;
	  // Bond force constant between 100 and 1000
	  if (row[col] < 100.0 || row[col] > 1000.0)
	    row[col] = 900.0*( (double)rand()/(double)RAND_MAX ) + 100.0;
          break;
        }
        ++params_found;
      }
      if (parm_data->bond_data[count].DO_FIT_REQ==YES)
      {
        if (params_found==col-k_offset)
        {
	  // Bond length between 0 and 3 Angstroms
	  if (row[col] <= 0.0 || row[col] > 3.0)
	    row[col] = 3.0*( (double)rand()/(double)RAND_MAX );
          row[col] += row[col]*delta;
          break;
        }
        ++params_found;
      }
    }
  }
  else if (col-k_offset<global_options->ANGLE_PARAMS+global_options->BOND_PARAMS)
  {
    params_found=0;
    for (count=0;count<parm_data->unique_angles_found;++count)
    {
      if (parm_data->angle_data[count].DO_FIT_KT==YES)
      {
        if (params_found+global_options->BOND_PARAMS==col-k_offset)
        {
	  // Angle force constant between 0 and 170
	  if (row[col] <= 0.0 || row[col] > 170.0)
	    row[col] = 170.0*( (double)rand() / (double)RAND_MAX );
          row[col] += row[col]*delta;
          break;
        }
        ++params_found;
      }
      // Theta needs to stay between zero and PI
      if (parm_data->angle_data[count].DO_FIT_THEQ==YES)
      {
        if (params_found+global_options->BOND_PARAMS==col-k_offset)
        {
          row[col] += row[col]*delta;
          
          if (row[col] < 0.0 || row[col] > PI)
          {
            row[col] = ((double)rand()/(double)RAND_MAX) * PI;
            if (do_mutate==NO) printf("Fixed angle %i theta!\n", col);
          }
          break;
        }
        ++params_found;
      }
    }
  }
  else if (col-k_offset<global_options->DIHEDRAL_PARAMS+global_options->ANGLE_PARAMS+global_options->BOND_PARAMS)
  {
    params_found=0;
    for (count=0;count<parm_data->unique_dihedrals_found;++count)
    {
      // Kr based on initial value
      if (parm_data->dihedral_data[count].DO_FIT_KP==YES)
      {
        if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset)
        {
	  row[col] += row[col]*delta;
	  // Dihedral force constant between 0 and 30
	  if (row[col] < -30.0 || row[col] > 30.0)
	    row[col] = -30.0 + 60.0*( (double)rand() / (double)RAND_MAX );
          break;
        }
        ++params_found;
      }
      // Np positive or negative, integer in range 1-5
      if (parm_data->dihedral_data[count].DO_FIT_NP==YES)
      {
        if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset)
        {
          row[col] += row[col]*delta;
          
          if (row[col] < -5.0 || row[col] > 5.0)
          {
            row[col] = (double)((rand() % 6) + 1);
            if (rand()&1) row[col]*=-1.0;
          }
	  
	  // ensure that it is an integer
	    row[col] = (double)((int)(row[col]+0.5));
          break;
        }
        ++params_found;
      }
      // Phase between 0 and PI, based off initial value
      if (parm_data->dihedral_data[count].DO_FIT_PHASE==YES)
      {
        if (params_found+global_options->BOND_PARAMS+global_options->ANGLE_PARAMS==col-k_offset)
        {
          row[col] += row[col]*delta;
          
          if (row[col] < 0.0 || row[col] > PI)
          {
            if (do_mutate==NO) printf("Fixed dihedral %i phase!\n", col);
            row[col] = ((double)rand()/(double)RAND_MAX) * PI;
          }
          break;
        }
        ++params_found;
      }
    }
  }
  return SUCCESS;
}


/* Allocates a more traditional, non-contiguous 2D array for the genetic algorithm
 * We don't use alloc_2d_double because that gets a contiguous block of memory and
 * we specifically need non-contiguous since we will be swapping the rows around
 * in sorting each generation, and free-ing the contiguous matrix later is a pain
 * since it is now unclear which was the original first row.
 */
 double **alloc_data_matrix(int rows, int cols) {
   double **dm = (double **)calloc(rows, sizeof(double *));
   if (dm==NULL) return NULL;
   
   int i;
   for (i=0; i<rows; ++i) {
     dm[i] = (double *)calloc(cols, sizeof(double));
     if (dm[i]==NULL) return NULL;
   }
   return dm;
 }
 
 void free_data_matrix(double **dm, int rows) {
   int i;
   for (i=0; i<rows; ++i)
     free(dm[i]);
   free(dm);
 }
     













