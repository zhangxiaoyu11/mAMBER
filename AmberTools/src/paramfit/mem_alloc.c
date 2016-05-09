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

/*mem_alloc.c*/
/*Memory allocation routines*/

#include <stdlib.h>
#include <stdio.h>

#include "function_def.h"

double **alloc_2D_double(int nrows, int ncolumns)
{
  /*Allocates a 2d_double_array consisting of a series of pointers pointing to each
    row that are then allocated to be ncolumns long each.*/

  /*Uses calloc - slower but all locations will be zeroed*/
  
  /*Tries to keep contents contiguous - thus reallocation is difficult!*/

  /*Returns the pointer **array. Returns NULL on error*/
  int i;
  
  double **array = (double **)calloc(nrows,sizeof(double *));
  if (array==NULL)
    return NULL;
  array[0] = (double *)calloc(nrows * ncolumns,sizeof(double));
  if (array[0]==NULL)
     return NULL;
  
  for (i = 1; i < nrows; ++i)
    array[i] = array[0] + i * ncolumns;

  return array;
}

void double_2D_array_free(double **array)
{
  /*Frees the memory previously allocated by alloc_2D_double*/
  free(array[0]);
  free(array);
}

void global_unlock(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data)
{
  free(global_options->energy_filename);    
  free(global_options->mdcrd_filename);     
  free(global_options->prmtop_filename);       
  global_options->energy_filename = NULL;
  global_options->mdcrd_filename = NULL;
  global_options->prmtop_filename = NULL;
  
  if (global_options->job_control_filename) {
    free(global_options->job_control_filename); 
    global_options->job_control_filename = NULL;
  }
  if (global_options->WRITE_FRCMOD)
  {
    free(global_options->WRITE_FRCMOD);
    global_options->WRITE_FRCMOD=NULL;
  }
  if (global_options->QMHEADER)
  {
    free(global_options->QMHEADER);
    global_options->QMHEADER=NULL;
  }
  if (global_options->WRITE_ENERGY)
  {
    free(global_options->WRITE_ENERGY);
    global_options->WRITE_ENERGY=NULL;
  }
  if (global_options->PARAMETER_FILE_NAME)
  {
    free(global_options->PARAMETER_FILE_NAME);
    global_options->PARAMETER_FILE_NAME = NULL;
  }
  
  if (global_options->QMFILEOUTSTART) {
    free(global_options->QMFILEOUTSTART);
    global_options->QMFILEOUTSTART = NULL;
  }
  if (global_options->QMFILEOUTEND) {
    free(global_options->QMFILEOUTEND);
    global_options->QMFILEOUTEND = NULL;
  }

if (global_options->RUNTYPE==FIT)
{
    free_coords(global_options, coords_data);
    free(parm_data->bond_data);
    free(parm_data->angle_data);
    free(parm_data->dihedral_data);
    
    parm_data->bond_data = NULL;
    parm_data->angle_data = NULL;
    parm_data->dihedral_data = NULL;
    if (global_options->FUNC_TO_FIT==AMBER_FORCES)
      free(parm_data->fit_atom);
  }

  if ( parm_data->NHB > 0 )
  {
      free(parm_data->bg);  
      free(parm_data->ag); 
      
      parm_data->bg = NULL;
      parm_data->ag = NULL;
  }

  if ( parm_data->NEXT > 0 )
  {
    free(parm_data->natex); 
    parm_data->natex = NULL;
  }
  if (parm_data->MPHIA > 0)
  {
    free(parm_data->pdihedral); 
    parm_data->pdihedral = NULL;
  }
  if (parm_data->NPHIH > 0)
  {
    free(parm_data->pdihedralH); 
    parm_data->pdihedralH = NULL;
  }

  if (parm_data->MTHETS > 0)
  {
    free(parm_data->pangle);
    parm_data->pangle = NULL;
  }
  if (parm_data->NTHETH > 0 )
  {
    free(parm_data->pangleH);
    parm_data->pangleH = NULL;
  }
  if (parm_data->NBONA > 0 )
  {
    free(parm_data->pbond);
    parm_data->pbond = NULL;
  }
  if (parm_data->NBONH > 0)
  {
    free(parm_data->pbondH);
    parm_data->pbondH = NULL;
  }

  if (parm_data->NTYPES > 0)
  {
    free(parm_data->cn2);
    free(parm_data->cn1);
    
    parm_data->cn2 = NULL;
    parm_data->cn1 = NULL;
  }

  if (parm_data->NATYP > 0)
  {
    free(parm_data->solty);   
    parm_data->solty = NULL;
  }
  if (parm_data->MPTRA > 0)
  {
    free(parm_data->phase);  
    free(parm_data->pn);    
    free(parm_data->pk);
  
    parm_data->phase = NULL;
    parm_data->pn = NULL;
    parm_data->pk = NULL;
  }
  if (parm_data->MUMANG)
  {
    free(parm_data->teq);
    free(parm_data->tk);
  
    parm_data->teq = NULL;
    parm_data->tk = NULL;
  }
  if (parm_data->MUMBND)
  {
    free(parm_data->req);
    free(parm_data->rk);
  
    parm_data->req = NULL;
    parm_data->rk = NULL;
  }
  if (parm_data->NTOTRS > 0)
  {
    free(parm_data->residue);
    parm_data->residue = NULL;
  }  
  if (parm_data->NTYPES > 0)
  {
      free(parm_data->nno);
      parm_data->nno = NULL;
  }
  if ( parm_data->NTOTAT > 0)
  {
      free(parm_data->atom);   
      parm_data->atom = NULL;
  }
  
  free(parm_data->title);
  parm_data->title = NULL;
    
}

coords_struct* alloc_coords(global_options_struct *global_options)
{ 
    coords_struct *coords_data = (struct _coords_struct*) malloc(global_options->NSTRUCTURES*sizeof(coords_struct));
    if (coords_data == NULL)
    {
      malloc_failure_char("alloc_coords", "coords_data", global_options->NSTRUCTURES*sizeof(coords_struct));
      return NULL;
    }
    int i,j;
    
    for (i=0; i<global_options->NSTRUCTURES; ++i)
    {
      coords_data[i].mem_allocated = 0;
      coords_data[i].x_coord = (double *) malloc(global_options->NATOMS*sizeof(double));
      if (coords_data[i].x_coord == NULL)
      {
        malloc_failure_char("alloc_coords", "coords_data->x_coord",(global_options->NATOMS*sizeof(double)));
        return NULL;
      }
      coords_data[i].mem_allocated+=(global_options->NATOMS*sizeof(double));
      coords_data[i].y_coord = (double *) malloc(global_options->NATOMS*sizeof(double));
      if (coords_data[i].y_coord == NULL)
      {
        malloc_failure_char("alloc_coords", "coords_data->y_coord",(global_options->NATOMS*sizeof(double)));
        return NULL;
      }
      coords_data[i].mem_allocated+=(global_options->NATOMS*sizeof(double));
      coords_data[i].z_coord = (double *) malloc(global_options->NATOMS*sizeof(double));
      if (coords_data[i].z_coord == NULL)
      {
	      malloc_failure_char("alloc_coords", "coords_data->z_coord",(global_options->NATOMS*sizeof(double)));
        return NULL;
      }
      coords_data[i].mem_allocated+=(global_options->NATOMS*sizeof(double));
      
      // If fitting forces, allocate space for the force on each atom
      if (global_options->FUNC_TO_FIT==AMBER_FORCES) {
        coords_data[i].force = (force_struct *) malloc(global_options->NATOMS*sizeof(force_struct));
        if (coords_data[i].force == NULL) {
          malloc_failure_char("alloc_coords", "coords_data->force", (global_options->NATOMS*sizeof(force_struct)));
          return NULL;
        }
      }
      else
        coords_data[i].force = NULL;
    }
    return coords_data;
}

void free_coords(global_options_struct *global_options, coords_struct *coords_data)
{
  int i;
  for (i=0; i<global_options->NSTRUCTURES; ++i)
  {
    free(coords_data[i].x_coord);
    free(coords_data[i].y_coord);
    free(coords_data[i].z_coord);
    if (global_options->FUNC_TO_FIT==AMBER_FORCES)
      free(coords_data[i].force);
  }
  free(coords_data);
  coords_data = NULL;
}
