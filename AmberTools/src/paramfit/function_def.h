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

/*function_def.h*/

#include "prmtop_params.h"
#include "forces.h"

#include <stdio.h>

/*Contains constant declarations, function definitions etc.*/

/*Return codes*/
/*
      Note, greater than ABORT = not a failure, an exit code
      Less than ABORT is a real failure if not zero
      0 = SUCCESS - success
      -1 = FAILURE - Generic failure
      -2 = NOT_IMPLEMENTED not implemented
      -101 = UNKNOWN_OPT failure unknown option
      -201 = TOO_MANY_OPT failure too many options
      -301 = ALLOC_FAIL allocation failure
      -401 = FILE_OPEN_FAIL File open failure - Read
      -1000 = ABORT
      -1001 = CMD_HELP_REQ User requested command line help
      -1002 = HIST_REQ User requested program history
      -1003 = HELP_REQ User requested a help topic, print it and return 6
      -2001 = INVALID_FORMAT
      -2002 = INVALID_DATA
      -2003 = INVALID_LINE
      -3001 = EXCEEDEDMAXITERATIONS - EXCEEDED MAXIMUM ITERATIONS
      -3002 = MINSTATIC Function no longer varying
      -3003 = DATA_OVERFLOW
      -3004 = UNKNOWN_ELEMENT
      
Note lower than -1000 are abort codes, not error codes
*/
#define _GNU_SOURCE

#define SUCCESS  -0
#define FAILURE  -1
#define NOT_IMPLEMENTED -2
#define UNKNOWN_OPT    -101
#define TOO_MANY_OPT   -201
#define ALLOC_FAIL     -301
#define FILE_OPEN_FAIL -401
#define FILE_READ_FAIL -501
#define ABORT        -1000
#define CMD_HELP_REQ -1001
#define HIST_REQ     -1002
#define HELP_REQ     -1003
#define INVALID_FORMAT -2001
#define INVALID_DATA   -2002
#define INVALID_LINE   -2003
#define EXCEEDEDMAXITERATIONS -3001
#define MINSTATIC -3002
#define DATA_OVERFLOW -3003
#define UNKNOWN_ELEMENT -3004

#define OFF 0
#define ON  1
#define DEBUG 2
#define WARN 3

#define BONDS 1
#define ANGLES 2
#define DIHEDRALS 3

/*MAXIMUM VALUES FOR ITERATIONS ETC*/
#define NSIMPLEX_INNER_PER_DIM 25 /*This is multiplied by the number of dimensions in order to work out how many inner loops to run*/
#define NSIMPLEX_OUTER_MAX 100000
#define RAND_RATIO    0.2  /*Ratio by which to multiply the random number received from rand() when
adjusting simplex vertices by gamma*/                                                                

/*MINIMISING FUNCTIONS AVAILABLE*/
typedef enum
{
  SUM_SQUARES_AMBER_STANDARD,
  AMBER_FORCES,
  DIHEDRAL_LEAST_SQUARES
} function_t;

/*MINIMISATION ALGORITHMS AVAILABLE*/
typedef enum
{
  SIMPLEX,
  GENETIC,
  BOTH,
  NONE
} algorithm_t;

/*QM FILE FORMATS AVAILABLE*/
typedef enum
{
  GAUSSIAN,
  ADF,
  GAMESS
} qm_format_t ;

typedef enum
{
  YES  = 1,
  TRUE = 1,
  NO   = 0,
  FALSE = 0
} bool_t;

typedef enum
{
  READ,
  WRITE
} readwrite_t;

typedef enum
{
  DEFAULT,                               // fit a default set of parameters
  LOAD,                                  // read in parameters from a previously created file
  SAVE,                                  // save parameters to an output file
  K_ONLY                                 // only fit K
} parameter_mode_t;

typedef enum
{
  CREATE_INPUT,                          // write input files for quantum package
  FIT,                                   // do the actual fitting
  SET_PARAMS                             // set which parameters should be fit
} runtype_t;

typedef enum
{
  HARTREE,
  KCALMOL,
  KJMOL
} energy_t;

typedef enum
{
  HARTREE_BOHR,
  KCALMOL_ANGSTROM
} force_t;

/*VERBOSITY LEVELS*/
typedef enum
{
  LOW,
  MEDIUM,
  HIGH
} verbosity_t;

/*DEFINITIONS OF STRUCTURES*/

/*global_options_struct - This structure contains all of the global options used by the program
                         for things like number of QMTerms to expect, settings filenames etc */
typedef struct _global_options_struct
{
   int mem_allocated;                 // Updated with number of bytes allocated for pointer inside this structure
   
   /* Command line switches */
   verbosity_t VERBOSITY;             // How verbose to be - low, medium, or high
   char *job_control_filename;        // Contains the path and filename of the control file 
   char *prmtop_filename;             // Contains the path and filename of the topology file 
   char *mdcrd_filename;              // Contains the path and filename of the coordinate file 
   char *energy_filename;             // Contains the path and filename of the energy file which contains the energies to be fitted to

   /* General options */
   char *PARAMETER_FILE_NAME;         // Filename with parameters to be fit 
   char *WRITE_ENERGY;                // Filename, if any, to save final qm and md energies of structures to 
   char *WRITE_FRCMOD;                // Filename, if any, to save ffrcmod to 
   char *WRITE_PRMTOP;                // Filename, if any, to save a new prmtop to 
   int RANDOM_SEED;                   // for duplicating runs if necessary, for debugging usually 
   runtype_t RUNTYPE;
   qm_format_t QMFILEFORMAT;
   energy_t QM_ENERGY_UNITS;          // Unit of energy to expect from QM data file 
   force_t QM_FORCE_UNITS;            // Unit of force from QM data file
   parameter_mode_t PARAMETERS_TO_FIT;
   algorithm_t ALGORITHM;             // Fitting routine to be used 
   function_t FUNC_TO_FIT;            // The function to be used to fit to the energy surface 
   bool_t K_FIT;                      // whether or not to fit the K parameter

   /* Things to be fit */
   double K;                          // intrinsic difference between quantum and classical
   int BOND_PARAMS;                   // stores number of parameters of each type to be fit
   int ANGLE_PARAMS;
   int DIHEDRAL_PARAMS;
   int NDIHEDRALS;                    // number of terms to give fitted dihedrals, at minimum
   
   /* Genetic algorithm options */
   int NOPTIMIZATIONS;                // number of parameter sets per generation
   int MAX_GENERATIONS;               // maximum number of generations to run
   int GENERATIONS_TO_CONVERGE;       // number generations in a row without a change to end algorithm
   double SEARCH_SPACE;               // distance away from initial parameter set to search
   
   /* Simplex algorithm options */
   double CONV_LIMIT;                 // convergence limit
   double BONDFC_dx;                  // simplex step sizes
   double BONDEQ_dx;
   double ANGLEFC_dx;
   double ANGLEEQ_dx;
   double DIHEDRALBH_dx;              // ! note dihedral BH in prmtop file is actually Vn/2
   double DIHEDRALN_dx;
   double DIHEDRALG_dx;
   double K_dx;
   
   /* Bounds checking options */
   bool_t CHECK_BOUNDS;               // whether or not to check bounds
   double ANGLE_LIMIT;                // converged angle theta must be this close to values in input structures
   double BOND_LIMIT;                 // converged bond length must be this close to values in input structures
   int DIHEDRAL_SPAN;                 // each dihedral must be spanned by this many input structures
   bool_t SCATTERPLOTS;               // whether or not to write scatter plots with input and output equilibrium parameters
   
   /* Molecule information options */
   int NDIMENSIONS;                   // number of dimensions of fit, will be calculated
   int NSTRUCTURES;                   // number of input structures
   int NATOMS;                        // number of atoms in the system
   double SCNB;                       // 1-4 scaling factors for use with standard amber force field equation
   double SCEE; 

   /* Quantum input file creation options */
   char *QMFILEOUTSTART;              // filename for QM input files- format is startNNNend where NNN is structure number
   char *QMFILEOUTEND;
   char *QMHEADER;                    // stores the location of a file to go at the beginning of qm input
   int QM_SYSTEM_CHARGE;              // integral charge of the system
   int QM_SYSTEM_MULTIPLICITY;        // integral multiplicity of the system
   
   /* Dihedral least squares fit options */
   bool_t FIT_PHASE;
   
} global_options_struct;

typedef struct _coords_struct
{
   int mem_allocated;                 // updated for amount of memory allocated per structure
   double *x_coord;                   // x, y, and z coordinates of each atom in the structure
   double *y_coord;
   double *z_coord;
   double energy;                     // qm energy for this coordinate set in Kcal/mol
   force_struct *force;               // will contain forces on each atom if fitting forces
} coords_struct;

typedef struct _bounds_struct
{
  // these hold the values at which bond length, angles, and dihedrals are defined
  double **bond_lengths;
  double **angle_thetas;
  double **dihedral_thetas;
  int mem_allocated;
} bounds_struct;

/*Externs for global variables*/

/*Function defs*/
void print_program_info(void);
void print_program_history(void);
int set_default_options(global_options_struct *global_options);
void process_retval(int err_code, verbosity_t VERBOSITY);
void malloc_failure_char(char *routine, char *var_name, int chars_requested);
void malloc_failure_int(char *routine, char *var_name, int ints_requested);
void malloc_failure_short_int(char *routine, char *var_name, int short_ints_requested);
void malloc_failure_double(char *routine, char *var_name, int doubles_requested);
void file_open_failure(char *routine, char *var_name);
int process_command_line(int argc, char **argv, global_options_struct *global_options);
void command_line_help(char *cmd_line_options[]);
double **alloc_2D_double(int nrows, int ncolumns);
void global_unlock(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data);
int process_job_control_setting(char *setting_line, int length, int *number_settings, global_options_struct *global_options);
int read_job_control_file(global_options_struct *global_options);
void print_job_control_summary(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data);
void print_close_line_box(int no_spaces);
void print_open_line_box(int *i);
int check_for_valid_filename(const char *data_string, const int length);
int read_prmtop(global_options_struct *global_options, parm_struct *parm_data);
int s_getline(char *line, int max, FILE *fp);
int find_flag( FILE *fptr, char *label );
int name_copy( FILE *fptr, char *stringp);
int read_mdcrd(global_options_struct *global_options, coords_struct *coords_data);
int create_input(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data);
int write_input_gaussian(global_options_struct *global_options, parm_struct *parm_data, coords_struct *current_struct, int num, FILE *fptr);
int find_atomic_number_from_parm(parm_struct *parm_data, int atom);
void print_atomic_number_as_symbol(FILE *fptr, int atomic_number);
int read_qm_energy(global_options_struct *global_options, coords_struct *coords_data);
double calc_r_squared(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data);
int unObfuscateAtom( int at );
int ObfuscateAtom( int at ); 
double eval_amber_std_for_single_struct(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data);
void print_parameter_summary(global_options_struct *global_options, parm_struct *parm_data);
int process_prmtop(global_options_struct *global_options, parm_struct *parm_data);
double calc_bond_length(double bond1x, double bond1y, double bond1z, double bond2x, double bond2y, double bond2z);
double calc_angle_radians(double atom1x, double atom1y, double atom1z, double atom2x, double atom2y, double atom2z,
                          double atom3x, double atom3y, double atom3z);
double calc_dihedral_radians(double atom1x, double atom1y, double atom1z, double atom2x, double atom2y, double atom2z,
                          double atom3x, double atom3y, double atom3z, double atom4x, double atom4y, double atom4z);
void calc_fit_dimensions(global_options_struct *global_options, parm_struct *parm_data);
int modify_params_scratch_data(global_options_struct *global_options, parm_struct *parm_data, double *parameters, readwrite_t MODE);
int calculate_no_fit_params(parm_struct *parm_data, short int MODE);
void double_2D_array_free(double **array);

int write_input_parameters(global_options_struct *global_options, parm_struct *parm_data);
int read_input_parameters(global_options_struct *global_options, parm_struct *parm_data);
int read_parameter_file(global_options_struct *global_options, parm_struct *parm_data);
int calculate_structure_diversity(global_options_struct *global_options, bounds_struct *bounds_data, parm_struct *parm_data, coords_struct *coords_data);
void clean_up_bounds(bounds_struct *bounds_data);
int check_dihedrals(global_options_struct *global_options, parm_struct *parm_data, bounds_struct *bounds_data);
int check_bonds(global_options_struct *global_options, parm_struct *parm_data, bounds_struct *bounds_data);
int check_angles(global_options_struct *global_options, parm_struct *parm_data, bounds_struct *bounds_data);
int write_frcmod(global_options_struct *global_options, parm_struct *parm_data);
void handle_sigint(int param);
void print_backtrace(int signal);
coords_struct* alloc_coords(global_options_struct *global_options);
void free_coords(global_options_struct *global_options, coords_struct *coords_data);
int dihedral_types_equal(dihedral_data_struct *first, dihedral_data_struct *second);
void print_dihedral(dihedral_data_struct *hi);
int write_input_adf(global_options_struct *global_options, parm_struct *parm_data, coords_struct *current_struct, int num, FILE *fptr);
int write_input_gamess(global_options_struct *global_options, parm_struct *parm_data, coords_struct *current_struct, int num, FILE *fptr);
int write_energy(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data, int generation);
int compare_energy(const void *a, const void *b); // this is for qsort.
int not_enough_dihedrals(parm_struct *parm_data, int n);
int write_prmtop(global_options_struct *global_options, parm_struct *parm_data);
int write_mdcrd(global_options_struct *global_options, coords_struct *coords_data);
int do_fit(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data);
double eval_sum_squares_amber_std(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data);
int minimise_function_simplex(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data);  
void check_range(global_options_struct *global_options, parm_struct *parm_data);

/* Forces algorithms */
int read_gaussian_forces(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data);
int eval_amber_forces_single_struct(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data, force_struct *forces, int structure);
double eval_sum_amber_forces(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data); 
void print_forces(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data, int atom);

/* Genetic algorithm functions */
int minimise_function_genetic(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data);
int do_mutation(global_options_struct *global_options, parm_struct *parm_data, double *row, int col, bool_t do_mutate);
double **alloc_data_matrix(int rows, int cols);
void free_data_matrix(double **dm, int rows);

/* Wizard functions */
int job_control_wizard(global_options_struct *global_options);
int get_option(int min, int max);

/* New dihedral fittings */
int dihedral_least_squares(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data, bool_t fit_phase);
int conduct_dihedral_least_squares(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data, 
                                       bool_t fit_phase);
