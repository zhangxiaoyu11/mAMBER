#ifndef ChargeFitDataStructures
#define ChargeFitDataStructures

#include "TopologyDS.h"

struct FittingPoint {
  int ix;          // Grid index in X
  int iy;          // Grid index in Y
  int iz;          // Grid index in Z
  int flagged;     // Flag to indicate that this point is eligible (0) or
                   //   ineligible (1) as fitting data
  double minr;     // Minimum distance from the molecule (includes distance to
                   //   extra points that may have been added at run-time)
};
typedef struct FittingPoint fitpt;

struct RestraintData {
  int natm;        // The number of atoms in this restraint
  double target;   // The target value of the charge for this restraint
  double mult;     // The strength with which to retain this restraint
  int* atoms;      // The list of atoms involved in this restraint
  double* wt;      // The weights by which each atom's value contributes to
                   //   the restraint ("atom A must have twice the charge of
                   //   atom B")
  char* maskstr;   // The mask string which leads to the list of atoms
};
typedef struct RestraintData nail;

struct FittingSet {
  int ngrd;          // The number of grids available for the fitting
  int nfitpt;        // The number of fitting points to seek out
  int exclusive;     // Implies that any particular data point can only be used
                     //   in the fitting or testing sets, respectively
  int nqeq;          // Number of charge equalization restraints
  int nqmin;         // Number of charge minimization restraints
  int nqsum;         // Number of charge group sum restraints
  int tpcount;       // The number of unique topologies to deal with in this
                     //   fit
  int q2fit;         // The number of unique charges to fit
  int DispAllDP;     // Flag to toggle display of ALL dipole moments, for all
                     //   fitted conformations (default 0, do not display)
  int verbose;       // Flag to print output relating to fitting run progress
                     //   (default 1, print output)
  int* fitpthist;    // Histogram of fitting point distances (bin width
                     //   fhistbin)
  int* tpidx;        // The topology identification numbers 
  double fitprob;    // Probability that a certain data point (defined as a 
                     //   a point in the electrostatic potential grid) will
                     //   wind up in the fitting set
  double testprob;   // Probability that a data point will wind up in the test
                     //   set
  double flimit;     // Proximity limit of fitting points from the same grid
  double peps;       // Probe epsilon value (default eps for TIP4P-Ew oxygen)
  double psig;       // Probe sigma value (default sigma for TIP4P-Ew oxygen)
  double prbarm;     // Probe arm (by default, the O-H distance in TIP waters)
  double stericlim;  // Energetic cutoff at which the probe is no longer
                     //   tolerated near the surface of the molecule
  double qminwt;     // The stiffness by which individual charges which are to
                     //   be restrained are held to zero
  double fhistbin;   // Bin width for fitting point distance-to-molecule
                     //   histogram
  double Rc;         // The cutoff beyond which the acceptance probability for
                     //   new fitting points drops below 1 (default 3.0)
  double Rmax;       // The cutoff beyond which the acceptance probability for
                     //   new fitting points is zero (default 6.0)
  double* totalq;    // Total charge on each system, enforced during fit
  double* wt;        // Weight assinged to data points extracted from each grid
                     //   in the overall fit
  cmat gname;        // Names of all grid files
  cmat auxgname;     // Names of all auxiliary grid files
  cmat tpname;       // Names of topology files for individual grids
  cmat eprule;       // Names of extra point rules files for individual grids
  nail* qeq;         // Charge equalization restraints
  nail* qmin;        // Charge minimzation restraints
  nail* qsum;        // Charge group sum restraints
  prmtop* TPbank;    // Topology bank for systems describing all potential
                     //   grids
  long long int MaxMem;    // The maximum memory that the machine is assumed to
                           // have available.  Default 1GB.
  char epext[MAXNAME];     // Extension of new extra points file to be written
                           //   which will modify topologies in future
                           //   simulations
  char confext[MAXNAME];   // For printing conformations of each system to
                           //   PDB files
  char histfile[MAXNAME];  // For printing the extra points distribution
                           //   (distance from solute to extra point)
};
typedef struct FittingSet fset;

#endif
