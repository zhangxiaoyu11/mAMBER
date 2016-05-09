#ifndef ATOMIC_NUMBER_H
#define ATOMIC_NUMBER_H

#include "sff.h"

// Atomic number 0 has mass 0 and is a lone pair or extra point

#define NUM_ATOMS 119

const static int ATOMIC_MASS[NUM_ATOMS] = 
   {
     0.0, 1.0079, 4.0026, 6.941, 9.0122, 10.811, 12.0107,
     14.0067, 15.9994, 18.9984, 20.1797, 22.9898, 24.3050,
     26.9815, 28.0855, 30.9738, 32.065, 35.453, 39.948,
     39.0983, 40.078, 44.9559, 47.867, 50.9415, 51.9961,
     54.9380, 55.845, 58.9331, 58.6934, 63.546, 65.409,
     69.723, 72.64, 74.9216, 78.96, 79.904, 83.798,
     85.4678, 87.62, 88.9059, 91.224, 92.9064, 95.94,
     98., 101.07, 102.9055, 106.42, 107.8682, 112.411,
     114.818, 118.710, 121.760, 127.60, 126.9045, 131.293,
     132.9055, 137.327, 138.9055, 140.116, 140.9077, 144.242,
     145., 150.36, 151.964, 157.25, 158.9254, 162.500,
     164.9303, 167.259, 168.9342, 173.04, 174.967, 178.49,
     180.9479, 183.84, 186.207, 190.23, 192.217, 195.084,
     196.9666, 200.59, 204.3833, 207.2, 208.9804, 209.,
     210., 222., 223., 226., 227., 232.0381,
     231.0359, 238.0289, 237., 244., 243., 247.,
     247., 251., 252., 257., 258., 259.,
     262., 261., 262., 266., 264., 277.,
     268., 281., 272., 285., 284., 289.,
     288., 292., 291., 294.
   };

// This routine will return the atomic number of 
int get_atomic_number(REAL_T mass);
#endif
