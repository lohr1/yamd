/*
 * This file contains miscellaneous helper functions. Some compute values while
 * others help with filesystem operations.
 */

#ifndef MY_MD_CODE_HELPERS_H
#define MY_MD_CODE_HELPERS_H

#include "atoms.h"


// Returns Kinetic Energy of atoms
double kinetic_energy(Atoms& atoms);
double kinetic_energy_local(Atoms& atoms, int nb_local);

// Returns Temp
double temp(std::double_t KE, int nb_atoms);

#endif // MY_MD_CODE_HELPERS_H
