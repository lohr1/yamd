//
// Created by robin on 05.07.23.
//

#ifndef MY_MD_CODE_HELPERS_H
#define MY_MD_CODE_HELPERS_H

#include "atoms.h"

// Returns Kinetic Energy of atoms
double kinetic_energy(Atoms &atoms);

// Returns Temp
double temp(std::double_t KE, int nb_atoms);

#endif // MY_MD_CODE_HELPERS_H
