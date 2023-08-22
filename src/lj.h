//
// Created by robin on 24.05.23.
//

#ifndef MY_MD_CODE_LJ_DIRECT_SUMMATION_H
#define MY_MD_CODE_LJ_DIRECT_SUMMATION_H

#include "atoms.h"
double lj_cutoff_summation(Atoms &atoms, double cutoff = 5.0, double epsilon = 1.0, double sigma = 1.0);

#endif // MY_MD_CODE_LJ_DIRECT_SUMMATION_H
