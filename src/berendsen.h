//
// Created by robin on 05.07.23.
//

#ifndef MY_MD_CODE_BERENDSEN_H
#define MY_MD_CODE_BERENDSEN_H

#include "atoms.h"

void berendsen_thermostat(Atoms &atoms, double temperature, double timestep, double relaxation_time);

#endif // MY_MD_CODE_BERENDSEN_H
