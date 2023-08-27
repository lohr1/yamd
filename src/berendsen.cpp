//
// Created by robin on 05.07.23.
//

#include "atoms.h"
#include "helpers.h"

void berendsen_thermostat(Atoms &atoms, double target_temp, double timestep, double relaxation_time){
    // Directly modifies velocities in atoms.velocities
    double curr_temp = temp(kinetic_energy(atoms),atoms.nb_atoms());
    // Scaling factor lambda:
    double lambda = sqrt(1 + (target_temp /curr_temp - 1) * timestep / relaxation_time);
    // Rescale velocities
    atoms.velocities *= lambda;
}