//
// Created by robin on 05.07.23.
//

#include "atoms.h"

void berendsen_thermostat(Atoms &atoms, double target_temp, double timestep, double relaxation_time){
    // Directly modifies velocities in atoms.velocities
    // Current temp T = (sum m*v^2)/(3 * N * k_B)
    double mass = 1.;
    Eigen::Array3Xd arr = mass * atoms.velocities.square();
    double curr_temp = arr.sum() / (3 * atoms.nb_atoms() * (8.617333e-5));
    // Scaling factor lambda:
    double lambda = sqrt(1 + (target_temp /curr_temp - 1) * timestep / relaxation_time);
    // Rescale velocities
    atoms.velocities *= lambda;
}
