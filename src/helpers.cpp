#include "atoms.h"

double kinetic_energy(Atoms &atoms){
    // Returns kinetic energy of atoms. atoms.velocities should be instantiated
    double mass = 1.0; // Until we start using different atoms

    Eigen::Array3Xd KE_arr = 0.5 * mass * atoms.velocities.square();
    return KE_arr.sum();
}


