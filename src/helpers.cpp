#include "atoms.h"

double kinetic_energy(Atoms &atoms){
    // Returns kinetic energy of atoms. Velocities should be present in atoms.
    double mass = 1.0; // Until we start using different atoms

    Eigen::ArrayXd squared_velocities = atoms.velocities.colwise().squaredNorm();
    Eigen::ArrayXd ke_per_particle = 0.5 * mass * squared_velocities;

    return ke_per_particle.sum();
}

double temp(double KE, int nb_atoms){
    return 2 * KE / (3 * nb_atoms * 8.617333e-5);
}


