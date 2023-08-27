#include "atoms.h"

// Value for boltzmann constant in eV/K
const double k_b = 8.617333e-5;

double kinetic_energy(Atoms &atoms){
    // Returns kinetic energy of atoms if mass of all particles is 1. Velocities should be present in atoms.
    double mass = 1.0;

    Eigen::ArrayXd squared_velocities = atoms.velocities.colwise().squaredNorm();
    Eigen::ArrayXd ke_per_particle = 0.5 * mass * squared_velocities;

    return ke_per_particle.sum();
}

double kinetic_energy(Atoms& atoms, double mass){
    /*
     * Returns kinetic energy of atoms if all particles are of identical mass
     * using EAM units (mass in g/mol, time in 10.18 fs).
     */
    Eigen::ArrayXd squared_velocities = atoms.velocities.colwise().squaredNorm();
    Eigen::ArrayXd ke_per_particle = 0.5 * mass * squared_velocities;

    return ke_per_particle.sum();
}

double temp(double KE, int nb_atoms){
    // Returns temperature given KE and nb_atoms
    return 2 * KE / (3 * nb_atoms * k_b);
}

double average_velocity(double temp){
    // Returns magnitude of average velocity in 3 dimensions given temperature
    return sqrt(3 * k_b * temp);
}
