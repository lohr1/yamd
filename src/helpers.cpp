#include "atoms.h"

// Value for boltzmann constant in eV/K
const double k_b = 8.617333e-5;

double kinetic_energy_local(Atoms& atoms, int nb_local){
    /*
     * Returns kinetic energy of local atoms.
     */
    Eigen::ArrayXd squared_velocities = atoms.velocities.colwise().squaredNorm();
    Eigen::ArrayXd ke_per_particle = 0.5 * atoms.masses * squared_velocities;
    ke_per_particle.conservativeResize(nb_local);
    return ke_per_particle.sum();
}

double kinetic_energy(Atoms& atoms){
    /*
     * Returns total kinetic energy of atoms.
     */
    Eigen::ArrayXd squared_velocities = atoms.velocities.colwise().squaredNorm();
    Eigen::ArrayXd ke_per_particle = 0.5 * atoms.masses * squared_velocities;

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
