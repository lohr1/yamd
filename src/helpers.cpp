#include "atoms.h"

// Value for boltzmann constant in eV/K
const double k_b = 8.617333e-5;

double kinetic_energy(Atoms &atoms){
    // Returns kinetic energy of atoms. Velocities should be present in atoms.
    double mass = 1.0; // Until we start using different atoms

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
