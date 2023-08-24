#include "lj_direct_summation.h"

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma){
    // Calculates Lennard-Jones potential due to current positions in atoms.
    // Resulting forces are then stored in atoms, net potential of system is returned.

    // Forces depend only on current position, so neglect existing forces:
    atoms.forces.setZero();

    // Loop through atoms, calc potential and force on each one due to LJ pot. Update force.
    // Variables for calculation:
    double net_pot = 0;
    int n = atoms.nb_atoms();
    double r, f_mag;  // For distance and force magnitudes, resp.
    Eigen::Array3d pos_i, pos_j, r_ij, F_ij;  // Arrays to store components

    // Loop through all pairs of atoms just once:
    for(int i = 0; i < n - 1; i++){
        for(int j = i + 1; j < n; j++){
            // Get positions of atoms i and j
            pos_i = atoms.positions.col(i);
            pos_j = atoms.positions.col(j);

            // Find difference between these positions (distance vector)
            r_ij = pos_i - pos_j;

            // Find norm (magnitude) of distance vector for calculations
            Eigen::Vector3d dist_vec{r_ij};
            r = dist_vec.norm();

            // Lennard-Jones potential accumulation
            net_pot += 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));

            // Magnitude of force on atom i due to j:
            f_mag = (24*epsilon / r) * (2 * pow(sigma / r, 12) - pow(sigma/r, 6));

            // Create force vector with appropriate direction and magnitude:
            F_ij = f_mag * dist_vec /r;

            // Update forces for both atoms (Newton's third law)
            atoms.forces.col(i) += F_ij;
            atoms.forces.col(j) -= F_ij;
        }
    }
    // Return potential
    return net_pot;
}