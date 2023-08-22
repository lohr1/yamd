#include "lj_direct_summation.h"
#include "neighbors.h"

double lj_cutoff(double r, double cutoff, double epsilon, double sigma){
    // Returns value of shifted LJ potential with distance r and cutoff cutoff
    double V_r = 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
    double V_rc = 4 * epsilon * (pow(sigma / cutoff, 12) - pow(sigma / cutoff, 6));
    return V_r - V_rc;
}

double lj_cutoff_summation(Atoms &atoms, double cutoff, double epsilon, double sigma) {
    // Sums only over neighbor pairs of atoms, which are determined by cutoff.
    // LJ potential is shifted to be continuous at cutoff.
    // Modifies forces directly in atoms. Returns total potential.
    atoms.forces.setZero(); // Forces depend only on current position
    double net_pot = 0;
    int n = atoms.nb_atoms();
    double r;
    Eigen::Array3d pos_i, pos_j, r_ij, f_ij;

    // Loop through atoms, calc potential and force on each one due to LJ pot. Update force.

    // Create neighborlist
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    // Loop over each pair only once
    for (auto [i, j] : neighbor_list) {
        if (i < j) {
            // Find difference in position, r_ij, between atoms i and j
            pos_i = atoms.positions.col(i);
            pos_j = atoms.positions.col(j);
            r_ij = pos_i - pos_j;
            // Convert to Eigen Vector to calc norm, r
            Eigen::Vector3d dist_vec{r_ij};
            r = dist_vec.norm();
            // Use r to calc LJ potential (with cutoff) and add to running sum of net potential
            net_pot += lj_cutoff(r, cutoff, epsilon, sigma);
            // Force is same as in normal LJ potential:
            f_ij = (4 * epsilon *
                    ((-12 * pow(sigma, 12) / pow(r, 13)) +
                     6 * pow(sigma, 6) / pow(r, 7))) *
                   dist_vec / r;
            // Update forces for both atoms i and j:
            atoms.forces.col(j) += f_ij;
            atoms.forces.col(i) -= f_ij;
        }
    }
    // Return potential
    return net_pot;
}