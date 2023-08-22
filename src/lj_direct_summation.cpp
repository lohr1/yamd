#include "lj_direct_summation.h"

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma){
    // Directly modifies forces in atoms. Returns potential energy.
    atoms.forces.setZero(); // Forces depend only on current position
    // Loop through atoms, calc potential and force on each one due to LJ pot. Update force.
    double net_pot = 0;

    int n = atoms.nb_atoms();
    double r;
    Eigen::Array3d pos_i, pos_j, r_ij, f_ij;
    for(int i = 0; i < n - 1; i++){
        pos_i = atoms.positions.col(i);
        for(int j = i + 1; j < n; j++){
            pos_j = atoms.positions.col(j);
            r_ij = pos_i - pos_j;
            Eigen::Vector3d dist_vec{r_ij};
            r = dist_vec.norm();
            //r = r_ij.matrix().norm();
            net_pot += 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
            // Force on atom i due to j:
            //f_ij = (24*epsilon / r_ij) * (-2 * pow(sigma / r_ij, 12) + pow(sigma/r_ij, 6));
            f_ij = (4 * epsilon * ((-12 * pow(sigma, 12)/pow(r, 13)) + 6 * pow(sigma, 6) / pow(r, 7))) *
                   dist_vec /r;

            atoms.forces.col(j) += f_ij;
            atoms.forces.col(i) -= f_ij;
        }
    }
    // Return potential
    return net_pot;
}