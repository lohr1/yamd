#include "lj_direct_summation.h"
#include "../../cmake-build-debug/_deps/eigen3-src/Eigen/Core"
#include <iostream>
//
// Created by robin on 24.05.23.
//
double lj_direct_summation(Atoms &atoms, double epsilon, double sigma){
    // Directly modifies forces in atoms. Returns potential energy.
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
//
//#include "lj_direct_summation.h"
//
//double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
//    double const sigmaPower12 = pow(sigma, 12);
//    double const sigmaPower6 = pow(sigma, 6);
//
//    atoms.forces.setZero();
//
//    double potential_energy = 0;
//
//    // loop over all pairs of atom
//    for (int i = 0; i < atoms.nb_atoms() - 1; ++i) {
//        for (int k = i + 1; k < atoms.nb_atoms(); ++k) {
//            // compute pair term (derivative of potential energy with respect to atom k)
//            const Eigen::Vector3d distance_vector{atoms.positions.col(i) - atoms.positions.col(k)};
//            const double distance = distance_vector.norm();
//            const Eigen::Array3d term{4 * epsilon * ((-12 * sigmaPower12) / pow(distance, 13) + (6 * sigmaPower6) / pow(distance, 7) ) * (distance_vector / distance)};
//
//            // add it to both forces arrays
//            atoms.forces.col(k) += term;
//            atoms.forces.col(i) -= term;
//
//            // compute potential energy (leave out constant epsilon * 4 part)
//            potential_energy += (sigmaPower12/pow(distance, 12) - sigmaPower6/pow(distance, 6));
//        }
//    }
//
//    return potential_energy * epsilon * 4;
//}