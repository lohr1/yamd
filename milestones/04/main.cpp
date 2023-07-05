//
// Created by robin on 24.05.23.
//

#include "lj_direct_summation.h"
#include "verlet.h"
#include "xyz.h"
#include <fstream>
#include <iostream>

int main(int argc, char *argv[]) {
    auto[names, positions, velocities]{read_xyz_with_velocities("/home/robin/School/yamd/milestones/04/lj54.xyz")};
    Atoms atoms{positions, velocities};

    // Integrate
    double mass = 1.0;
    double epsilon = 1.0;
    double sigma = 1.0;

    double expr = sqrt(mass * sigma * sigma / epsilon);
    double time_tot = 100 * expr;
    double time_step = 0.001 * expr;
    int nb_steps = time_tot / time_step;


    Eigen::Array3Xd KE_arr = 0.5 * mass * atoms.velocities.square();
    double KE = KE_arr.sum();
    double PE, E_tot1, E_tot2;

    // Calc forces and store potential energy
    PE = lj_direct_summation(atoms, epsilon, sigma);
    E_tot1 = KE + PE;

    // To store Energies (row 0) and diff's (row 1)
    Eigen::Array2Xd Energies{2, nb_steps};

    // Variables for XYZ output
    double out_thresh = expr;
    std::ofstream traj("traj.xyz");

    for (int i = 0; i < nb_steps; ++i) {
        Energies(0, i) = E_tot1;

        // Verlet predictor step (changes pos and vel)
        verlet_step1(atoms, time_step);

        // Update forces with new positions (also save PE here since verlet_step2 doesn't change positions)
        PE = lj_direct_summation(atoms, epsilon, sigma);

        // Verlet step 2 updates velocities assumes the new forces are present in params:
        verlet_step2(atoms, time_step);

        // Recalculate energy
        KE_arr = 0.5 * mass * atoms.velocities.square();
        KE = KE_arr.sum();
        E_tot2 = PE + KE;
        double E_diff = E_tot2 - E_tot1;

        if(i % 50 == 0) {
            std::cout << E_tot1 << ' ' << E_diff << std::endl;
        }
        Energies(1, i) = E_diff;
        E_tot1 = E_tot2;

        // XYZ output
        if(i > out_thresh) {
            write_xyz(traj, atoms);
            out_thresh += expr;
        }
    }

    traj.close();

    // Writing Energies to file
    std::ofstream outfile("Energies.txt");
    if (outfile.is_open()) {
        outfile << Energies << std::endl;
        double tot_diff = Energies(0, 0) - Energies(0, Energies.cols()-1);
        outfile << "Total Energy difference: " << tot_diff << std::endl;
        std::cout << "Total Energy difference: " << tot_diff << std::endl;
        std::cout << "Made it." << std::endl;
        outfile.close();
    } else {
        std::cerr << "Error opening the file." << std::endl;
    }
    return 0;
}
