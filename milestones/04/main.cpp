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
    double PE, E_tot;

    // Calc forces and store potential energy
    PE = lj_direct_summation(atoms, epsilon, sigma);
    E_tot = KE + PE;

    // To store Total energy as function of time
    Eigen::ArrayXd Energies(nb_steps+1);

    // Variables for XYZ output
    double iter_out = expr / time_step; // Output position every iter_out iterations
    double out_thresh = iter_out; // Threshold to keep track
    std::ofstream traj("traj.xyz");

    for (int i = 0; i < nb_steps; ++i) {
        // Store total Energy
        Energies(i) = E_tot;

        // Verlet predictor step (changes pos and vel)
        verlet_step1(atoms, time_step);

        // Update forces with new positions (also save PE here since verlet_step2 doesn't change positions)
        PE = lj_direct_summation(atoms, epsilon, sigma);

        // Verlet step 2 updates velocities assumes the new forces are present in params:
        verlet_step2(atoms, time_step);

        // Recalculate energy
        KE_arr = 0.5 * mass * atoms.velocities.square();
        KE = KE_arr.sum();
        E_tot = PE + KE;

        // XYZ output
        if(i > out_thresh) {
            write_xyz(traj, atoms);
            out_thresh += iter_out;
        }
    }

    // Store final energy
    Energies(nb_steps) = E_tot;

    traj.close();

    // Writing Energies to file with time_step in filename
    std::string ts_string = std::to_string(time_step);
    std::string filename = "ts_" + ts_string + "_energy.txt";
    std::ofstream outfile(filename);
    if (outfile.is_open()) {
        outfile << Energies << std::endl;
        outfile.close();
    } else {
        std::cerr << "Error opening the file." << std::endl;
    }
    return 0;
}
