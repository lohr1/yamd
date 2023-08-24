/*
 * This program was used to investigate the effect of different timesteps on
 * total energy for the provided equilibrated system "lj54.xyz"
 */

#include "lj_direct_summation.h"
#include "verlet.h"
#include "xyz.h"
#include "helpers.h"
#include <fstream>
#include <iostream>

int main(int argc, char *argv[]) {
    // Directory to store data:
    std::string out_dir = "/home/robin/School/HPC/Data/04/time_step/1e-4/";
    // lj54.xyz input file:
    std::string xyz_file_path = "/home/robin/School/yamd/milestones/04/lj54.xyz";

    auto[names, positions, velocities]{read_xyz_with_velocities(xyz_file_path)};
    Atoms atoms{positions, velocities};

    // Params
    double mass = 1.0;
    // LJ params
    double epsilon = 1.0;
    double sigma = 1.0;
    // Time params
    double t_unit = sqrt(mass * sigma * sigma / epsilon);
    double time_tot = 100 * t_unit;
    double time_step = 1e-4 * t_unit;
    int nb_steps = time_tot / time_step;

    // Variables to monitor energies:
    double KE = kinetic_energy(atoms);
    double PE = lj_direct_summation(atoms, epsilon, sigma);  // Note: this stores initial forces already in atoms
    double E_tot = KE + PE;

    // Array to store total energy throughout entire simulation
    Eigen::ArrayXd Energies(nb_steps+1);

    // Variables to regulate XYZ output
    double iter_out =
        t_unit / time_step; // Output position every iter_out iterations
    double out_thresh = iter_out; // Threshold to keep track
    std::ofstream trajectory_file(out_dir + "traj.xyz");

    for (int i = 0; i < nb_steps; ++i) {
        // Store total Energy
        Energies(i) = E_tot;

        // Verlet predictor step (changes pos and vel)
        verlet_step1(atoms, time_step);

        // Update forces with new positions (store PE here since verlet_step2 doesn't change positions)
        PE = lj_direct_summation(atoms, epsilon, sigma);

        // Verlet step 2 updates velocities, assuming the new forces are present in params:
        verlet_step2(atoms, time_step);

        // Recalculate energy
        KE = kinetic_energy(atoms);
        E_tot = PE + KE;

        // XYZ output
        if(i > out_thresh) {
            write_xyz(trajectory_file, atoms);
            out_thresh += iter_out;
        }
    }

    // Store final energy
    Energies(nb_steps) = E_tot;

    trajectory_file.close();

    // Writing Energies array to file with timestep in file name
    std::string ts_string = std::to_string(time_step);
    std::string filename = "ts_" + ts_string + "_energy.txt";
    std::ofstream outfile(out_dir + filename);
    if (outfile.is_open()) {
        outfile << Energies << std::endl;
        outfile.close();
    } else {
        std::cerr << "Error opening the file." << std::endl;
    }
    return 0;
}
