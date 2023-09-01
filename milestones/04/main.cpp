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

void simulate_ts(std::string ts_str, double ts){
    std::string out_dir = "/home/robin/School/HPC/Data/04/ReportData/";

    // Filepath to store energy data:
    std::string e_filepath = out_dir + ts_str + "_energy.csv";
    // Filepath for trajectory:
    std::string trajectory_filepath = out_dir + ts_str + "_trajectory.xyz";

    // lj54.xyz input file:
    std::string xyz_file_path = "/home/robin/School/yamd/milestones/04/lj54.xyz";

    auto[names, positions, velocities]{read_xyz_with_velocities(xyz_file_path)};
    Atoms atoms{positions, velocities};

    // Params
    double mass = 1.0;
    atoms.masses.setConstant(mass);
    // LJ params
    double epsilon = 1.0;
    double sigma = 1.0;
    // Time params
    double t_unit = sqrt(mass * sigma * sigma / epsilon);
    double time_tot = 100 * t_unit;
    double time_step = ts * t_unit;
    int nb_steps = time_tot / time_step;

    // Variables to monitor energies:
    double KE = kinetic_energy(atoms);
    double PE = lj_direct_summation(atoms, epsilon, sigma);  // Note: this stores initial forces already in atoms
    double E_tot = KE + PE;

    // Array to store total energy throughout entire simulation (before + after)
    Eigen::ArrayXd Energies(nb_steps+1);

    // Variables to regulate XYZ output
    double iter_out =
        t_unit / time_step; // Output position every iter_out iterations
    double out_thresh = iter_out; // Threshold to keep track
    std::ofstream trajectory_file(trajectory_filepath);

    std::cout << "Running lj54 sim with timestep: " << time_step << std::endl;
    std::cout << "Storing data in: " << out_dir <<std::endl;

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

    // Open files for data
    std::ofstream energy_out_file(e_filepath);

    // Write headers
    energy_out_file << "Time/(LJ),TotalEnergy/(epsilon)" << std::endl;

    // Write data to the file
    for (int i = 0; i < nb_steps + 1; ++i) {
        energy_out_file << (i)*time_step << "," << Energies(i) << std::endl;
    }

    // Close the data files
    energy_out_file.close();
}

int main(int argc, char *argv[]) {
    std::string time_step_str[] = {"1e-3", "5e-3", "1e-2", "5e-2"};
    double time_step[] = {1e-3, 5e-3, 1e-2, 5e-2};
    for(int i = 0; i < 4; i++){
        simulate_ts(time_step_str[i], time_step[i]);
    }
    return 0;
}
