/*
 * This program was used to equilibrate Mackay icosahedra and then test the effect
 * of different time steps on the total energy.
 *
 * Simulation units:
 * Length: Angstroms
 * Energy: eV
 * Mass: g_mol
 *
 * -> 1 Time unit in code is 10.18 fs
 */

#include "ducastelle.h"
#include "neighbors.h"
#include "verlet.h"
#include "xyz.h"
#include "helpers.h"
#include "berendsen.h"

#include <fstream>
#include <iostream>



void equilibrate_EAM(Atoms& atoms, double target_temp, const std::string& dir){

    // Generate neighbor list
    double cutoff = 10.0; // Achieves approx. 5th nearest neighbor
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    // Propagate system
    double time_tot = 200; // times 10.18 fs to get real time
    double real_time_step = 1; // in fs
    double time_step = real_time_step/10.18;
    double real_relax_const_pico = 0.1; // relaxation constant in picoseconds
    double tau = (real_relax_const_pico * 1000)/10.18; // relaxation constant converted to sim units (see top comment)
    int nb_steps = time_tot / time_step;
    int nb_atoms = atoms.nb_atoms();

    std::cout << atoms.velocities.col(0);

    // Calc potential energy with ducastelle. Resulting forces stored in atoms.
    double PE = ducastelle(atoms, neighbor_list);

    // Calc KE and total E
    double KE = kinetic_energy(atoms);
    double E = PE + KE;

    // Array to monitor total energy
    Eigen::ArrayXd Energy(nb_steps+1);
    Eigen::ArrayXd Temp (nb_steps+1);

    // Variables for XYZ output
    double iter_out = time_tot / 100; // Output position every iter_out iterations
    double out_thresh = iter_out; // Threshold to count iter_out's

    // Trajectory file:
    std::ofstream traj(dir + "trajectory.xyz");
    // Write initial frame
    write_xyz(traj, atoms);

    for (int i = 0; i < nb_steps; ++i) {
        // Monitor values
        Energy(i) = E;
        Temp(i) = temp(KE, nb_atoms);

        // Verlet predictor step (changes pos and vel)
        verlet_step1(atoms, time_step);

        // Update forces with new positions
        PE = ducastelle(atoms, neighbor_list);

        // Verlet step 2 updates velocities, assuming the new forces are present:
        verlet_step2(atoms, time_step);

        // Berendsen velocity rescaling
        berendsen_thermostat(atoms, target_temp, time_step, tau);

        // Calc new KE and total E
        KE = kinetic_energy(atoms);
        E = PE + KE;

        // XYZ output
        if(i > out_thresh) {
            write_xyz(traj, atoms);
            out_thresh += iter_out;
        }
    }

    // Save final values for monitoring
    Energy(nb_steps) = E;
    Temp(nb_steps) = temp(KE, nb_atoms);

    traj.close();

    // Open files for data
    std::ofstream outputFile(dir + "energy_data.csv");
    std::ofstream temp_outputFile(dir + "temp_data.csv");


    // Write headers
    outputFile << "Time(fs),TotalEnergy" << std::endl;
    temp_outputFile << "Time(fs),Temp(K)" << std::endl;


    // Write data to the file
    for (int i = 0; i < nb_steps+1; ++i) {
        outputFile << (i) * real_time_step << "," << Energy(i) << std::endl;
        temp_outputFile << (i) * real_time_step << "," << Temp(i) << std::endl;
    }

    // Close the data files
    outputFile.close();
    temp_outputFile.close();

    // Write final configuration to XYZ file:
    std::ofstream final_xyz_file(dir + "final_state.xyz");
    write_xyz(final_xyz_file, atoms);
    final_xyz_file.close();
}

int main(int argc, char *argv[]) {
    // Path to xyz file for this run:
    std::string xyz_file = "/home/robin/School/HPC/Data/07/Equilibration/cluster_3871/final_state.xyz";
    // Directory to store data from this run:
    std::string dir = "/home/robin/School/HPC/Data/07/Equilibration/cluster_3871/2nd(1fs)/";

    // There are 2 starting cases - with or without velocities.
    // So I use an if block and must therefore declare these var's outside.
    std::vector<std::string> names;
    Eigen::Array3Xd positions;
    Eigen::Array3Xd velocities;

    if (xyz_file.find("cluster_") != std::string::npos) {
        // We are working with an initial icosahedron (no velocities)
        std::tie(names, positions) = read_xyz(xyz_file);
    }
    else{
        // We are "warm" starting - so read velocities
        std::tie(names, positions, velocities) = read_xyz_with_velocities(xyz_file);
    }

    if(positions.cols() != velocities.cols()){
        velocities = Eigen::Array3Xd(3, positions.cols());
    }

    // Initialize atoms pos and vel
    Atoms atoms{positions, velocities};

    // Target temperature for equilibration
    double target_temp = 300;

    equilibrate_EAM(atoms,target_temp,dir);

    return 0;
}
