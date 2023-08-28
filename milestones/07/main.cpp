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

const double Au_molar_mass = 196.96657; // g/mol, since time is in units of 10.18 fs


void gen_ET_point(Atoms& atoms, double dQ, const std::string& dir){
    /*
     * Adds dQ eV to atoms by rescaling velocities.
     * Lets system relax and returns average temperature.
     */
    int nb_atoms = atoms.nb_atoms();
    double dQ_per_atom = dQ / nb_atoms;  // Energy to add to each atom

}

void test_eq(Atoms& atoms, double t_tot, double real_ts, const std::string& dir){
    /*
     * Runs simulation with thermostat turned off and saves
     * temp data/trajectory.
     */
    // Generate neighbor list
    double cutoff = 10.0; // Achieves approx. 5th nearest neighbor
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    // Propagate system
    double time_tot = t_tot; // times 10.18 fs to get real time
    double real_time_step = real_ts; // in fs
    double time_step = real_time_step/10.18;
    int nb_steps = time_tot / time_step;
    int nb_atoms = atoms.nb_atoms();

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
    std::ofstream traj(dir + "therm_off_trajectory.xyz");
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

        // Calc new KE and total E
        KE = kinetic_energy(atoms);
        E = PE + KE;

        // XYZ output
        if(i > out_thresh) {
            write_xyz(traj, atoms);
            out_thresh += iter_out;
        }
        
        // Update neighbor_list with new positions
        neighbor_list.update(atoms,cutoff);
    }

    // Save final values for monitoring
    Energy(nb_steps) = E;
    Temp(nb_steps) = temp(KE, nb_atoms);

    traj.close();

    // Open files for data
    std::ofstream outputFile(dir + "therm_off_energy_data.csv");
    std::ofstream temp_outputFile(dir + "therm_off_temp_data.csv");


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
}


void equilibrate_EAM(Atoms& atoms, double target_temp, double time_eq, double real_ts, double tau_eq_pico, const std::string& dir){
    /*
     * Brings atoms object to target temp (EAM potential).
     * Writes trajectory, data, and param files to dir.
     *
     * Note: temp in K, time in 10.18 fs, mass in g/mol
     */
    // Generate neighbor list
    double cutoff = 10.0; // Achieves approx. 5th nearest neighbor
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    // Propagate system    
    double time_tot = time_eq; // times 10.18 fs to get real time
    double real_time_step = real_ts; // in fs
    double time_step = real_time_step/10.18;
    double real_tau_pico = tau_eq_pico; // relaxation constant in picoseconds
    double real_tau_femto = real_tau_pico * 1000;
    double tau = real_tau_femto/10.18; // relaxation constant converted to sim units (see top comment)
    int nb_steps = time_tot / time_step;
    int nb_atoms = atoms.nb_atoms();

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
    
    // Write params to a text file
    std::ofstream paramFile(dir + "params.txt");
    paramFile << "time_tot: " << time_tot <<  std::endl;
    paramFile << "real_time_step: " << real_time_step <<  std::endl;
    paramFile << "time_step: " << time_step <<  std::endl;
    paramFile << "real_tau_pico: " << real_tau_pico <<  std::endl;
    paramFile << "tau: " << tau <<  std::endl;
    paramFile << "nb_steps: " << nb_steps << std::endl;
    paramFile.close();
}


int main(int argc, char *argv[]) {
    // Path to xyz file for this run:
    std::string xyz_file = "/home/robin/School/HPC/Data/07/Equilibration/cluster_3871/final_state.xyz";
    // Directory to store data from this run:
    std::string dir = "/home/robin/School/HPC/Data/07/Equilibration/cluster_3871/longer_eq/";

    double time_eq = 800; // Times 10.18 fs for real time
    double tau_eq_pico = 1; // Times 1000 for fs
    double real_time_step = 1; // fs
    double target_temp = 300; // K (for equilibration only)


    // Time for cluster to relax after energy increase
    // Also time over which avg temp is calculated
    double real_tau_relax = 200; 


    // For initial cluster:

//     auto[names, positions]{read_xyz(xyz_file)};
//     Atoms atoms{positions};
//     atoms.velocities.setRandom();
//     atoms.velocities *= 0.4;
//     atoms.masses.setConstant(Au_molar_mass);

     // Otherwise:

    auto[names, positions, velocities]{read_xyz_with_velocities(xyz_file)};
    Atoms atoms{positions, velocities};
    atoms.masses.setConstant(Au_molar_mass);


    equilibrate_EAM(atoms,target_temp, time_eq, real_time_step, tau_eq_pico, dir);

    // Now test equilibration by running without thermostat

    double t_therm_off = 300; // Times 10.18 fs for real time
    test_eq(atoms,t_therm_off,real_time_step, Au_molar_mass, dir);
    return 0;
}
