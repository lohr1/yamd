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
#include <experimental/filesystem>

const double Au_molar_mass = 196.96657; // g/mol, since time is in units of 10.18 fs

void simulate_ts(std::string ts_str, double ts){
    /*
     * Given a time step (ts) in femtoseconds, simulates the equilibrated
     * cluster_923 at 300K to investigate effect on total energy.
     */
    std::string out_dir = "/home/robin/School/HPC/Data/07/ReportData/TimeStep/";
    std::string file_prefix = out_dir + ts_str + "fs";


    // Filepath to store energy data:
    std::string e_filepath = out_dir + ts_str + "_energy.csv";
    // Filepath for trajectory:
    std::string trajectory_filepath = out_dir + ts_str + "_trajectory.xyz";

    // Equilibrated cluster_923 input file:
    std::string xyz_file_path = "/home/robin/School/HPC/Data/07/ReportData/Equilibrated_clusters/923.xyz";

    auto[names, positions, velocities]{read_xyz_with_velocities(xyz_file_path)};
    Atoms atoms{positions, velocities};
    atoms.masses.setConstant(Au_molar_mass);

    // Simulate for 50000 femtoseconds
    double time_fs = 50000;

    int num_frames = 100; // Number of trajectory frames for output

    // Generate neighbor list
    double cutoff = 6.0; // Achieves approx. 5th nearest neighbor
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    // Propagate system
    double time_tot = time_fs / 10.18; // times 10.18 fs to get real time
    double time_step = ts / 10.18;

    int nb_steps = time_tot / time_step;
    int nb_atoms = atoms.nb_atoms();

    // Calc potential energy with ducastelle. Resulting forces stored in atoms.
    double PE = ducastelle(atoms, neighbor_list);

    // Calc KE and total E
    double KE = kinetic_energy(atoms);
    double E = PE + KE;

    // Array to monitor total energy
    Eigen::ArrayXd Energy(nb_steps+1);

    // Variables for XYZ output
    double iter_out = nb_steps / num_frames; // Output position every iter_out iterations
    double out_thresh = iter_out; // Threshold to count iter_out's

    // Trajectory file:
    std::ofstream traj(trajectory_filepath);
    // Write initial frame
    write_xyz(traj, atoms);

    for (int i = 0; i < nb_steps; ++i) {
        // Monitor values
        Energy(i) = E;

        // Verlet predictor step (changes pos and vel)
        verlet_step1(atoms, time_step);

        // Update neighbor_list with new positions (Ok to turn off for stable solid)
        //neighbor_list.update(atoms, cutoff);

        // Update forces with new positions
        PE = ducastelle(atoms, neighbor_list);

        // Verlet step 2 updates velocities, assuming the new forces are present:
        verlet_step2(atoms, time_step);

        // Calc new KE and total E
        KE = kinetic_energy(atoms);
        E = PE + KE;

        // XYZ output
        if(i > out_thresh) {
            std::cout << "ts=" << ts << ", i=" << i << "(/" << nb_steps << ") Writing frame..." << std::endl;
            write_xyz(traj, atoms);
            out_thresh += iter_out;
        }


    }

    // Save final values for monitoring
    Energy(nb_steps) = E;

    traj.close();

    // Open files for data
    std::ofstream outputFile(e_filepath);

    // Write headers
    outputFile << "Time(fs),TotalEnergy" << std::endl;

    // Write data to the file
    for (int i = 0; i < nb_steps+1; ++i) {
        outputFile << (i)*ts << "," << Energy(i) << std::endl;
    }
    outputFile.close();
}


void add_energy(Atoms& atoms, double dQ){
    /*
     * Adds dQ energy units to atoms by rescaling the velocities by
     * a factor of lambda = sqrt((KE + dQ)/KE)
     */
    std::cout << "Adding energy\n";
    double KE = kinetic_energy(atoms);
    double lambda = sqrt((KE + dQ)/KE);
    atoms.velocities *= lambda;
}

std::pair<double, double> relax(Atoms& atoms, NeighborList& neighbor_list, std::string traj_dir, std::string cluster_name, double cutoff, double time_relax, double ts, int iteration){
    /*
     * Relaxes cluster for time_relax. Puts trajectory in traj_dir with relax_iteration number in filename.
     * Returns pair containing (average Energy, average Temp) over final 500 iterations of this relaxation.
     */
    int nb_steps = time_relax / ts;
    double time_step = ts / 10.18;
    double PE, KE, E;
    std::string trajectory_file = traj_dir + std::to_string(iteration) + "_trajectory.xyz";
    std::ofstream traj(trajectory_file);

    int num_frames = 100;
    double iter_out = time_relax / num_frames;
    double out_thresh = iter_out;

    // Arrays for computing averages towards end of relaxation
    int avg_over_iter = 500;
    Eigen::ArrayXd Energy(avg_over_iter);
    Eigen::ArrayXd Temp(avg_over_iter);

    for (int i = 0; i < nb_steps; ++i) {
        // Verlet predictor step (changes pos and vel)
        verlet_step1(atoms, time_step);

        // Update neighbor_list with new positions (Ok to turn off for stable solid)
        neighbor_list.update(atoms, cutoff);

        // Update forces with new positions
        PE = ducastelle(atoms, neighbor_list);

        // Verlet step 2 updates velocities, assuming the new forces are present:
        verlet_step2(atoms, time_step);

        // XYZ output
        if(i > out_thresh) {
            std::cout << cluster_name << " is relaxing..." << " at " << "i=" << i << "(/" << nb_steps << std::endl;
            std::cout << "Writing to " << trajectory_file << std::endl;
            write_xyz(traj, atoms);
            out_thresh += iter_out;
        }

        if(nb_steps - i <= avg_over_iter){
            // Begin taking data to compute averages
            KE = kinetic_energy(atoms);
            E = PE + KE;

            int index = nb_steps - i - avg_over_iter;
            index *= -1;
            Energy(index) = E;
            Temp(index) = temp(KE, atoms.nb_atoms());
        }

    }
    double avg_E = Energy.sum() / avg_over_iter;
    double avg_T = Temp.sum() / avg_over_iter;
    return std::make_pair(avg_E, avg_T);
}

void energy_vs_temp(std::string cluster_name, int iterations, double dQ, double time_relax, double ts){
    /*
     * Generates iterations number of (Energy, Temp) data points for the given
     * cluster using given dQ
     */
    // Initialize atoms
    std::string xyz_filepath = "/home/robin/School/HPC/Data/07/ReportData/Equilibrated_clusters/" + cluster_name + ".xyz";
    //std::string xyz_filepath = "/home/robin/School/HPC/Data/07/ReportData/Adding_Energy2/cluster_3871/final_state.xyz";

    auto[names, positions, velocities]{read_xyz_with_velocities(xyz_filepath)};
    Atoms atoms{positions, velocities};
    atoms.masses.setConstant(Au_molar_mass);

    // Initialize arrays for data
    Eigen::ArrayXd Energy(iterations+1);
    Eigen::ArrayXd Temp(iterations+1);

    // Generate neighbor list
    double cutoff = 6.0;
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    // Generate initial data point
    double PE = ducastelle(atoms, neighbor_list);
    double KE = kinetic_energy(atoms);
    double E = PE + KE;
    double T = temp(KE, atoms.nb_atoms());

    std::string out_dir = "/home/robin/School/HPC/Data/07/ReportData/Adding_Energy2/" + cluster_name + "/";

    // Directory for trajectory files
    std::string traj_dir = out_dir + "trajectories/";

    // Energy vs. Temp file
    std::ofstream E_vs_T_file(out_dir + "energy_temp.csv");

    std::pair<double, double> ET_point;
    for(int i = 0; i < iterations; i++){
         // Store data from before
        Energy(i) = E;
        Temp(i) = T;

        // Add energy
        add_energy(atoms, dQ);

        // Let relax (writes new traj file, returns pair(avgE, avgT))
        ET_point = relax(atoms,neighbor_list,traj_dir,cluster_name,cutoff,time_relax,ts,i);
        std::cout << "Completed an iteration. (Energy, Temp): (" << ET_point.first << ", "
                  << ET_point.second << ")" << std::endl;
        E = ET_point.first;
        T = ET_point.second;
    }
    // Store final data point
    Energy(iterations) = E;
    Temp(iterations) = T;

    // Write headers
    E_vs_T_file << "Temp(K),TotalEnergy" << std::endl;


    // Write data to the file
    for (int i = 0; i < iterations+1; ++i) {
        E_vs_T_file << Temp(i) << "," << Energy(i) << std::endl;
    }

    // Close the data files
    E_vs_T_file.close();

    // Write final configuration to XYZ file:
    std::ofstream final_xyz_file(out_dir + "final_state.xyz");
    write_xyz(final_xyz_file, atoms);
    final_xyz_file.close();
}


int main(int argc, char *argv[]) {
    // Test time steps on equilibrated (300K) cluster_923
//    std::string time_steps[] = {"1", "2", "5", "10"}; // Timesteps in fs
//    double time_step[] = {1., 2., 5., 10.};
//    for(int i = 0; i < 4; i++){
//        simulate_ts(time_steps[i], time_step[i]);
//    }
//    simulate_ts("20", 20.);

    std::string cluster_names[] = {"cluster_923", "cluster_3871", "cluster_6525", "cluster_10179"};
    double ts = 5; // time step in femtoseconds

    double time_fs = 5000; // Time for system to relax after energy added

    int iterations = 27;
    double dQ = 400;
    energy_vs_temp(cluster_names[3],iterations,dQ,time_fs,ts);


    // Code for equilibrating several clusters
//    double time_fs = 10000;
//    double ts = 5;
//    double target_temp = 300;
//    double tau_fs = 100;
//    // To equilibrate clusters, run this until satisfactory equilibration:
//    for(int i = 0; i < 4; i++){
//        equilibrate_cluster(cluster_names[i],time_fs,ts,target_temp,tau_fs);
//    }


    // Now test equilibration by running without thermostat

//    double t_therm_off = sim_time_eq; // Times 10.18 fs for real time
//    test_eq(atoms,t_therm_off,real_time_step, out_dir);
    

    return 0;
}
