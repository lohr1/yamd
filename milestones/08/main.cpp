/*
 * Here we set up MPI
 */

#include "ducastelle.h"
#include "neighbors.h"
#include "verlet.h"
#include "xyz.h"
#include "helpers.h"
#include "berendsen.h"

#include <fstream>
#include <iostream>

#include "mpi.h"

void propagate(Atoms& atoms, double time_tot_fs, double time_step_fs, double tau_pico,){
    /*
     * Evolves system through time
     */
    // Generate neighbor list
    double cutoff = 10.0; // Achieves approx. 5th nearest neighbor
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    // Propagate system
    double time_tot = time_tot_fs / 10.18; // times 10.18 fs to get real time
    double time_step = time_step_fs /10.18;
    double tau_femto = tau_pico * 1000;
    double tau =
        tau_femto /10.18; // relaxation constant converted to sim units (see top comment)
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
        outputFile << (i) *time_step_fs << "," << Energy(i) << std::endl;
        temp_outputFile << (i) *time_step_fs << "," << Temp(i) << std::endl;
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
    paramFile << "time_step_fs: " << time_step_fs <<  std::endl;
    paramFile << "time_step: " << time_step <<  std::endl;
    paramFile << "tau_pico: " << tau_pico <<  std::endl;
    paramFile << "tau: " << tau <<  std::endl;
    paramFile << "nb_steps: " << nb_steps << std::endl;
    paramFile.close();
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    // Code goes here - MPI comm can be used


    MPI_Finalize();
}
