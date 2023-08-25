/*
 * Equilibrates atoms object to target temp (EAM potential).
 * Writes trajectory and data files to dir.
 */

#include "atoms.h"
#include "berendsen.h"
#include "neighbors.h"
#include "ducastelle.h"
#include "helpers.h"
#include "verlet.h"
#include "xyz.h"

#include <fstream>
#include <iostream>

void equilibrate(Atoms& atoms, double target_temp, const std::string& dir){

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
            outputFile << (i+1) * real_time_step << "," << Energy(i) << std::endl;
            temp_outputFile << (i+1) * real_time_step << "," << Temp(i) << std::endl;
        }

        // Close the data files
        outputFile.close();
        temp_outputFile.close();

        // Write final configuration to XYZ file just in case:
        std::ofstream final_xyz_file(dir + "final_state.xyz");
        write_xyz(final_xyz_file, atoms);
        final_xyz_file.close();
}