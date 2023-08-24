#include "ducastelle.h"
#include "neighbors.h"
#include "verlet.h"
#include "xyz.h"
#include "helpers.h"
#include "berendsen.h"

#include <fstream>
#include <iostream>

int main(int argc, char *argv[]) {
    // Using clusters provided, determine good time step for EAM potential
    // by monitoring total energy during simulation
    std::string xyz_file = "/home/robin/School/yamd/xyzs/cluster_923.xyz";

    // Directory to store data from run:
    std::string dir = "/home/robin/School/HPC/Data/07/Determine_timestep/no_initial_v/1fs/";
    auto[names, positions]{read_xyz(xyz_file)};

    // Initialize atoms pos and vel
    Atoms atoms{positions};
    //atoms.velocities.setRandom();  // units of Angstrom / 10.18 fs

    // Generate neighbor list
    double cutoff = 10.0; // Achieves approx. 5th nearest neighbor
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    // Propagate system
    double time_tot = 500; // times 10.18 fs to get real time
    double real_time_step = 1; // fs
    double time_step = real_time_step/10.18;
    //double relax_const = 1 * expr;  // No thermostat here as we want to conserve total energy
    int nb_steps = time_tot / time_step;
    int N = atoms.nb_atoms();

    // Calc potential energy with ducastelle. Resulting forces stored in atoms.
    double PE = ducastelle(atoms, neighbor_list);

    // Calc KE and total E
    double KE = kinetic_energy(atoms);
    double E = PE + KE;

    // Array to monitor total energy
    Eigen::ArrayXd Energy(nb_steps+1);


    // Variables for XYZ output
    double iter_out = time_tot / 100; // Output position every iter_out iterations
    double out_thresh = iter_out; // Threshold to count iter_out's

    // Trajectory file:
    std::ofstream traj(dir + "trajectory.xyz");
    // First frame
    write_xyz(traj, atoms);

    for (int i = 0; i < nb_steps; ++i) {
        // Monitor Energies
        Energy(i) = E;

        // Verlet predictor step (changes pos and vel)
        verlet_step1(atoms, time_step);

        // Update forces with new positions
        PE = ducastelle(atoms, neighbor_list);

        // Verlet step 2 updates velocities, assuming the new forces are present:
        verlet_step2(atoms, time_step);

        // Berendsen velocity rescaling - Not this time: we want to conserve total energy
        //berendsen_thermostat(atoms, 300, time_step, relax_const);

        // Calc new KE and total E
        KE = kinetic_energy(atoms);
        E = PE + KE;

        // XYZ output
        if(i > out_thresh) {
            write_xyz(traj, atoms);
            out_thresh += iter_out;
        }
    }

    // Save final E
    Energy(nb_steps) = E;

    traj.close();

    // Open a file for writing
    std::ofstream outputFile(dir + "energy_data.csv");

    // Write header
    outputFile << "Time(fs),TotalEnergy" << std::endl;

    // Write data to the file
    for (int i = 0; i <= nb_steps; ++i) {
        outputFile << i * real_time_step << "," << Energy(i) << std::endl;
    }

    // Close the file
    outputFile.close();

    return 0;
}
