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
#include "mpi_support.h"
#include "domain.h"

const double Au_molar_mass = 196.96657; // g/mol, since time is in units of 10.18 fs

inline void propagate_domain(Domain& domain, Atoms& atoms, double time_tot_fs, double time_step_fs, const std::string& dir){
    /*
     * Evolves system through time, making necessary adjustments for domain decomposition.
     *
     * Domain decomp should already be enabled before calling, and should be
     * disabled after calling.
     */
    // Generate neighbor list
    double cutoff = 10.0; // Achieves approx. 5th nearest neighbor
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    // Convert time variables to sim units
    double time_tot = time_tot_fs / 10.18; // times 10.18 fs to get real time
    double time_step = time_step_fs /10.18;
    int nb_steps = time_tot / time_step;

    // Update ghosts, saving number of local atoms
    int nb_local = atoms.nb_atoms();
    double border_width = 2 * cutoff;
    domain.update_ghosts(atoms, border_width);


    // Calc forces with ghost atoms, but save only PE of local atoms:
    double PE_local = ducastelle_domain_decomp(atoms, nb_local, neighbor_list);


    // Calc KE of local atoms
    double KE_local = kinetic_energy_local(atoms, nb_local);

    double E_local = PE_local + KE_local;


    // Calculate total energy with MPI
    double E_total = MPI::allreduce(E_local, MPI_SUM, MPI_COMM_WORLD);


    // Array to monitor total energy
    Eigen::ArrayXd Energy(nb_steps+1);

    // Process of rank 0 will handle file IO
    int rank = MPI::comm_rank(MPI_COMM_WORLD);

    // Variables for keeping track of when to save a frame
    double iter_out, out_thresh;
    std::ofstream trajectory_file;

    // Create trajectory output file
    // Recover replicated state:
    domain.disable(atoms);
    if(rank == 0) {
        iter_out =
            time_tot / 100; // Output position every iter_out iterations
        out_thresh = iter_out; // Threshold to count iter_out's

        // Create trajectory file
        trajectory_file = std::ofstream(dir + "trajectory.xyz");
        // Write initial frame
        write_xyz(trajectory_file, atoms);
    }

    // Go back into decomposed state and repopulate ghosts
    domain.enable(atoms);
    domain.update_ghosts(atoms, border_width);

    Eigen::ArrayXd pos_before;


    for (int i = 0; i < nb_steps; ++i) {

        // Monitor total energy
        Energy(i) = E_total;

        // Verlet predictor step (changes pos and vel)
        verlet_step1(atoms, time_step);

//        int frame = 3;
//        if(i== frame - 1){
//            pos_before = atoms.positions.row(0);
//        }
//        if(i == frame){
//            Eigen::ArrayXd x_vels = atoms.positions.row(0);
//            for(int j = 0; j < nb_local; j++){
//                std::cout << "Atom: " << j
//                          << " Delta X: " << x_vels(j) - pos_before(j) << std::endl;
//            }
//            break;
//        }


        // Update forces with new positions
        PE_local = ducastelle_domain_decomp(atoms, nb_local, neighbor_list);


        // Verlet step 2 updates velocities, assuming the new forces are present:
        verlet_step2(atoms, time_step);

//        Eigen::ArrayXd x_vels = atoms.velocities.row(0);
//        for(int i = 0; i < nb_local; i++){
//            if(x_vels(i) > 400)
//                std::cout << "Atom: " << i << " X-Velocity: " << x_vels << std::endl;
//        }

        // Calc new KE_local and E_local
        KE_local = kinetic_energy_local(atoms, nb_local);
        E_local = PE_local + KE_local;


        // Calc new E_total
        E_total = MPI::allreduce(E_local, MPI_SUM, MPI_COMM_WORLD);

//        // XYZ output
        if(i > out_thresh) {
            // Write a frame
            domain.disable(atoms);
            if(rank == 0) {
                write_xyz(trajectory_file, atoms);

                // Increment the output threshold counter
                out_thresh += iter_out;
            }
            domain.enable(atoms);
            domain.update_ghosts(atoms, border_width);
        }
        // Exchange atoms between subdomains after each step (Removes ghosts)
        domain.exchange_atoms(atoms);

        // Update num local atoms
        nb_local = atoms.nb_atoms();
        // Repopulate ghosts
        domain.update_ghosts(atoms,border_width);


        // Update neighborlist
        neighbor_list.update(atoms,cutoff);

    }

    std::cout << "Left loop" << std::endl;

    // Save final energy
    Energy(nb_steps) = E_total;

    if(rank == 0) {
        trajectory_file.close();

        // Open files for data
        std::ofstream outputFile(dir + "energy_data.csv");

        // Write headers
        outputFile << "Time(fs),TotalEnergy" << std::endl;

        // Write data to the file
        for (int i = 0; i < nb_steps + 1; ++i) {
            outputFile << (i)*time_step_fs << "," << Energy(i) << std::endl;
        }

        // Close the data files
        outputFile.close();

        // Write final configuration to XYZ file
        // in case we want to start where we left off
        std::ofstream final_xyz_file(dir + "final_state.xyz");
        domain.disable(atoms);
        write_xyz(final_xyz_file, atoms);
        final_xyz_file.close();
        domain.enable(atoms);

        // Write params to a text file
        std::ofstream paramFile(dir + "params.txt");
        paramFile << "time_tot: " << time_tot << std::endl;
        paramFile << "time_step_fs: " << time_step_fs << std::endl;
        paramFile << "time_step: " << time_step << std::endl;
        paramFile << "nb_steps: " << nb_steps << std::endl;
        paramFile.close();
    }
}

int main(int argc, char *argv[]) {
    // Initial setup which doesn't require MPI:
    std::string xyz_file = "/home/robin/School/yamd/xyzs/cluster_3871.xyz";
    std::string out_dir = "/home/robin/School/HPC/Data/08/DomainDecomp/cluster_3871/";

    // Time parameters in real units
    double time_fs = 1e5;
    double timestep_fs = 1;

    // Init atoms
//    auto[names, positions, velocities]{read_xyz_with_velocities(xyz_file)};
//    Atoms atoms{positions, velocities};
//    atoms.masses.setConstant(Au_molar_mass);

    auto[names, positions]{read_xyz(xyz_file)};
    Atoms atoms{positions};
    atoms.masses.setConstant(Au_molar_mass);

    double cutoff = 10.0;


    // Define domain dimensions with atoms positions extrema plus wiggle room to keep cluster from interacting with itself
    Eigen::Array3d domain_lengths(3);
    for(int i = 0; i < 3; i++){
        domain_lengths(i) = atoms.positions.row(i).maxCoeff() - atoms.positions.row(i).minCoeff() + 10.1;
    }

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Init domain
    Domain domain(MPI_COMM_WORLD, domain_lengths, {1, 1, 1}, {1, 1, 1});

    // Decompose atoms into subdomains
    domain.enable(atoms);

    // Run simulation
    propagate_domain(domain, atoms, time_fs, timestep_fs, out_dir);

    domain.disable(atoms);
    MPI_Finalize();
    return 0;
}
