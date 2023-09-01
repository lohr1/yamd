/*
 * Investigation of Stress vs. Strain on gold nano-wires
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

double volume(Atoms& atoms, double side_length_shared) {
    /*
     * Calculate rough volume of whisker by first calculating rectangular volume
     * containing whisker, then subtracting the appropriate volumes (the "holes"
     * whose cross-sections (orthogonal to the length) are triangles).
     *
     * Param side_length_shared is the length of the side of the whisker which
     * is concurrent with the long side of the rectangle in each cross-section.
     *
     * Note: this method assumes the whisker is oriented correctly with
     * respect to the x, y, and z axes, which can be easily seen in ovito.
     */
    Eigen::Array3d side_lengths(3);
    for (int i = 0; i < 3; i++) {
        side_lengths(i) = atoms.positions.row(i).maxCoeff() -
                          atoms.positions.row(i).minCoeff();
    }
    double prism_volume = side_lengths.prod();
    // Now find the dimensions of the triangles in each cross section
    double triangle_base, triangle_height;
    // Long side of rectangular cross section is shorter than length of prism and
    // longer than short side of cross section:
    std::sort(side_lengths.begin(), side_lengths.end());
    double rect_cross_sec_length = side_lengths(1);
    triangle_base = 0.5 * (rect_cross_sec_length - side_length_shared);
    triangle_height = 0.5 * side_lengths(0);
    double triangle_area = 0.5 * triangle_base * triangle_height;
    return prism_volume - 4 * (triangle_area * side_lengths(2));
}

Eigen::Array33d stress_tensor(Atoms& atoms, double side_length_shared){
    // Note: contribution from temperature is negligible
    // FORCES/POS MUST BE UP TO DATE!!!
    int nb_atoms = atoms.nb_atoms();

    // Sum of outer products for all i < j
    Eigen::Matrix3Xd sum_outer_prod(3, 3);
    Eigen::VectorXd r_ij(3), f_ij(3);
    sum_outer_prod.setZero();
    for(int i = 0; i < nb_atoms - 1; i++){
        for(int j = i + 1; j < nb_atoms; j++){
            r_ij = atoms.positions.col(i) - atoms.positions.col(j);
            f_ij = atoms.forces.col(i) - atoms.forces.col(j);
            sum_outer_prod += r_ij * f_ij.transpose();
        }
    }
    double volume = volume(atoms, side_length_shared);
    Eigen::Array33d cauchy_stress = (1 / volume) * sum_outer_prod;
    return cauchy_stress;
}

void stretch_whisker(Domain& domain, Atoms& atoms, double side_length_shared, double initial_length, double cutoff, double time_tot_fs, double time_step_fs, const std::string& dir){
    /*
     * Evolves system through time, making necessary adjustments for domain decomposition.
     *
     *
     *
     * Domain decomp should be enabled before calling, and should be
     * disabled after calling.
     */

    // For XYZ output
    int num_frames = 500;

    // Convert time variables to sim units
    double time_tot = time_tot_fs / 10.18; // times 10.18 fs to get real time
    double time_step = time_step_fs /10.18;
    int nb_steps = time_tot / time_step;

    // Update ghosts, saving number of local atoms
    int nb_local = atoms.nb_atoms();
    double border_width = 2 * cutoff;
    domain.update_ghosts(atoms, border_width);

    // Create neighbor list including ghosts
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

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

    iter_out =
        nb_steps / num_frames; // Output position every iter_out iterations
    out_thresh = iter_out; // Threshold to count iter_out's

    // Create trajectory output file
    // Recover replicated state:
    domain.disable(atoms);
    if(rank == 0) {
        // Create trajectory file
        trajectory_file = std::ofstream(dir + "trajectory.xyz");
        // Write initial frame
        write_xyz(trajectory_file, atoms);
    }

    // Go back into decomposed state and repopulate ghosts
    domain.enable(atoms);
    domain.update_ghosts(atoms, border_width);

    for (int i = 0; i < nb_steps; ++i) {

        // Monitor total energy
        Energy(i) = E_total;

        // Verlet predictor step (changes pos and vel)
        verlet_step1(atoms, time_step);

        // Pos changed, so exchange atoms (Removes ghosts)
        domain.exchange_atoms(atoms);

        // Update num local atoms
        nb_local = atoms.nb_atoms();

        // Repopulate ghosts for PE calc
        domain.update_ghosts(atoms, border_width);

        // Update neighbor_list with new ghosts
        neighbor_list.update(atoms,cutoff);

        // Calc PE and forces with new neighbor_list, saving just PE for atoms in subdomain (nb_local)
        PE_local = ducastelle_domain_decomp(atoms, nb_local, neighbor_list);

        // Update velocities according to new forces:
        verlet_step2(atoms, time_step);

        // Calc new KE_local and E_local
        KE_local = kinetic_energy_local(atoms, nb_local);
        E_local = PE_local + KE_local;


        // Calc new E_total
        E_total = MPI::allreduce(E_local, MPI_SUM, MPI_COMM_WORLD);

        // XYZ output
        if(i > out_thresh) {
            // Write a frame
            domain.disable(atoms);
            if(rank == 0) {
                std::cout << "i= " << i << ". writing frame" << std::endl;
                write_xyz(trajectory_file, atoms);
            }
            // Increment the output threshold counter
            out_thresh += iter_out;

            domain.enable(atoms); // Ghosts removed (They are updated properly in next iteration)
        }
    }

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

    // Number of processors
    int p = 4;
    // Files:
    std::string xyz_file = "/home/robin/School/yamd/xyzs/whisker_small.xyz";
    std::string out_dir = "/home/robin/School/HPC/Data/09/whisker_small/Domain_fromXYZ/";

    // For calculating volume of whisker

    // Time parameters in real units
    double time_fs = 7e4;
    double timestep_fs = 1;

    // Cutoff for neighbor_list
    double cutoff = 6.0;

    // Init atoms
//    auto[names, positions, velocities]{read_xyz_with_velocities(xyz_file)};
//    Atoms atoms{positions, velocities};
//    atoms.masses.setConstant(Au_molar_mass);

    auto[names, positions]{read_xyz(xyz_file)};
    Atoms atoms{positions};
    atoms.masses.setConstant(Au_molar_mass);

//    // Set positions data to start from 0 in each dimension
//    double min;
//    for(int d = 0; d < 3 ; d++){
//        min = atoms.positions.row(d).minCoeff();
//        atoms.positions.row(d) -= min;
//    }

    // Define domain dimensions
    double wiggle_room = 20.0;  // For x and y dimensions only (whisker should continue in z dimension)
//    Eigen::Array3d domain_lengths(3);
//    for(int i = 0; i < 3; i++){
//        domain_lengths(i) =
//            atoms.positions.row(i).maxCoeff() - atoms.positions.row(i).minCoeff();
//        domain_lengths(i) = atoms.positions.row(i).maxCoeff();
//        if(i!=2)
//            domain_lengths(i) += wiggle_room;
//    }

//    // Center xy positions in domain
//    double shift = wiggle_room / 2;
//    for(int j = 0; j < 2; j++){
//        atoms.positions.row(j) += shift;
//    }

    Eigen::Array3d domain_lengths{40.38993934, 40.79999999, 144.24978336};

    // Calculate initial length of whisker
    double initial_length = domain_lengths(2);

    std::cout << "Initial length: " << initial_length << std::endl;

    std::cout << "Stretching whisker : " <<
        xyz_file << "\nStoring data in: " << out_dir <<
        "\nInitializing MPI..." << std::endl;


    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Init domain
    Domain domain(MPI_COMM_WORLD, domain_lengths,
                  {1, 1, p}, {0, 0, 1});

    // Decompose atoms into subdomains
    domain.enable(atoms);

    // Run simulation
    stretch_whisker(domain, atoms, initial_length, cutoff, time_fs, timestep_fs, out_dir);

    domain.disable(atoms);
    MPI_Finalize();
    return 0;
}
