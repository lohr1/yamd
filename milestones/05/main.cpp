
#include "lj_direct_summation.h"
#include "verlet.h"
#include "xyz.h"
#include "helpers.h"
#include "gnuplot-iostream.h"

#include <fstream>
#include <iostream>

int main(int argc, char *argv[]) {
    auto[names, positions]{read_xyz("/home/robin/School/yamd/milestones/05/lattice.xyz")};

    Atoms atoms{positions};
    atoms.velocities.setRandom();

    // Integrate
    double mass = 1.0;
    double epsilon = 1.0;
    double sigma = 1.0;

    double expr = sqrt(mass * sigma * sigma / epsilon);
    double time_tot = 100 * expr;
    double time_step = 0.001 * expr;
    int nb_steps = time_tot / time_step;

    // Calc forces, store in atoms. Also calc energies.
    double PE = lj_direct_summation(atoms, epsilon, sigma);
    double KE = kinetic_energy(atoms);

    // To monitor Energies (First index potential)
    Eigen::Array2Xd Energies(2, nb_steps);


    // Variables for XYZ output
    double out_thresh = expr;
    std::ofstream traj("lattice_traj.xyz");

    for (int i = 0; i < nb_steps; ++i) {

        // Verlet predictor step (changes pos and vel)
        verlet_step1(atoms, time_step);

        // Update forces with new positions
        lj_direct_summation(atoms, epsilon, sigma);

        // Verlet step 2 updates velocities assumes the new forces are present in params:
        verlet_step2(atoms, time_step);

        // XYZ output
        if(i > out_thresh) {
            write_xyz(traj, atoms);
            out_thresh += expr;
        }
    }

    traj.close();

//    // Writing Energies to file
//    std::ofstream outfile("Energies.txt");
//    if (outfile.is_open()) {
//        outfile << Energies << std::endl;
//        double tot_diff = Energies(0, 0) - Energies(0, Energies.cols()-1);
//        outfile << "Total Energy difference: " << tot_diff << std::endl;
//        std::cout << "Total Energy difference: " << tot_diff << std::endl;
//        std::cout << "Made it." << std::endl;
//        outfile.close();
//    } else {
//        std::cerr << "Error opening the file." << std::endl;
//    }
    return 0;
}
