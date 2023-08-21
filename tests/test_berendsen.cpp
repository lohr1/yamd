//
// Created by robin on 05.07.23.
//

#include "test_berendsen.h"
#include "berendsen.h"
#include "types.h"
#include "helpers.h"

#include "xyz.h"



// Note: Everyone else implemented simple calc test... then Simulation goes in milestones

TEST(BerendsenRescaleTest, Forces) {
    // Checks Berendsen implementation against results of other scaling factor from Eq. 9 of Ch4 Lecture Notes.
    // We expect the resulting velocities to be near.
    auto[names, positions]{read_xyz("/home/robin/School/yamd/xyzs/lattice_4_1.1.xyz")};
    Atoms atoms{positions};
    atoms.velocities.setRandom();
    Velocities_t initial_v = atoms.velocities;

    // Params
    double time_step = 0.001;
    double relax_const = time_step * 1000;
    double target_temp = 300;

    // Rescale velocities with Berendsen implementation
    berendsen_thermostat(atoms, target_temp, time_step, relax_const);

    // Now perform the same calculation but using lambda from other side of Eq. 9
    double mass = 1.;
    Eigen::Array3Xd mv2 = (0.5) * mass * initial_v.square();
    double curr_temp = temp(mv2.sum(), atoms.nb_atoms());
    double lambda = sqrt(target_temp/curr_temp + (1 - target_temp/curr_temp) * exp(-time_step/relax_const));
    // Scale velocities
    Velocities_t scaled_v = lambda * initial_v;

    // Check that magnitudes of resulting velocities are near
    Eigen::ArrayXd mag1 = atoms.velocities.colwise().norm();
    Eigen::ArrayXd mag2 = scaled_v.colwise().norm();

    for(int i = 0; i < atoms.nb_atoms(); i++){
        EXPECT_NEAR(mag1(i), mag2(i), 1e-5);
    }
}