/*
 * Test file for verlet integrator of single and multiple atoms.
 */
#include "verlet.h"
#include "types.h"
#include <gtest/gtest.h>
#include <vector>

// Single atom test
TEST(VerletTest, SingleAtom) {
    srand((unsigned int) time(0)); // Seed random num generator

    // Initialize pos, vel, and forces for single atom
    // using random doubles in [-1, 1]
    double x = std::rand() * 2.0 / RAND_MAX - 1;
    double y = std::rand() * 2.0 / RAND_MAX - 1;
    double z = std::rand() * 2.0 / RAND_MAX - 1;
    double vx = std::rand() * 2.0 / RAND_MAX - 1;
    double vy = std::rand() * 2.0 / RAND_MAX - 1;
    double vz = std::rand() * 2.0 / RAND_MAX - 1;
    double fx = std::rand() * 2.0 / RAND_MAX - 1;
    double fy = std::rand() * 2.0 / RAND_MAX - 1;
    double fz = std::rand() * 2.0 / RAND_MAX - 1;

    // Save initials for analysis (force is constant)
    std::vector<double> init_pos = {x, y, z};
    std::vector<double> init_vel = {vx, vy, vz};

    // Run verlet integration for nb_steps
    int nb_steps = 1000;
    double timestep = .01;
    for (int i = 0; i < nb_steps; ++i) {
        verlet_step1(x, y, z, vx, vy ,vz, fx, fy, fz, timestep);
        verlet_step2(vx, vy, vz, fx, fy, fz, timestep);
    }
    //Compute pos and vel analytically
    double t = (nb_steps) * timestep;
    // NOTE: Assumes constant force and m=1.
    double anal_x = fx * 0.5 * t*t + init_vel[0] * t + init_pos[0];
    double anal_y = fy * 0.5 * t*t + init_vel[1] * t + init_pos[1];
    double anal_z = fz * 0.5 * t*t + init_vel[2] * t + init_pos[2];

    double anal_vx = fx * t + init_vel[0];
    double anal_vy = fy * t + init_vel[1];
    double anal_vz = fz * t + init_vel[2];

    EXPECT_NEAR(x, anal_x, 1e-6);
    EXPECT_NEAR(y, anal_y, 1e-6);
    EXPECT_NEAR(z, anal_z, 1e-6);
    EXPECT_NEAR(vx, anal_vx, 1e-6);
    EXPECT_NEAR(vy, anal_vy, 1e-6);
    EXPECT_NEAR(vz, anal_vz, 1e-6);
}

// Multiple atom test
TEST(VerletTest, MultipleAtom) {
    // Test for multiple atoms using Eigen Arrays
    Positions_t pos, i_pos;
    Velocities_t vel, i_vel;
    Forces_t frc;
    int num_atoms = 3;

    // Initialize
    pos.setRandom(3, num_atoms);
    vel.setRandom(3, num_atoms);
    frc.setRandom(3, num_atoms);

    // Save initial's
    i_pos = pos;
    i_vel = vel;

    // Integrate
    int nb_steps = 10000;
    double timestep = .1;
    for (int i = 0; i < nb_steps; ++i) {
        verlet_step1(pos, vel, frc, timestep);
        //Calc new forces - constant force
        verlet_step2(vel, frc, timestep);
    }

    // Calc Analytical Solution
    // NOTE: Assumes constant force and m=1.
    double t = (nb_steps) * timestep;
    Positions_t anal_pos = frc * 0.5 * t * t + i_vel * t + i_pos;
    Velocities_t anal_vel = frc * t + i_vel;
    for(int i = 0; i < 3 ; i++){
        for(int j = 0; j < num_atoms; j++){
            EXPECT_NEAR(anal_pos(i, j), pos(i, j), 1e-6);
            EXPECT_NEAR(anal_vel(i, j), vel(i, j), 1e-6);
        }
    }
}