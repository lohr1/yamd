/*
 * This file contains functions for verlet integrations steps of
 * single atoms and multiple atoms.
 */
#include "verlet.h"
#include "types.h"

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep) {
    // Predictor step - mass is assumed to be 1.
    vx += 0.5 * fx * timestep;
    vy += 0.5 * fy * timestep;
    vz += 0.5 * fz * timestep;
    x += vx * timestep;
    y += vy * timestep;
    z += vz * timestep;
}

void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep) {
    // Assuming new forces in params
    vx += 0.5 * fx * timestep;
    vy += 0.5 * fy * timestep;
    vz += 0.5 * fz * timestep;
}

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double timestep){
    // Predictor step
    velocities += 0.5 * forces * timestep;
    positions += velocities * timestep;
}

void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double timestep){
    // Assuming new forces in params
    velocities += 0.5 * forces * timestep;
}

void verlet_step1(Atoms &atoms, double timestep){
    verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);
}
void verlet_step2(Atoms &atoms, double timestep){
    verlet_step2(atoms.velocities, atoms.forces, timestep);
}
