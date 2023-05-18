//
// Created by robin on 08.05.23.
//

#ifndef __VERLET_H
#define __VERLET_H

#include <Eigen/Core>
void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep);
void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep);

void verlet_step1(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                  const Eigen::Array3Xd &forces, double timestep);
void verlet_step2(Eigen::Array3Xd &velocities, const Eigen::Array3Xd &forces, double timestep);

#endif  // __VERLET_H