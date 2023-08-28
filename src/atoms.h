//
// Created by robin on 23.05.23.
//

#ifndef MY_MD_CODE_ATOMS_H
#define MY_MD_CODE_ATOMS_H

#include "types.h"

class Atoms {
  public:
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Masses_t masses;

    Atoms(const Atoms& atoms){
        positions = atoms.positions;
        velocities = atoms.velocities;
        forces = atoms.forces;
        masses = atoms.masses;
    }

    Atoms(size_t i) : positions{3, i}, velocities{3, i}, forces{3, i}, masses{3, i} {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        masses.setZero();
    }

    Atoms(const Positions_t &p) :
          positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{3, p.cols()} {
        velocities.setZero();
        forces.setZero();
        masses.setZero();
    }

    Atoms(const Positions_t &p, const Velocities_t &v) :
          positions{p}, velocities{v}, forces{3, p.cols()}, masses{3, p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
        masses.setZero();
    }

    size_t nb_atoms() const {
        return positions.cols();
    }

    void resize(int size){
        positions.conservativeResize(Eigen::NoChange,size);
        velocities.conservativeResize(Eigen::NoChange,size);
        forces.conservativeResize(Eigen::NoChange,size);
        masses.conservativeResize(Eigen::NoChange,size);
    }
    
    void setRandom(){
        positions.setRandom();
        velocities.setRandom();
        forces.setRandom();
        masses.setRandom();
    }

    /*
     * Returns an atoms object of just the first num atoms
     */
    Atoms just_the_first(int num){
        Atoms cp{*this};
        cp.resize(num);
        return cp;
    }
};


#endif // MY_MD_CODE_ATOMS_H
