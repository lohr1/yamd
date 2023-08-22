#include "lj.h"

int main(int argc, char *argv[]) {
    Atoms atoms{10};
    atoms.positions.setRandom();

    double cutoff = 5.0;
    double epsilon = 1.0;
    double sigma = 1.0;
    lj_cutoff_summation(atoms, cutoff, epsilon, sigma);
}