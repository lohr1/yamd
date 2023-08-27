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

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    // Code goes here - MPI comm can be used
    MPI_Finalize();
}
