#include <gtest/gtest.h>
#include "atoms.h"

TEST(AtomsResizeTest, ResizeTest){
    Atoms atoms{10};
    atoms.setRandom();

    // Copy:
    Atoms atoms_cp{atoms};

    int new_size = 5;
    atoms.resize(new_size);
    ASSERT_EQ(atoms.nb_atoms(), new_size);

    for(int i = 0; i < new_size; i++){
        for(int j = 0; j < 3; j++){
            ASSERT_EQ(atoms.positions(j, i),atoms_cp.positions(j, i));
        }
    }
}