#ifndef LBM_HPP
#define LBM_HPP

#include "lbm.h"
#include <vector>

class LBM {
private:
    LatticeType lattice_type;
    CollisionType collision_type;
    BoundaryType boundary_type;
    LBMData data;
    int max_steps;

public:
    LBM(LatticeType lt, CollisionType ct, BoundaryType bt, int steps = 100);
    ~LBM();
    void run();
    void visualize();
};

#endif
