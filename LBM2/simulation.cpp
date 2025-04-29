#include "lbm.hpp"

LBM::LBM(LatticeType lt, CollisionType ct, BoundaryType bt, int steps)
    : lattice_type(lt), collision_type(ct), boundary_type(bt), max_steps(steps) {
    init_lbm_data(&data, lt);
}

LBM::~LBM() {
    free_lbm_data(&data);
}

void LBM::run() {
    for (int step = 0; step < max_steps; step++) {
        compute_macro(&data);
        collide(&data, collision_type);
        stream(&data);
        apply_boundary(&data, boundary_type);
        if (step % 10 == 0) visualize();
    }
}
