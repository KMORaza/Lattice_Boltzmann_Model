#include "lbm.hpp"
#include <iostream>
#include <cmath>

void LBM::visualize() {
    std::cout << "\nVelocity magnitude (z=0 slice):\n";
    for (int x = 0; x < data.nx; x++) {
        for (int y = 0; y < data.ny; y++) {
            int idx = (x * data.ny + y) * data.nz;
            double mag = std::sqrt(data.ux[idx] * data.ux[idx] + data.uy[idx] * data.uy[idx] + data.uz[idx] * data.uz[idx]);
            std::cout << mag << " ";
        }
        std::cout << "\n";
    }
}
