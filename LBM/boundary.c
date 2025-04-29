#include "lbm.h"
#include <stdio.h>
#include <stdlib.h>

void apply_boundary(LBMData *data, BoundaryType bt) {
    int nx = data->nx, ny = data->ny, nz = data->nz;
    int q = data->lattice.q;
    
    if (bt == BOUNCE_BACK) {
        // Bounce-back on all boundaries
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int z = 0; z < nz; z++) {
                    if (x == 0 || x == nx-1 || y == 0 || y == ny-1) {
                        int i = (x * ny + y) * nz + z;
                        for (int j = 0; j < q; j++) {
                            int j_opp = data->lattice.opp[j];
                            data->f[i * q + j] = data->f[i * q + j_opp];
                        }
                    }
                }
            }
        }
    } else if (bt == VELOCITY) {
        // Zou-He velocity boundary at x=0 (ux=VELOCITY_INLET, uy=0), bounce-back elsewhere
        for (int y = 0; y < ny; y++) {
            for (int z = 0; z < nz; z++) {
                int i = (0 * ny + y) * nz + z; // x=0
                double rho = 0.0;
                // Compute rho from known distributions (exclude f[1,5,8])
                for (int j = 0; j < q; j++) {
                    if (j != 1 && j != 5 && j != 8) {
                        rho += data->f[i * q + j];
                    }
                }
                rho /= (1.0 - VELOCITY_INLET);
                // Set unknown distributions (f[1,5,8])
                data->f[i * q + 1] = data->f[i * q + 3] + (2.0/3.0) * rho * VELOCITY_INLET;
                data->f[i * q + 5] = data->f[i * q + 7] + 0.5 * (data->f[i * q + 4] - data->f[i * q + 2]) + (1.0/6.0) * rho * VELOCITY_INLET;
                data->f[i * q + 8] = data->f[i * q + 6] + 0.5 * (data->f[i * q + 2] - data->f[i * q + 4]) + (1.0/6.0) * rho * VELOCITY_INLET;
            }
        }
        // Bounce-back on x=nx-1, y=0, y=ny-1
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                for (int z = 0; z < nz; z++) {
                    if (x == nx-1 || y == 0 || y == ny-1) {
                        int i = (x * ny + y) * nz + z;
                        for (int j = 0; j < q; j++) {
                            int j_opp = data->lattice.opp[j];
                            data->f[i * q + j] = data->f[i * q + j_opp];
                        }
                    }
                }
            }
        }
        int debug_index = 0;
        printf("After boundary, f at index %d: ", debug_index);
        for (int j = 0; j < q; j++) {
            printf("%f ", data->f[debug_index * q + j]);
        }
        printf("\n");
    }
}
