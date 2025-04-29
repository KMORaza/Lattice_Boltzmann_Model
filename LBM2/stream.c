#include "lbm.h"
#include <stdio.h>
#include <stdlib.h>

void stream(LBMData *data) {
    int nx = data->nx, ny = data->ny, nz = data->nz;
    int q = data->lattice.q;
    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {
            for (int z = 0; z < nz; z++) {
                int i = (x * ny + y) * nz + z;
                for (int j = 0; j < q; j++) {
                    // Validate lattice velocity
                    if (data->lattice.c[j][2] != 0 && data->lattice.q == 9) {
                        printf("Invalid z-velocity for D2Q9 at dir=%d: c=(%f, %f, %f)\n",
                               j, data->lattice.c[j][0], data->lattice.c[j][1], data->lattice.c[j][2]);
                        exit(1);
                    }
                    int x_new = x + (int)data->lattice.c[j][0];
                    int y_new = y + (int)data->lattice.c[j][1];
                    int z_new = z + (int)data->lattice.c[j][2]; 
                    if (x_new >= 0 && x_new < nx && y_new >= 0 && y_new < ny && z_new >= 0 && z_new < nz) {
                        int i_new = (x_new * ny + y_new) * nz + z_new;
                        data->f_new[i_new * q + j] = data->f[i * q + j];
                    } else {
                        printf("Out-of-bounds streaming at x=%d, y=%d, z=%d, dir=%d: x_new=%d, y_new=%d, z_new=%d\n",
                               x, y, z, j, x_new, y_new, z_new);
                    }
                }
            }
        }
    }
    double *temp = data->f;
    data->f = data->f_new;
    data->f_new = temp;
    int debug_index = 0;
    printf("After streaming, f at index %d: ", debug_index);
    for (int j = 0; j < q; j++) {
        printf("%f ", data->f[debug_index * q + j]);
    }
    printf("\n");
}
