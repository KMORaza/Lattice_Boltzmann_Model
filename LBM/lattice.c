#include "lbm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void init_lattice(Lattice *lattice, LatticeType type) {
    memset(lattice->c, 0, sizeof(lattice->c));
    memset(lattice->w, 0, sizeof(lattice->w));
    memset(lattice->opp, 0, sizeof(lattice->opp));
    if (type == D2Q9) {
        lattice->q = 9;
        // Velocities: {x, y, z}
        lattice->c[0][0] = 0; lattice->c[0][1] = 0; lattice->c[0][2] = 0; // Rest
        lattice->c[1][0] = 1; lattice->c[1][1] = 0; lattice->c[1][2] = 0; // East
        lattice->c[2][0] = 0; lattice->c[2][1] = 1; lattice->c[2][2] = 0; // North
        lattice->c[3][0] = -1; lattice->c[3][1] = 0; lattice->c[3][2] = 0; // West
        lattice->c[4][0] = 0; lattice->c[4][1] = -1; lattice->c[4][2] = 0; // South
        lattice->c[5][0] = 1; lattice->c[5][1] = 1; lattice->c[5][2] = 0; // Northeast
        lattice->c[6][0] = -1; lattice->c[6][1] = 1; lattice->c[6][2] = 0; // Northwest
        lattice->c[7][0] = -1; lattice->c[7][1] = -1; lattice->c[7][2] = 0; // Southwest
        lattice->c[8][0] = 1; lattice->c[8][1] = -1; lattice->c[8][2] = 0; // Southeast
        // Weights
        lattice->w[0] = 4.0/9.0;
        lattice->w[1] = lattice->w[2] = lattice->w[3] = lattice->w[4] = 1.0/9.0;
        lattice->w[5] = lattice->w[6] = lattice->w[7] = lattice->w[8] = 1.0/36.0;
        // Opposite directions
        lattice->opp[0] = 0; // Rest
        lattice->opp[1] = 3; // East <-> West
        lattice->opp[2] = 4; // North <-> South
        lattice->opp[3] = 1; // West <-> East
        lattice->opp[4] = 2; // South <-> North
        lattice->opp[5] = 7; // Northeast <-> Southwest
        lattice->opp[6] = 8; // Northwest <-> Southeast
        lattice->opp[7] = 5; // Southwest <-> Northeast
        lattice->opp[8] = 6; // Southeast <-> Northwest
        printf("D2Q9 Lattice Velocities:\n");
        for (int j = 0; j < lattice->q; j++) {
            printf("Dir %d: c=(%f, %f, %f), w=%f, opp=%d\n",
                   j, lattice->c[j][0], lattice->c[j][1], lattice->c[j][2], lattice->w[j], lattice->opp[j]);
        }
    } else if (type == D3Q15) {
        lattice->q = 15;
        // Initialize D3Q15 (not needed for this case)
    } else if (type == D3Q19) {
        lattice->q = 19;
        // Initialize D3Q19 (not needed)
    } else if (type == D3Q27) {
        lattice->q = 27;
        // Initialize D3Q27 (not needed)
    }
}

void init_lbm_data(LBMData *data, LatticeType lt) {
    init_lattice(&data->lattice, lt);
    data->nx = NX;
    data->ny = NY;
    data->nz = NZ;
    int size = NX * NY * NZ;
    int q = data->lattice.q;
    data->f = (double *)calloc(size * q, sizeof(double));
    data->f_new = (double *)calloc(size * q, sizeof(double));
    data->rho = (double *)calloc(size, sizeof(double));
    data->ux = (double *)calloc(size, sizeof(double));
    data->uy = (double *)calloc(size, sizeof(double));
    data->uz = (double *)calloc(size, sizeof(double));
    // Initialize equilibrium distribution
    for (int i = 0; i < size; i++) {
        data->rho[i] = 1.0;
        data->ux[i] = 0.0;
        data->uy[i] = 0.0;
        data->uz[i] = 0.0;
        for (int j = 0; j < q; j++) {
            double ciu = data->lattice.c[j][0] * data->ux[i] + 
                        data->lattice.c[j][1] * data->uy[i] + 
                        data->lattice.c[j][2] * data->uz[i];
            double u2 = data->ux[i] * data->ux[i] + data->uy[i] * data->uy[i] + data->uz[i] * data->uz[i];
            data->f[i * q + j] = data->lattice.w[j] * data->rho[i] * 
                                (1.0 + 3.0 * ciu + 4.5 * ciu * ciu - 1.5 * u2);
        }
    }
}
void free_lbm_data(LBMData *data) {
    free(data->f);
    free(data->f_new);
    free(data->rho);
    free(data->ux);
    free(data->uy);
    free(data->uz);
}
