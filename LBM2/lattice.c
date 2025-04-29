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
    } else if (type == D3Q15) {
        lattice->q = 15;
        // Velocities: {x, y, z}
        lattice->c[0][0] = 0; lattice->c[0][1] = 0; lattice->c[0][2] = 0; // Rest
        lattice->c[1][0] = 1; lattice->c[1][1] = 0; lattice->c[1][2] = 0; // East
        lattice->c[2][0] = -1; lattice->c[2][1] = 0; lattice->c[2][2] = 0; // West
        lattice->c[3][0] = 0; lattice->c[3][1] = 1; lattice->c[3][2] = 0; // North
        lattice->c[4][0] = 0; lattice->c[4][1] = -1; lattice->c[4][2] = 0; // South
        lattice->c[5][0] = 0; lattice->c[5][1] = 0; lattice->c[5][2] = 1; // Up
        lattice->c[6][0] = 0; lattice->c[6][1] = 0; lattice->c[6][2] = -1; // Down
        lattice->c[7][0] = 1; lattice->c[7][1] = 1; lattice->c[7][2] = 1; // Northeast-Up
        lattice->c[8][0] = -1; lattice->c[8][1] = -1; lattice->c[8][2] = -1; // Southwest-Down
        lattice->c[9][0] = 1; lattice->c[9][1] = -1; lattice->c[9][2] = 1; // Southeast-Up
        lattice->c[10][0] = -1; lattice->c[10][1] = 1; lattice->c[10][2] = -1; // Northwest-Down
        lattice->c[11][0] = 1; lattice->c[11][1] = 1; lattice->c[11][2] = -1; // Northeast-Down
        lattice->c[12][0] = -1; lattice->c[12][1] = -1; lattice->c[12][2] = 1; // Southwest-Up
        lattice->c[13][0] = 1; lattice->c[13][1] = -1; lattice->c[13][2] = -1; // Southeast-Down
        lattice->c[14][0] = -1; lattice->c[14][1] = 1; lattice->c[14][2] = 1; // Northwest-Up
        // Weights
        lattice->w[0] = 2.0/9.0;
        lattice->w[1] = lattice->w[2] = lattice->w[3] = lattice->w[4] = lattice->w[5] = lattice->w[6] = 1.0/9.0;
        lattice->w[7] = lattice->w[8] = lattice->w[9] = lattice->w[10] = lattice->w[11] = lattice->w[12] = lattice->w[13] = lattice->w[14] = 1.0/72.0;
        // Opposite directions
        lattice->opp[0] = 0; // Rest
        lattice->opp[1] = 2; // East <-> West
        lattice->opp[2] = 1; // West <-> East
        lattice->opp[3] = 4; // North <-> South
        lattice->opp[4] = 3; // South <-> North
        lattice->opp[5] = 6; // Up <-> Down
        lattice->opp[6] = 5; // Down <-> Up
        lattice->opp[7] = 8; // Northeast-Up <-> Southwest-Down
        lattice->opp[8] = 7; // Southwest-Down <-> Northeast-Up
        lattice->opp[9] = 10; // Southeast-Up <-> Northwest-Down
        lattice->opp[10] = 9; // Northwest-Down <-> Southeast-Up
        lattice->opp[11] = 12; // Northeast-Down <-> Southwest-Up
        lattice->opp[12] = 11; // Southwest-Up <-> Northeast-Down
        lattice->opp[13] = 14; // Southeast-Down <-> Northwest-Up
        lattice->opp[14] = 13; // Northwest-Up <-> Southeast-Down
    } else if (type == D3Q19) {
        lattice->q = 19;
        // Velocities: {x, y, z}
        lattice->c[0][0] = 0; lattice->c[0][1] = 0; lattice->c[0][2] = 0; // Rest
        lattice->c[1][0] = 1; lattice->c[1][1] = 0; lattice->c[1][2] = 0; // East
        lattice->c[2][0] = -1; lattice->c[2][1] = 0; lattice->c[2][2] = 0; // West
        lattice->c[3][0] = 0; lattice->c[3][1] = 1; lattice->c[3][2] = 0; // North
        lattice->c[4][0] = 0; lattice->c[4][1] = -1; lattice->c[4][2] = 0; // South
        lattice->c[5][0] = 0; lattice->c[5][1] = 0; lattice->c[5][2] = 1; // Up
        lattice->c[6][0] = 0; lattice->c[6][1] = 0; lattice->c[6][2] = -1; // Down
        lattice->c[7][0] = 1; lattice->c[7][1] = 1; lattice->c[7][2] = 0; // Northeast
        lattice->c[8][0] = -1; lattice->c[8][1] = -1; lattice->c[8][2] = 0; // Southwest
        lattice->c[9][0] = 1; lattice->c[9][1] = -1; lattice->c[9][2] = 0; // Southeast
        lattice->c[10][0] = -1; lattice->c[10][1] = 1; lattice->c[10][2] = 0; // Northwest
        lattice->c[11][0] = 1; lattice->c[11][1] = 0; lattice->c[11][2] = 1; // East-Up
        lattice->c[12][0] = -1; lattice->c[12][1] = 0; lattice->c[12][2] = -1; // West-Down
        lattice->c[13][0] = 1; lattice->c[13][1] = 0; lattice->c[13][2] = -1; // East-Down
        lattice->c[14][0] = -1; lattice->c[14][1] = 0; lattice->c[14][2] = 1; // West-Up
        lattice->c[15][0] = 0; lattice->c[15][1] = 1; lattice->c[15][2] = 1; // North-Up
        lattice->c[16][0] = 0; lattice->c[16][1] = -1; lattice->c[16][2] = -1; // South-Down
        lattice->c[17][0] = 0; lattice->c[17][1] = 1; lattice->c[17][2] = -1; // North-Down
        lattice->c[18][0] = 0; lattice->c[18][1] = -1; lattice->c[18][2] = 1; // South-Up
        // Weights
        lattice->w[0] = 1.0/3.0;
        lattice->w[1] = lattice->w[2] = lattice->w[3] = lattice->w[4] = lattice->w[5] = lattice->w[6] = 1.0/18.0;
        lattice->w[7] = lattice->w[8] = lattice->w[9] = lattice->w[10] = lattice->w[11] = lattice->w[12] = 
        lattice->w[13] = lattice->w[14] = lattice->w[15] = lattice->w[16] = lattice->w[17] = lattice->w[18] = 1.0/36.0;
        // Opposite directions
        lattice->opp[0] = 0; // Rest
        lattice->opp[1] = 2; // East <-> West
        lattice->opp[2] = 1; // West <-> East
        lattice->opp[3] = 4; // North <-> South
        lattice->opp[4] = 3; // South <-> North
        lattice->opp[5] = 6; // Up <-> Down
        lattice->opp[6] = 5; // Down <-> Up
        lattice->opp[7] = 8; // Northeast <-> Southwest
        lattice->opp[8] = 7; // Southwest <-> Northeast
        lattice->opp[9] = 10; // Southeast <-> Northwest
        lattice->opp[10] = 9; // Northwest <-> Southeast
        lattice->opp[11] = 12; // East-Up <-> West-Down
        lattice->opp[12] = 11; // West-Down <-> East-Up
        lattice->opp[13] = 14; // East-Down <-> West-Up
        lattice->opp[14] = 13; // West-Up <-> East-Down
        lattice->opp[15] = 16; // North-Up <-> South-Down
        lattice->opp[16] = 15; // South-Down <-> North-Up
        lattice->opp[17] = 18; // North-Down <-> South-Up
        lattice->opp[18] = 17; // South-Up <-> North-Down
    } else if (type == D3Q27) {
        lattice->q = 27;
        // Velocities: {x, y, z}
        lattice->c[0][0] = 0; lattice->c[0][1] = 0; lattice->c[0][2] = 0; // Rest
        lattice->c[1][0] = 1; lattice->c[1][1] = 0; lattice->c[1][2] = 0; // East
        lattice->c[2][0] = -1; lattice->c[2][1] = 0; lattice->c[2][2] = 0; // West
        lattice->c[3][0] = 0; lattice->c[3][1] = 1; lattice->c[3][2] = 0; // North
        lattice->c[4][0] = 0; lattice->c[4][1] = -1; lattice->c[4][2] = 0; // South
        lattice->c[5][0] = 0; lattice->c[5][1] = 0; lattice->c[5][2] = 1; // Up
        lattice->c[6][0] = 0; lattice->c[6][1] = 0; lattice->c[6][2] = -1; // Down
        lattice->c[7][0] = 1; lattice->c[7][1] = 1; lattice->c[7][2] = 0; // Northeast
        lattice->c[8][0] = -1; lattice->c[8][1] = -1; lattice->c[8][2] = 0; // Southwest
        lattice->c[9][0] = 1; lattice->c[9][1] = -1; lattice->c[9][2] = 0; // Southeast
        lattice->c[10][0] = -1; lattice->c[10][1] = 1; lattice->c[10][2] = 0; // Northwest
        lattice->c[11][0] = 1; lattice->c[11][1] = 0; lattice->c[11][2] = 1; // East-Up
        lattice->c[12][0] = -1; lattice->c[12][1] = 0; lattice->c[12][2] = -1; // West-Down
        lattice->c[13][0] = 1; lattice->c[13][1] = 0; lattice->c[13][2] = -1; // East-Down
        lattice->c[14][0] = -1; lattice->c[14][1] = 0; lattice->c[14][2] = 1; // West-Up
        lattice->c[15][0] = 0; lattice->c[15][1] = 1; lattice->c[15][2] = 1; // North-Up
        lattice->c[16][0] = 0; lattice->c[16][1] = -1; lattice->c[16][2] = -1; // South-Down
        lattice->c[17][0] = 0; lattice->c[17][1] = 1; lattice->c[17][2] = -1; // North-Down
        lattice->c[18][0] = 0; lattice->c[18][1] = -1; lattice->c[18][2] = 1; // South-Up
        lattice->c[19][0] = 1; lattice->c[19][1] = 1; lattice->c[19][2] = 1; // Northeast-Up
        lattice->c[20][0] = -1; lattice->c[20][1] = -1; lattice->c[20][2] = -1; // Southwest-Down
        lattice->c[21][0] = 1; lattice->c[21][1] = -1; lattice->c[21][2] = 1; // Southeast-Up
        lattice->c[22][0] = -1; lattice->c[22][1] = 1; lattice->c[22][2] = -1; // Northwest-Down
        lattice->c[23][0] = 1; lattice->c[23][1] = 1; lattice->c[23][2] = -1; // Northeast-Down
        lattice->c[24][0] = -1; lattice->c[24][1] = -1; lattice->c[24][2] = 1; // Southwest-Up
        lattice->c[25][0] = 1; lattice->c[25][1] = -1; lattice->c[25][2] = -1; // Southeast-Down
        lattice->c[26][0] = -1; lattice->c[26][1] = 1; lattice->c[26][2] = 1; // Northwest-Up
        // Weights
        lattice->w[0] = 8.0/27.0;
        lattice->w[1] = lattice->w[2] = lattice->w[3] = lattice->w[4] = lattice->w[5] = lattice->w[6] = 2.0/27.0;
        lattice->w[7] = lattice->w[8] = lattice->w[9] = lattice->w[10] = lattice->w[11] = lattice->w[12] = 
        lattice->w[13] = lattice->w[14] = lattice->w[15] = lattice->w[16] = lattice->w[17] = lattice->w[18] = 1.0/54.0;
        lattice->w[19] = lattice->w[20] = lattice->w[21] = lattice->w[22] = lattice->w[23] = lattice->w[24] = 
        lattice->w[25] = lattice->w[26] = 1.0/216.0;
        // Opposite directions
        lattice->opp[0] = 0; // Rest
        lattice->opp[1] = 2; // East <-> West
        lattice->opp[2] = 1; // West <-> East
        lattice->opp[3] = 4; // North <-> South
        lattice->opp[4] = 3; // South <-> North
        lattice->opp[5] = 6; // Up <-> Down
        lattice->opp[6] = 5; // Down <-> Up
        lattice->opp[7] = 8; // Northeast <-> Southwest
        lattice->opp[8] = 7; // Southwest <-> Northeast
        lattice->opp[9] = 10; // Southeast <-> Northwest
        lattice->opp[10] = 9; // Northwest <-> Southeast
        lattice->opp[11] = 12; // East-Up <-> West-Down
        lattice->opp[12] = 11; // West-Down <-> East-Up
        lattice->opp[13] = 14; // East-Down <-> West-Up
        lattice->opp[14] = 13; // West-Up <-> East-Down
        lattice->opp[15] = 16; // North-Up <-> South-Down
        lattice->opp[16] = 15; // South-Down <-> North-Up
        lattice->opp[17] = 18; // North-Down <-> South-Up
        lattice->opp[18] = 17; // South-Up <-> North-Down
        lattice->opp[19] = 20; // Northeast-Up <-> Southwest-Down
        lattice->opp[20] = 19; // Southwest-Down <-> Northeast-Up
        lattice->opp[21] = 22; // Southeast-Up <-> Northwest-Down
        lattice->opp[22] = 21; // Northwest-Down <-> Southeast-Up
        lattice->opp[23] = 24; // Northeast-Down <-> Southwest-Up
        lattice->opp[24] = 23; // Southwest-Up <-> Northeast-Down
        lattice->opp[25] = 26; // Southeast-Down <-> Northwest-Up
        lattice->opp[26] = 25; // Northwest-Up <-> Southeast-Down
    }

    // Debug output for lattice initialization
    printf("%s Lattice Velocities (q=%d):\n", 
           type == D2Q9 ? "D2Q9" : type == D3Q15 ? "D3Q15" : type == D3Q19 ? "D3Q19" : "D3Q27", 
           lattice->q);
    for (int j = 0; j < lattice->q; j++) {
        printf("Dir %d: c=(%f, %f, %f), w=%f, opp=%d\n",
               j, lattice->c[j][0], lattice->c[j][1], lattice->c[j][2], lattice->w[j], lattice->opp[j]);
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