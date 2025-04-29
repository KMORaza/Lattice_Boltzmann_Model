#include "lbm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void compute_macro(LBMData *data) {
    int size = data->nx * data->ny * data->nz;
    for (int i = 0; i < size; i++) {
        double rho = 0.0, ux = 0.0, uy = 0.0, uz = 0.0;
        for (int q = 0; q < data->lattice.q; q++) {
            double fi = data->f[i * data->lattice.q + q];
            rho += fi;
            ux += fi * data->lattice.c[q][0];
            uy += fi * data->lattice.c[q][1];
            uz += fi * data->lattice.c[q][2];
        }
        data->rho[i] = rho;
        // Prevent division by zero
        if (rho < 1e-10) {
            printf("Zero density at index %d, clamping rho to 1e-10\n", i);
            rho = 1e-10;
        }
        data->ux[i] = ux / rho;
        data->uy[i] = uy / rho;
        data->uz[i] = uz / rho;
        if (isinf(data->ux[i]) || isinf(data->uy[i]) || isinf(data->uz[i]) || isnan(data->ux[i]) || isnan(data->uy[i]) || isnan(data->uz[i])) {
            printf("Numerical error at index %d: rho=%f, ux=%f, uy=%f, uz=%f\n", i, rho, data->ux[i], data->uy[i], data->uz[i]);
            exit(1);
        }
    }
}
void collide(LBMData *data, CollisionType ct) {
    int size = data->nx * data->ny * data->nz;
    int q = data->lattice.q;
    if (data->lattice.q != 9 || ct == SRT) {
        // SRT Collision
        for (int i = 0; i < size; i++) {
            double rho = data->rho[i], ux = data->ux[i], uy = data->uy[i], uz = data->uz[i];
            for (int j = 0; j < q; j++) {
                double ciu = data->lattice.c[j][0] * ux + data->lattice.c[j][1] * uy + data->lattice.c[j][2] * uz;
                double u2 = ux * ux + uy * uy + uz * uz;
                double feq = data->lattice.w[j] * rho * (1.0 + 3.0 * ciu + 4.5 * ciu * ciu - 1.5 * u2);
                data->f[i * q + j] += OMEGA * (feq - data->f[i * q + j]);
            }
        }
    } else if (ct == MRT && data->lattice.q == 9) {
        // MRT for D2Q9
        double s[9] = {0.0, 1.19, 1.4, 0.0, 1.2, 0.0, 1.2, 1.98, 1.98}; // Relaxation rates
        for (int i = 0; i < size; i++) {
            double rho = data->rho[i], ux = data->ux[i], uy = data->uy[i];
            double m[9], meq[9], f[9];
            // Transform to moment space
            for (int j = 0; j < 9; j++) f[j] = data->f[i * q + j];
            m[0] = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8]; // Density
            m[1] = -4*f[0] - f[1] - f[2] - f[3] - f[4] + 2*f[5] + 2*f[6] + 2*f[7] + 2*f[8]; // Energy
            m[2] = 2*f[0] + f[1] + f[2] + f[3] + f[4] - f[5] - f[6] - f[7] - f[8]; // Energy squared
            m[3] = f[1] - f[3]; // x-momentum
            m[4] = -f[1] + f[3]; // x-flux
            m[5] = f[2] - f[4]; // y-momentum
            m[6] = -f[2] + f[4]; // y-flux
            m[7] = f[5] + f[6] - f[7] - f[8]; // xx-yy stress
            m[8] = f[5] - f[6] + f[7] - f[8]; // xy stress
            // Equilibrium moments
            double u2 = ux * ux + uy * uy;
            meq[0] = rho;
            meq[1] = -2*rho + 3*rho*u2;
            meq[2] = rho - 3*rho*u2;
            meq[3] = rho*ux;
            meq[4] = -rho*ux;
            meq[5] = rho*uy;
            meq[6] = -rho*uy;
            meq[7] = rho*(ux*ux - uy*uy);
            meq[8] = rho*ux*uy;
            // Relaxation
            for (int j = 0; j < 9; j++) m[j] -= s[j] * (m[j] - meq[j]);
            // Transform back to distribution space
            data->f[i*q + 0] = (m[0] - 4*m[1] + 2*m[2])/9.0;
            data->f[i*q + 1] = (m[0] - m[1] + m[2] + 3*m[3] - 3*m[4])/9.0;
            data->f[i*q + 2] = (m[0] - m[1] + m[2] + 3*m[5] - 3*m[6])/9.0;
            data->f[i*q + 3] = (m[0] - m[1] + m[2] - 3*m[3] + 3*m[4])/9.0;
            data->f[i*q + 4] = (m[0] - m[1] + m[2] - 3*m[5] + 3*m[6])/9.0;
            data->f[i*q + 5] = (m[0] + 2*m[1] - m[2] + 3*m[7] + 3*m[8])/36.0;
            data->f[i*q + 6] = (m[0] + 2*m[1] - m[2] + 3*m[7] - 3*m[8])/36.0;
            data->f[i*q + 7] = (m[0] + 2*m[1] - m[2] - 3*m[7] + 3*m[8])/36.0;
            data->f[i*q + 8] = (m[0] + 2*m[1] - m[2] - 3*m[7] - 3*m[8])/36.0;
        }
    } else if (ct == TRT && data->lattice.q == 9) {
        // TRT for D2Q9
        double omega_plus = OMEGA, omega_minus = 1.0; // Symmetric and antisymmetric rates
        for (int i = 0; i < size; i++) {
            double rho = data->rho[i], ux = data->ux[i], uy = data->uy[i];
            for (int j = 0; j < q; j++) {
                double ciu = data->lattice.c[j][0] * ux + data->lattice.c[j][1] * uy;
                double u2 = ux * ux + uy * uy;
                double feq = data->lattice.w[j] * rho * (1.0 + 3.0 * ciu + 4.5 * ciu * ciu - 1.5 * u2);
                int j_opp = data->lattice.opp[j];
                double f_plus = 0.5 * (data->f[i*q + j] + data->f[i*q + j_opp]);
                double f_minus = 0.5 * (data->f[i*q + j] - data->f[i*q + j_opp]);
                double feq_plus = 0.5 * (feq + data->lattice.w[j_opp] * rho * (1.0 + 3.0 * (data->lattice.c[j_opp][0] * ux + data->lattice.c[j_opp][1] * uy) + 4.5 * (data->lattice.c[j_opp][0] * ux + data->lattice.c[j_opp][1] * uy) * (data->lattice.c[j_opp][0] * ux + data->lattice.c[j_opp][1] * uy) - 1.5 * u2));
                double feq_minus = 0.5 * (feq - data->lattice.w[j_opp] * rho * (1.0 + 3.0 * (data->lattice.c[j_opp][0] * ux + data->lattice.c[j_opp][1] * uy) + 4.5 * (data->lattice.c[j_opp][0] * ux + data->lattice.c[j_opp][1] * uy) * (data->lattice.c[j_opp][0] * ux + data->lattice.c[j_opp][1] * uy) - 1.5 * u2));
                data->f[i*q + j] = f_plus - omega_plus * (f_plus - feq_plus) + f_minus - omega_minus * (f_minus - feq_minus);
            }
        }
    } else if (ct == ENTROPIC && data->lattice.q == 9) {
        // Entropic for D2Q9
        for (int i = 0; i < size; i++) {
            double rho = data->rho[i], ux = data->ux[i], uy = data->uy[i];
            double H = 0.0, H_eq = 0.0, alpha = 1.0;
            double f[9], feq[9];
            for (int j = 0; j < 9; j++) {
                f[j] = data->f[i*q + j];
                double ciu = data->lattice.c[j][0] * ux + data->lattice.c[j][1] * uy;
                double u2 = ux * ux + uy * uy;
                feq[j] = data->lattice.w[j] * rho * (1.0 + 3.0 * ciu + 4.5 * ciu * ciu - 1.5 * u2);
                H += f[j] * log(f[j] / data->lattice.w[j]);
                H_eq += feq[j] * log(feq[j] / data->lattice.w[j]);
            }
            // Find alpha to minimize entropy difference
            double alpha_min = 0.0, alpha_max = 2.0;
            for (int iter = 0; iter < 20; iter++) {
                alpha = 0.5 * (alpha_min + alpha_max);
                double H_alpha = 0.0;
                for (int j = 0; j < 9; j++) {
                    double f_alpha = f[j] + alpha * (feq[j] - f[j]);
                    H_alpha += f_alpha * log(f_alpha / data->lattice.w[j]);
                }
                if (H_alpha > H_eq) alpha_max = alpha;
                else alpha_min = alpha;
            }
            // Apply relaxation
            for (int j = 0; j < 9; j++) {
                data->f[i*q + j] += alpha * (feq[j] - f[j]);
            }
        }
    }
}
