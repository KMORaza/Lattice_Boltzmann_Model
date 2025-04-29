#ifndef LBM_H
#define LBM_H
#define NX 100
#define NY 100
#define NZ 1 // Set to 1 for D2Q9 (2D lattice)
#define OMEGA 1.0
#define VELOCITY_INLET 0.001
#define DENSITY_OUTLET 1.0

typedef enum { D2Q9, D3Q15, D3Q19, D3Q27 } LatticeType;
typedef enum { SRT, MRT, TRT, ENTROPIC } CollisionType;
typedef enum { BOUNCE_BACK, VELOCITY, PRESSURE, PERIODIC, INLET_OUTLET, OPEN } BoundaryType;

typedef struct {
    double c[27][3]; // Velocity vectors
    double w[27];    // Weights
    int q;           // Number of directions
    int opp[27];     // Opposite directions
} Lattice;
typedef struct {
    Lattice lattice;
    double *f, *f_new; // Distribution functions
    double *rho;       // Density
    double *ux, *uy, *uz; // Velocities
    int nx, ny, nz;
} LBMData;

#ifdef __cplusplus
extern "C" {
#endif

void init_lattice(Lattice *lattice, LatticeType type);
void init_lbm_data(LBMData *data, LatticeType lt);
void compute_macro(LBMData *data);
void collide(LBMData *data, CollisionType ct);
void stream(LBMData *data);
void apply_boundary(LBMData *data, BoundaryType bt);
void free_lbm_data(LBMData *data);

#ifdef __cplusplus
}
#endif

#endif
