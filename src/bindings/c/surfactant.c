#include "surfactant.h"

void calculate_surface_tension_c(real *volume_mean, real *volume_change, real *frequency, real *volumes, &
	real *radii, real *area, real *dA, real *surf_concentration, real *gamma_star, real *gamma_max, real *bulk_c, &
	real *k_a, real *k_d, real *m2, real *sigma, real *sigma_hat);

void calculate_surface_tension(real volume_mean, real volume_change, real frequency, real volumes, &
	real radii, real area, real dA, real surf_concentration, real gamma_star, real gamma_max, real bulk_c, &
	real k_a, real k_d, real m2, real sigma, real sigma_hat)
{
  calculate_surface_tension_c(&volume_mean, &volume_change, &frequency, &volumes, &radii, &area, &dA, &
    &surf_concentration, &gamma_star, &gamma_max, &bulk_c, &k_a, &k_d, &m2, &sigma, &sigma_hat);
}