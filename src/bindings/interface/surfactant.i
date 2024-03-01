%module(package="aether") surfactant
%include symbol_export.h
%include surfactant.h

%{
#include "surfactant.h"
%}

void calculate_surface_tension(real volume_mean, real volume_change, real frequency, real volumes, &
	real radii, real area, real dA, real surf_concentration, real gamma_star, real gamma_max, real bulk_c, &
	real k_a, real k_d, real m2, real sigma, real sigma_hat);