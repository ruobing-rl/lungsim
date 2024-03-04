#include "surfactant.h"

void update_surface_tension_c(double *surf_concentration, double *surface_tension);

void update_surface_tension(double surf_concentration, double surface_tension)
{
  update_surface_tension_c(&surf_concentration, &surface_tension);
}