#include "surfactant.h"

//void update_surface_tension_c(double *surf_concentration, double *surface_tension);

//void update_surface_tension(double surf_concentration, double surface_tension)
//{
//  update_surface_tension_c(&surf_concentration, &surface_tension);
//}

void evaluate_surf_c ();
//(int *num_steps, double *dt, double *t, double *volumes, double *radii, double *area, double *dA, double *surf_concentration, double *surface_tension)
void evaluate_surf ()//(int num_steps, double dt, double t, double volumes, double radii, double area, double dA, double surf_concentration, double surface_tension)
{
  evaluate_surf_c ();
  //(&num_steps, &dt, &t, &volumes, &radii, &area, &dA, &surf_concentration, &surface_tension)
}