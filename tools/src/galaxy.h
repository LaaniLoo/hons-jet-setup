#ifndef galaxy_h
#define galaxy_h

#include "pluto.h"
#include "pluto_usr.h"
/* #include "read_external_pressure.h" */
#include "units.h"
#include "coordinates.h"

double get_stellar_density_hernquist(double stellar_mass, double x1, double x2, double x3);
double get_mass_injection_rate(double stellar_density);
double get_energy_injection_rate(double stellar_density);
double get_T0();
double get_cV();
double get_r_eff_elliptical(double stellar_mass);
double get_alpha(double time);

#endif
