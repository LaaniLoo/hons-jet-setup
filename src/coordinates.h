#ifndef coordinates_h
#define coordinates_h

#include "pluto.h"


/* Utilities for handling different coordinate systems, and converting between them */

/* Spherical coordinates */
double sph_radius(const double, const double, const double);
double sph_theta(const double, const double, const double);
double sph_phi(const double, const double, const double);

/* Cylindrical coordinates */
double cyl_radius(const double, const double, const double);
double cyl_z(const double, const double, const double);


void sph2cart(const double, const double, const double, double*, double*, double*);

#endif
