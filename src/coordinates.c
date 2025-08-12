#include "coordinates.h"

/* Spherical coordinates */

/* Spherical radius */
inline double sph_radius(const double x1, const double x2, const double x3) {
#if GEOMETRY == CARTESIAN
    return sqrt(EXPAND(x1*x1, + x2*x2, + x3*x3));
#elif GEOMETRY == POLAR
    return sqrt(x1*x1 + x3*x3);
#elif GEOMETRY == SPHERICAL
    return x1;
#elif GEOMETRY == CYLINDRICAL
    return sqrt(x1*x1 + x2*x2);
#endif
}

/* Spherical polar angle */
inline double sph_theta(const double x1, const double x2, const double x3) {
#if GEOMETRY == CARTESIAN
    return acos((D_SELECT(x3, x2, x3))/(sph_radius(x1, x2, x3)));
#elif GEOMETRY == POLAR
    return acos(x3/sph_radius(x1, x2, x3));
#elif GEOMETRY == SPHERICAL
    return x2;
#elif GEOMETRY == CYLINDRICAL
    return acos(x2/sph_radius(x1, x2, x3));
#endif
}

/* Spherical azimuth angle */
inline double sph_phi(const double x1, const double x2, const double x3) {
#if GEOMETRY == CARTESIAN
    double phi;
    phi = atan2(x2, x1);
    if (phi < 0.0)
        phi += 2.0*CONST_PI;
    return phi;
#elif GEOMETRY == POLAR
    return x3;
#elif GEOMETRY == SPHERICAL
    return x3;
#elif GEOMETRY == CYLINDRICAL
    print ("! Coordinate conversion: invalid or unsupported coordinate transformation.\n");
    QUIT_PLUTO(1);
#endif
}

/* Cylindrical radius */
inline double cyl_radius(const double x1, const double x2, const double x3) {
#if GEOMETRY == CARTESIAN
    return sqrt(x1*x1 + x2*x2);
#elif GEOMETRY == POLAR
    print ("! Coordinate conversion: invalid or unsupported coordinate transformation.\n");
    QUIT_PLUTO(1);
#elif GEOMETRY == SPHERICAL
    print ("! Coordinate conversion: invalid or unsupported coordinate transformation.\n");
    QUIT_PLUTO(1);
#elif GEOMETRY == CYLINDRICAL
    return x1;
#endif
}

/* Cylindrical Z */
inline double cyl_z(const double x1, const double x2, const double x3) {
#if GEOMETRY == CARTESIAN
    return x3;
#elif GEOMETRY == POLAR
    print ("! Coordinate conversion: invalid or unsupported coordinate transformation.\n");
    QUIT_PLUTO(1);
#elif GEOMETRY == SPHERICAL
    print ("! Coordinate conversion: invalid or unsupported coordinate transformation.\n");
    QUIT_PLUTO(1);
#elif GEOMETRY == CYLINDRICAL
    return x2;
#endif
}

inline void sph2cart(const double r, const double theta, const double phi,
                     double* x1, double* x2, double* x3) {
    *x1 = r * sin(theta) * cos(phi);
    *x2 = r * sin(theta) * sin(phi);
    *x3 = r * cos(theta);
}
