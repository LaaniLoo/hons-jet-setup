#include "environment_utilities.h"

inline double calculate_offset_radius(const double x1, const double x2, const double x3, void* bep) {

    /* -- Let's calculate the offset radius  -- */
#if DENSITY_PROFILE == ENV_MAKINO || DENSITY_PROFILE == ENV_KING || DENSITY_PROFILE == ENV_CONST
    AnalyticEnvironmentProperties* ep = (AnalyticEnvironmentProperties*) bep;
    double *offsets = ep->offset;
#elif DENSITY_PROFILE == ENV_TABULATED
    TabulatedEnvironmentProperties* tep = (TabulatedEnvironmentProperties*) bep;
    double *offsets = tep->offset;
#elif DENSITY_PROFILE == ENV_EXTERNAL
    double *offsets = ((ExternalEnvironmentProperties*) bep)->offset;
#endif

#if COMPONENTS == 2 && GEOMETRY == SPHERICAL
    double r_offset = offsets[0];
    double theta = sph_theta(x1, x2, x3);
    double r = sph_radius(x1, x2, x3);
    double r_prime = sqrt(r*r + r_offset*r_offset - 2*r*r_offset*cos(theta));
    return r_prime;
#elif GEOMETRY == CARTESIAN
    return sph_radius(x1 - offsets[0], x2 - offsets[1], x3 - offsets[2]);
#endif
}

