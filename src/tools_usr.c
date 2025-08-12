#include "tools_usr.h"

double get_smoothing(double R, double Cr, int n)
{
    return 1.0/cosh(pow(R/Cr, n));
}

double lorentz2speed(const double lorentz) {
    return sqrt(1. - 1. / (lorentz * lorentz)) * (CONST_c / UNIT_VELOCITY);
}

double speed2lorentz(const double speed) {
    return 1. / sqrt(1 - pow((speed) / (CONST_c / UNIT_VELOCITY), 2));
}

double clamp(double d, double min, double max) {
    const double t = d < min ? min : d;
    return t > max ? max : t;
}

double between(double a, double b, double c) {
    return (a <= b && b <= c) || (c <= b && b <= a);
}

double taub_eff_adiabat_index(double rho, double prs) {
    double theta, h;

    // For RHD and MRHD, c / unit_velocity = 1.
    // We include it explicitly for the non-relativistic case.
    theta = prs / (rho * pow(CONST_c / UNIT_VELOCITY, 2));

    // The enthalpy calculation follows Mignone+2007,
    // and also PLUTO/Src/EOS/Taub/eos.c
    h = 2.5*theta + sqrt(2.25*theta*theta + 1.0);

    return (h - 1)/(h - 1 - theta);
}
