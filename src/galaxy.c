#include "galaxy.h" 

double get_stellar_density_hernquist(double stellar_mass, double x1, double x2, double x3) {
    double a_scale, r_eff, dim_r, r;

    r = sph_radius(x1, x2, x3);
    if (r < 0.001)
        r = 0.001;
    r *= vn.length; /* This is now in cm */

    r_eff = get_r_eff_elliptical(stellar_mass); /* This is in kpc */
    a_scale=(r_eff * vn.length)/1.815; /* Convert r_eff to cm, use this to calculate the scale radius */
    dim_r = r / a_scale;

    return (stellar_mass * 2e33) / (2 * CONST_PI * pow(a_scale, 3)) * (1/(dim_r*pow((1+dim_r), 3)));
}

double get_mass_injection_rate(double stellar_density) {
    return stellar_density * get_alpha(GALAXY_AGE);
}

double get_T0() {
    double sn_rate, s;
    s = 1.3;

    sn_rate = SUPERNOVAE_FACTOR * pow((GALAXY_AGE)/13.0, -s);
    return sn_rate * T0_FACTOR / get_alpha(GALAXY_AGE);
}

double get_energy_injection_rate(double stellar_density) {
    double T0 = get_T0();
    return stellar_density * get_alpha(GALAXY_AGE) * T0 * get_cV();
}

double get_cV() {
    return CONST_kB / ((nonrel_gamma - 1) * MU_NORM * CONST_mp);
}

double get_r_eff_elliptical(double stellar_mass) {
    // Stellar mass is in units of solar mass, and we return the effective radius in kpc
    double MtoL_K=0.6; // mass-to-light ratio in the Ks band (McGaugh & Schombert 2014, see their Conclusion)
    double LK=stellar_mass / MtoL_K; // Ks-band luminosity (for MtoL<1 have lumore Lsun than Mstars)
    double MK=4.74-2.5*log10(LK); // approximate K-band magnitude, from stellar mass; 4.74 is the bol mag of the Sun
    double Mb=MK+4.0; // B-band optical magnitude, using the crude conversion of Graham & Worley (2008)
    double n_sersic=pow(10, (-(14.3+Mb)/9.4)); // Graham & Worley eqn 17
    double bn=1.9992*n_sersic-0.3271; // Graham & Worley between eqns 16-17
    double logReff=1.137+0.217*bn+0.1*Mb+0.5*log10(pow(bn,2*n_sersic)/(n_sersic*exp(bn)*tgamma(2*n_sersic))); // in kpc, Eqn 16 of Graham & Worley (2008); the gamma function is *complete* as shown in Graham & Colless (1997) eqn 12
    return pow(10, logReff);
}

double get_alpha(double time) {
    // time in gigayears, return rate in seconds^-1
    return 4.7e-20 * pow(time / 13.0, -1.3);
}
