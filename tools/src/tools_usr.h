/** 
 * @file tools_usr.h
 *
 * Some useful tools.
*/
#ifndef tools_usr_h
#define tools_usr_h

#include "pluto.h"

double get_smoothing(double R, double Cr, int n);

double lorentz2speed(const double lorentz);

double speed2lorentz(const double speed);

double clamp(double d, double min, double max);

double between(double, double, double);

/** Calculates effective adiabatic index using the TAUB equation of state.
 * This equation of state is based on the Synge EOS (Synge 1957) for a relativistic
 * perfect gas. The actual implementation is used is as given in Mignone & McKinnney 2007,
 * and is the same as is used in PLUTO.
 *
 * @param rho density
 * @param prs pressure
 * @return The effective adiabatic index
 */
double taub_eff_adiabat_index(double rho, double prs);
#endif
