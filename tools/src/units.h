#ifndef units_h
#define units_h

#include "pluto.h"

/* This structure holds the normalisation variables for our simulation.
 * They convert from code to cgs units, i.e.
 * <quantity_in_code_units> * <relevant_normalisation> = <quantity_in_cgs_units>
 */
typedef struct {
    double length;
    double density;
    double velocity;
    double temperature;
    double time;
    double area;
    double pressure;
    double power;
    double energy_flux;
    double energy;
    double gravitational_potential;
    double acceleration;
    double number_density;
    double mass;
} VariableNormalisation;

extern VariableNormalisation vn;

void SetNormalisation();
void PrintNormalisations();

#endif
