#ifndef environment_h
#define environment_h

#include "pluto.h"
#include "pluto_usr.h"
#include "environment_utilities.h"
#include "read_external_pressure.h"
#include "galaxy.h"
#include "units.h"
#include "coordinates.h"

double halo_primitives[NVAR];

void setup_analytic_environment(AnalyticEnvironmentProperties*);
void setup_external_environment(ExternalEnvironmentProperties*);
void PrintEnvironmentProperties(void*);

void set_default_environment(void*);

double density_at_position(const double x1, const double x2, const double x3, void*);
double temperature_at_position(const double x1, const double x2, const double x3, void*);
double pressure_at_position(const double x1, const double x2, const double x3, void*);
void grav_accel_at_position(const double x1, const double x2, const double x3, double *g, void*);
double potential_at_position(const double x1, const double x2, const double x3, void*);
void velocity_at_position(const double x1, const double x2, const double x3, double *v, void*);

double get_pressure(double rho, double T);
double get_temperature(double rho, double P);
double get_density(double P, double T);

double get_sound_speed(double T);

double nfw_gravitational_acceleration(const double x1, const double x2, const double x3);
double hernquist_gravitational_acceleration(const double stellar_mass, const double x1, const double x2, const double x3);

void set_external_env_paths(ExternalEnvironmentProperties*, const char*, const char*);
void set_tabulated_env_path(TabulatedEnvironmentProperties*, const char*);

#endif
