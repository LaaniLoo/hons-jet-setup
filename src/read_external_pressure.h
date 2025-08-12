#ifndef read_external_pressure_h
#define read_external_pressure_h

#include <stdio.h>
#include <stdlib.h>
#include "pluto.h"
#include "definitions_usr.h"
#include "coordinates.h"
#include "units.h"
#include "environment_utilities.h"

void ReadHotAndGravTable(TabulatedEnvironmentProperties*);

double GravitationalAcceleration(const double x1, const double x2, const double x3, TabulatedEnvironmentProperties*);

void HotHaloPrimitives(double *halo,
                       const double x1, const double x2, const double x3, TabulatedEnvironmentProperties*);

int hunter2(const double *arr, const int narr, const double val);

int hunter(const double *arr, const int narr, const double val);

double LinearInterpolate(double y1, double y2, double fc);

int sgn(double val);
 
double InterpolationWrapper(const double arg_arr[], const double val_arr[], const int narr, const double arg);

#endif
