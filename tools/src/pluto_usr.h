#ifndef pluto_usr_h
#define pluto_usr_h

#include "definitions_usr.h"

/** Environment types. */
#define  ENV_MAKINO    0
#define  ENV_KING      1
#define  ENV_TABULATED 2
#define  ENV_EXTERNAL  3

/** Jet density calculation types. */
#define DENSITY_LENGTH_SCALE 0
#define DENSITY_JET_POWER    1

/** Jet injection zone shapes. */
#define CONICAL_INJECTION_ZONE   0
#define SPHERICAL_INJECTION_ZONE 1

#define SEC2MYR (3.171e-14)
#define CM2KPC (1/(1e3*(CONST_pc)))

/** Galaxy age in gigayears. */
#define GALAXY_AGE 13.0

/** Supernovae at tn. */
#define SUPERNOVAE_FACTOR 0.18

/** T0 factor. */
#define T0_FACTOR                6.42e-11

/** Jet domain settings. */
#define JD_GRAD 0 ///< Set the jet domain based on the pressure gradient (originally in PLUTO)
#define JD_PRES 1 ///< Set the jet domain based on the difference between current and initial pressure

/** Moving injection region settings. */
#define MOVING_INJECTION_NONE 0
#define MOVING_INJECTION_LINEAR 1

#ifndef MOVING_INJECTION_REGION
  #define MOVING_INJECTION_REGION MOVING_INJECTION_NONE
#endif

extern double nonrel_gamma;

#endif
