#ifndef jet_h
#define jet_h

#include "pluto.h"
#include "pluto_usr.h"
#include "tools_usr.h"
#include "coordinates.h"
#include "units.h"
#include "environment.h"
#include "structs_usr.h"

/** Calculate the jet injection region properties.
 * The jet injection region properties are calculated based off the injection parameters
 * set in the PLUTO ini file. This method calculates injection geometry, jet density, and pressure.
 *
 * @param jirp The jet injection region properties structure to be populated
 * @param ep The environment properties structure. Either AnalyticEnvironmentProperties, TabulatedEnvironmentProperties, or ExternalEnvironmentProperties
 * @param jet_id Whether this is the main jet or counter jet
 */
void Jet_LoadInjectionRegionProperties(JetInjectionRegionProperties* jirp, void* ep, jet_identifier jet_id);

/** Calculates the jet intermittency properties.
 * The jet intermittency properties are calculated based off the intermittency parameters
 * set in the PLUTO ini file.
 *
 * @param jip The jet intermittency properties structure to be populated
 */
void Jet_LoadIntermittencyProperties(JetIntermittencyProperties* jip);

void Jet_PrintInjectionProperties(JetInjectionRegionProperties*, jet_identifier);
void Jet_PrintIntermittencyProperties(JetIntermittencyProperties*);

void Jet_Enable(JetIntermittencyProperties*);
void Jet_Disable(JetIntermittencyProperties*);

void Jet_UpdateIntermittency(JetIntermittencyProperties*);
void Jet_UpdateInjectionRegionPosition(JetInjectionRegionProperties* []);

void Jet_TransformToNozzleCoordinates(double, double, double, JetInjectionRegionProperties* [], JetTransformedCoordinates*);
void Jet_WithinNozzle(JetTransformedCoordinates*, int*, int*);
void Jet_FillQuantitiesArray(JetTransformedCoordinates*, JetIntermittencyProperties*, double []);

void Jet_NozzleToGrid(JetInjectionRegionProperties*, const double, const double, const double, double*, double*, double*);

/** Calculates the velocity components for a specific cell within the jet injection zone.
 * The components are calculated based on whether the jet is collimated, conical, or
 * using a spherical jet injection region. If necessary the velocity vectors are 
 * rotated to match the rotation of the jet cone.
 *
 * @param jirp Jet injection region properties structure
 * @param jtc Jet transformed coordinates structure
 * @param radial_velocity Velocity along the jet axis
 * @param vj State vector for the jet quantities
 */
void Jet_InjectionZoneVelocity(JetInjectionRegionProperties* jirp, JetTransformedCoordinates* jtc, double radial_velocity, double vj[]);

/** Calculates the density for a specific cell within the jet injection zone.
 *
 * @param jirp Jet injection region properties structure
 * @param jtc Jet transformed coordinates structure
 * @returns density The jet density within the injection zone
 */
double Jet_InjectionZoneDensity(JetInjectionRegionProperties* jirp, JetTransformedCoordinates* jtc);

/** Calculates the pressure for a specific cell within the jet injection zone.
 *
 * @param jirp Jet injection region properties structure
 * @param vj State vector for the jet quantities
 */
void Jet_InjectionZonePressure(JetInjectionRegionProperties* jirp, double vj[]);

int Jet_WithinCollimatedInjectionRegion(JetInjectionRegionProperties*, double x1, double x2, double x3, jet_identifier jet_id);

#endif
