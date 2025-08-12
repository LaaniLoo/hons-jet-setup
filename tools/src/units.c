#include "pluto.h"
#include "units.h"

/* Our global variables for normalisation */
VariableNormalisation vn;

void SetNormalisation() {
    /* Set up our normalisations in our normalisation structure.
     * These handle converting from code units to cgs units
     */

    /* Unit values */
    vn.length = UNIT_LENGTH;
    vn.density = UNIT_DENSITY;
    vn.velocity = UNIT_VELOCITY;
    vn.temperature = 1; /* Temperature is always in Kelvin */

    /* Derived normalisations */
    vn.time = vn.length / vn.velocity;
    vn.area = pow(vn.length, 2);
    vn.pressure = vn.density * pow(vn.velocity, 2);
    vn.power = vn.pressure * vn.velocity * vn.area;
    vn.energy_flux = vn.pressure * vn.velocity;
    vn.energy = vn.pressure * pow(vn.length, 3);
    vn.gravitational_potential = pow(vn.velocity, 2);
    vn.acceleration = vn.velocity / vn.time;
    vn.number_density = 1.0 / pow(vn.length, 3);
    vn.mass = vn.density * pow(vn.length, 3);

    print("> Calculated normalisations ");
}

void PrintNormalisations() {
    print("\n\nCode-to-CGS Normalisations:\n\n");
    print("length                  = %16e \n", vn.length);
    print("density                 = %16e \n", vn.density);
    print("velocity                = %16e \n", vn.velocity);
    print("temperature             = %16e \n", vn.temperature);
    print("time                    = %16e \n", vn.time);
    print("area                    = %16e \n", vn.area);
    print("pressure                = %16e \n", vn.pressure);
    print("power                   = %16e \n", vn.power);
    print("energy_flux             = %16e \n", vn.energy_flux);
    print("energy                  = %16e \n", vn.energy);
    print("gravitational_potential = %16e \n", vn.gravitational_potential);
    print("acceleration            = %16e \n", vn.acceleration);
    print("number_density          = %16e \n", vn.number_density);
    print("mass                    = %16e \n", vn.mass);
    print("\n\n");
}
