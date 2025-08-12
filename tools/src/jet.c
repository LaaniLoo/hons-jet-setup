#include "jet.h"
#include "coordinates.h"
#include "pluto_usr.h"
#include "structs_usr.h"
#include "tools_usr.h"

void Jet_LoadInjectionRegionProperties(JetInjectionRegionProperties* jirp, void* ep, jet_identifier jet_id) {
    double jet_angle_temp = 0;

    double jet_gm = 5./3.;
    double jet_gmm1 = jet_gm / (jet_gm - 1);

    double env_density_at_inlet = 0;

    /* Set up the jet opening angle */
    if (jet_id == MAIN_JET) 
        jet_angle_temp = g_inputParam[JET_OA_PRIMARY];
    else if (jet_id == COUNTER_JET) 
        jet_angle_temp = g_inputParam[JET_OA_SECONDARY];

    jirp->speed = g_inputParam[JET_SPD];
    jirp->power = g_inputParam[JET_PWR] / vn.power;
    jirp->half_opening_angle = jet_angle_temp * CONST_PI / 180.0;
    jirp->initial_width = g_inputParam[JET_INITIAL_RADIUS];
    jirp->injection_region_radius = g_inputParam[JET_INJECTION_HEIGHT];
    jirp->injection_rotation_angle = g_inputParam[JET_ROTATION_ANGLE] * CONST_PI / 180.0;
    jirp->pressure_factor = 1.0;
    jirp->chi = g_inputParam[JET_CHI];

    if (jirp->pressure_factor == 0) {
        print("JET_PRS_FACTOR cannot be 0! Exiting...\n");
        QUIT_PLUTO(1);
    }

    jirp->collimated = (jirp->half_opening_angle <= 1e-12) ? 1 : 0;

    /* Treat g_inputParam[JET_SPD] as lorentz factor if > speed of light */
    if (g_inputParam[JET_SPD] * UNIT_VELOCITY > CONST_c) {
        jirp->lorentz = g_inputParam[JET_SPD];
        jirp->speed = lorentz2speed(jirp->lorentz);

        print("> Treating JET_SPD as lorentz factor %f (%f c), as it is greater than the speed of light\n", jirp->lorentz, jirp->speed);
    }

    /* Load the jet offset parameters */
    jirp->initial_offset[0] = g_inputParam[JET_X1O];
    jirp->initial_offset[1] = g_inputParam[JET_X2O];
    jirp->initial_offset[2] = g_inputParam[JET_X3O];

    /* Set inlet speed parameters */
    /* jirp->inlet_speed[0] = 100.0 * (SEC2MYR * vn.time); */
    jirp->inlet_speed[0] = 0.0;
    jirp->inlet_speed[1] = 0.0;
    jirp->inlet_speed[2] = 0.0;

    /* Set up the jet injection cone parameters */
    if (jirp->collimated) {
        jirp->injection_cone_offset = 0;
    } else {
        if (jet_id == MAIN_JET)
            jirp->injection_cone_offset = jirp->initial_width / tan(jirp->half_opening_angle);
        else if (jet_id == COUNTER_JET)
            jirp->injection_cone_offset = jirp->initial_width / tan(CONST_PI - jirp->half_opening_angle);
#if JET_INJECTION_CAP == NO
        jirp->injection_cone_height = jirp->injection_region_radius * cos(jirp->half_opening_angle);
#endif
    }

    /* We also calculate the environment density at the inlet, which is used to caculate the 
     * eta parameter (ratio of jet density to environment density) */
    env_density_at_inlet = density_at_position(-jirp->offset[0]+0.001, -jirp->offset[1], -jirp->offset[2], ep);

    if (!jirp->collimated) {
        jirp->solid_angle = 2.0 * CONST_PI * (1.0-cos(jirp->half_opening_angle));
        jirp->L1 = 2.0 * sqrt(2.0) * sqrt((jirp->power) / (env_density_at_inlet * pow(jirp->speed, 3)));
        jirp->L1b = jirp->L1 * 0.5/sqrt(jirp->solid_angle);

#if JET_INJECTION_ZONE == CONICAL_INJECTION_ZONE
        print("> Using a conical injection zone\n");
        /* Because we offset our jets, the area actually changes. We need to account for this here: */
        jirp->area = jirp->solid_angle * pow(jirp->injection_region_radius + fabs(jirp->injection_cone_offset), 2);
#elif JET_INJECTION_ZONE == SPHERICAL_INJECTION_ZONE
        print("> Using a spherical injection zone\n");
        jirp->injection_cone_offset = 0.0;
        jirp->area = jirp->solid_angle * pow(jirp->injection_region_radius, 2);
#endif
    }
    else {
        jirp->area = CONST_PI * pow(jirp->initial_width, 2);
    }


    /* Calculate the jet density
     * Use either environment length scales,
     * or desired jet power */
#if JET_DENSITY_CALCULATION == DENSITY_LENGTH_SCALE
    if (jirp->collimated) {
        print("JET_DENSITY_CALCULATION == DENSITY_LENGTH_SCALE is only supported with uncollimated jets\n");
        QUIT_PLUTO(1);
    }
    jirp->density = density_at_position(jirp->L1b, 0, 0, ep);
    jirp->density_contrast = jirp->enthalpy_contrast = jirp->energy_density_constrast = jirp->density / env_density_at_inlet;

    /* In this case we get the pressure from the environment, scaled by the CHI parameter */
    jirp->pressure = pressure_at_position(-jirp->offset[0]+0.001, -jirp->offset[1], -jirp->offset[2], ep) * jirp->chi;
    print("> Setting density with length scale\n");
    print("  WARNING: Over-pressure factor is taken to be the CHI parameter\n");
#elif JET_DENSITY_CALCULATION == DENSITY_JET_POWER

    #if (PHYSICS == RHD || PHYSICS == RMHD) && JET_RELATIVISTIC_DENSITY == YES
        jirp->lorentz = speed2lorentz(jirp->speed);
        jirp->density = (jirp->power) / ((jirp->speed * jirp->area * pow(CONST_c / UNIT_VELOCITY, 2)) * (jirp->lorentz * (jirp->lorentz - 1) + (jirp->lorentz * jirp->lorentz) / (jirp->chi)));
        jirp->pressure = jirp->density * pow(CONST_c / UNIT_VELOCITY, 2) / (jirp->chi * jet_gmm1);
        jirp->enthalpy = pow(CONST_c / UNIT_VELOCITY, 2) + jet_gmm1 * jirp->pressure / jirp->density;
        jirp->density_contrast = jirp->density / env_density_at_inlet;
        jirp->enthalpy_contrast = jirp->density * jirp->enthalpy  / (env_density_at_inlet * pow(CONST_c / UNIT_VELOCITY, 2));
        jirp->energy_density_constrast = jirp->enthalpy_contrast * jirp->lorentz * jirp->lorentz;
        print("> Setting density with jet power (relativistic case)\n");
    #elif JET_RELATIVISTIC_DENSITY == NO
        jirp->density = 2.0* (jirp->power) / (pow(jirp->speed, 3) * jirp->area);
        jirp->density_contrast = jirp->enthalpy_contrast = jirp->energy_density_constrast = jirp->density / env_density_at_inlet;

        /* In this case we get the pressure from the environment, scaled by the CHI parameter */
        jirp->pressure = pressure_at_position(-jirp->offset[0]+0.001, -jirp->offset[1], -jirp->offset[2], ep) * jirp->chi;
        print("> Setting density with jet power (non-relativistic case)\n");
        print("  WARNING: Over-pressure factor is taken to be the CHI parameter\n");
    #elif JET_RELATIVISTIC_DENSITY == YES && !(PHYSICS == RHD || PHYSICS == RMHD)
        print("Cannot use relativistic density calculation in a non-relativistic simulation. Exiting...\n");
        QUIT_PLUTO(1);
    #endif
#endif
    print("> Loaded jet %i injection region properties\n", jet_id);
}

void Jet_LoadIntermittencyProperties(JetIntermittencyProperties* jip) {
    jip->delta_t = 0.05 / (SEC2MYR * vn.time);
    jip->start_time = g_inputParam[JET_START_TIME] / (SEC2MYR * vn.time);
    jip->end_time = g_inputParam[JET_END_TIME] / (SEC2MYR * vn.time);
    jip->total_active_time = g_inputParam[JET_ACTIVE_TIME] / (SEC2MYR * vn.time);
    jip->episodes = (int)round(g_inputParam[JET_EPISODES]);

    jip->active_time = jip->total_active_time / jip->episodes;
    jip->quiet_time = (jip->end_time - jip->total_active_time) / jip->episodes;
    jip->duty_cycle_time = jip->active_time + jip->quiet_time;

    jip->is_active = 0;
    jip->current_episode = 0;
    jip->have_set_timestep = 0;
    jip->intermittency_velocity_factor = 1.0;

    print("> Loaded jet intermittency properties\n");
}

void Jet_PrintInjectionProperties(JetInjectionRegionProperties* jirp, jet_identifier jet_id) {
    print("\n\nJet %i injection region properties:\n\n", jet_id);
    print("initial_width            = %16e [%16e kpc] \n", jirp->initial_width, jirp->initial_width * vn.length * CM2KPC);
    print("half_opening_angle       = %16e \n", jirp->half_opening_angle);
    print("power                    = %16e [%16e erg/s] \n", jirp->power, jirp->power * vn.power);
    print("speed                    = %16e [%16e c] \n", jirp->speed, jirp->speed * vn.velocity / CONST_c);
    print("lorentz factor           = %16e \n", jirp->lorentz);
    print("chi                      = %16e \n", jirp->chi);
    print("density                  = %16e [%16e g/cm^3] \n", jirp->density, jirp->density * vn.density);
    print("pressure                 = %16e [%16e barye] \n", jirp->pressure, jirp->pressure * vn.pressure);
    print("specific enthalpy        = %16e [%16e erg/g] \n", jirp->enthalpy, jirp->enthalpy * vn.gravitational_potential);
    print("proper density contrast  = %16e \n", jirp->density_contrast);
    print("enthalpy contrast        = %16e \n", jirp->enthalpy_contrast);
    print("energy density contrast  = %16e \n", jirp->energy_density_constrast);
    print("L1                       = %16e [%16e kpc] \n", jirp->L1, jirp->L1 * vn.length * CM2KPC);
    print("L1b                      = %16e [%16e kpc] \n", jirp->L1b, jirp->L1b * vn.length * CM2KPC);
    print("collimated               = %16i \n", jirp->collimated);
    print("solid_angle              = %16e [%16e kpc^2] \n", jirp->solid_angle, jirp->solid_angle * vn.area * CM2KPC * CM2KPC);
    print("area                     = %16e [%16e kpc^2] \n", jirp->area, jirp->area * vn.area * CM2KPC * CM2KPC);
    print("injection_rotation_angle = %16e \n", jirp->injection_rotation_angle);
    print("injection_region_radius  = %16e [%16e kpc] \n", jirp->injection_region_radius, jirp->injection_region_radius * vn.length * CM2KPC);
    print("injection_cone_height    = %16e [%16e kpc] \n", jirp->injection_cone_height, jirp->injection_cone_height * vn.length * CM2KPC);
    print("injection_cone_offset    = %16e [%16e kpc] \n", jirp->injection_cone_offset, jirp->injection_cone_offset * vn.length * CM2KPC);
    print("X1 offset                = %16e [%16e kpc] \n", jirp->offset[0], jirp->offset[0] * vn.length * CM2KPC);
    print("X2 offset                = %16e [%16e kpc] \n", jirp->offset[1], jirp->offset[1] * vn.length * CM2KPC);
    print("X3 offset                = %16e [%16e kpc] \n", jirp->offset[2], jirp->offset[2] * vn.length * CM2KPC);
    print("Pressure factor          = %16e \n", jirp->pressure_factor);
    print("\n\n");
}

void Jet_PrintIntermittencyProperties(JetIntermittencyProperties* jip) {
    print("\n\nJet intermittency properties:\n\n");
    print("delta_t           = %16e [%0.02f Myr] \n", jip->delta_t, jip->delta_t * vn.time * SEC2MYR);
    print("episodes          = %16i \n", jip->episodes);
    print("start_time        = %16e [%0.02f Myr] \n", jip->start_time, jip->start_time * vn.time * SEC2MYR);
    print("end_time          = %16e [%0.02f Myr] \n", jip->end_time, jip->end_time * vn.time * SEC2MYR);
    print("total_active_time = %16e [%0.02f Myr] \n", jip->total_active_time, jip->total_active_time * vn.time * SEC2MYR);
    print("active_time       = %16e [%0.02f Myr] \n", jip->active_time, jip->active_time * vn.time * SEC2MYR);
    print("quiet_time        = %16e [%0.02f Myr] \n", jip->quiet_time, jip->quiet_time * vn.time * SEC2MYR);
    print("duty_cycle_time   = %16e [%0.02f Myr] \n", jip->duty_cycle_time, jip->duty_cycle_time * vn.time * SEC2MYR);
    print("\n\n");
}

void Jet_Enable(JetIntermittencyProperties* jip) {
    if (!jip->is_active) {
        jip->is_active = 1;
        print("\n> Switched on jet outburst %i/%i\n\n", jip->current_episode+1, jip->episodes);
        
        /* Set global reset timestep flag */
        g_resetTimeStep = 1;
        jip->have_set_timestep = 1;
    }
}

void Jet_Disable(JetIntermittencyProperties* jip) {
    if (jip->is_active) {
        jip->is_active = 0;
        print("\n> Switched off jet outburst %i/%i\n\n", jip->current_episode+1, jip->episodes);

        jip->have_set_timestep = 0;
    }
}

void Jet_UpdateIntermittency(JetIntermittencyProperties* jip) {
    double jet_time;
    double next_jet_time;
    double next_g_time = g_time + g_dt;
    int ramp_up, ramp_down;

    /* Calculate the jet time relative to our start time 
     * Note that we wrap the jet_time between 0 & jip->duty_cycle_time.
     * For this to work, we need to add a separate check later on for
     * g_time >= jip->start_time to make sure we don't turn the jet on too early */
    jet_time = fmod(g_time - jip->start_time + jip->duty_cycle_time, jip->duty_cycle_time);
    next_jet_time = fmod(next_g_time - jip->start_time + jip->duty_cycle_time, jip->duty_cycle_time);

    ramp_up = (jet_time >= (jip->duty_cycle_time - jip->delta_t));
    /* Need to also check for ramping up at the start of the initial outburst, if it's not at the start of the simulation */
    ramp_up |= (jip->start_time - g_time < jip->delta_t && jip->start_time > g_time);
    ramp_down = (jet_time >= jip->active_time && jet_time <= jip->active_time + jip->delta_t);

    /* Set timestep reset flag if we're coming up on an injection */
    if (
            (
             (jip->start_time - next_g_time < jip->delta_t && jip->start_time > next_g_time) 
             || (next_jet_time >= (jip->duty_cycle_time - jip->delta_t))
            )
            && (!jip->have_set_timestep)
            ) {
        g_resetTimeStep = 1;
        jip->have_set_timestep = 1;
    }

    /* Check if we have started our injection cycle */
    if (jet_time >= 0 && g_time >= jip->start_time) {
        /* Calculate which jet episode we are in */
        jip->current_episode = (int)floor(jet_time / jip->duty_cycle_time);
        
        /* Set our velocity factor */
        jip->intermittency_velocity_factor = 1.0;

        /* Check if we have just switched off the jet */
        if (jet_time >= jip->active_time + jip->delta_t && jet_time <= (jip->duty_cycle_time - jip->delta_t)) {
            Jet_Disable(jip);
        }
        /* Check if we have just switched on the jet */
        else {
            Jet_Enable(jip);
        }

        /* Ramp the jet down */
        if (ramp_down) {
            jip->intermittency_velocity_factor = cos((CONST_PI / 2) * ((jet_time - jip->active_time) / jip->delta_t));
        }
    }

    /* Even if we haven't started the injection cycle, we should still check if we're ramping up to the very first outburst */
    if (ramp_up) {
        jip->intermittency_velocity_factor = (1 - cos((CONST_PI / 2) * ((jet_time - (jip->duty_cycle_time - jip->delta_t)) / jip->delta_t)));
        Jet_Enable(jip);
    }
}

void Jet_UpdateInjectionRegionPosition(JetInjectionRegionProperties* jirp[]) {
    int ii, jj;

#if MOVING_INJECTION_REGION == MOVING_INJECTION_LINEAR
    for (ii = 0; ii < 2; ii++) {
        for (jj = 0; jj < 3; jj++) {
            jirp[ii]->offset[jj] = jirp[ii]->initial_offset[jj] + jirp[ii]->inlet_speed[jj] * g_time;
        }
    }
#elif MOVING_INJECTION_REGION == MOVING_INJECTION_NONE
    for (ii = 0; ii < 2; ii++) {
        for (jj = 0; jj < 3; jj++) {
            jirp[ii]->offset[jj] = jirp[ii]->initial_offset[jj];
        }
    }
#endif
    return;
}

void Jet_TransformToNozzleCoordinates(double x1, double x2, double x3, JetInjectionRegionProperties* jirp[], JetTransformedCoordinates* jtc) {
    double x1_prime, x2_prime, x3_prime;
    double x1i, x2j, x3k;
    int is_main_jet;

    x1i = x1;
    x2j = x2;
    x3k = x3;
    /* Rotate if needed */
#if JET_INJECTION_ROTATED == YES
    x1i = x1*cos(jirp[0]->injection_rotation_angle) - x3*sin(jirp[0]->injection_rotation_angle);
    x3k = x1*sin(jirp[0]->injection_rotation_angle) + x3*cos(jirp[0]->injection_rotation_angle);
#endif

    /* Our final rotated, offset variables */
    x1_prime = x1i + jirp[0]->offset[0];
    x2_prime = x2j + jirp[0]->offset[1];
    x3_prime = x3k + jirp[0]->offset[2];

    /* calculate r, theta and velocity depending on the coordinate system in use */
    jtc->r = sph_radius(x1_prime, x2_prime, x3_prime);
    jtc->r_from_origin = sph_radius(x1i, x2j, x3k);
    jtc->phi = sph_phi(x1_prime, x2_prime, x3_prime);

    /* Let's find out which jet we are */
#if GEOMETRY == CARTESIAN
    is_main_jet = (SELECT(0, x2_prime, x3_prime) > 0);
#elif GEOMETRY == SPHERICAL
    is_main_jet = x2_prime < (CONST_PI * 0.5);
#endif

    if (is_main_jet) { // First jet
        jtc->jet_injection_properties = jirp[0];
        jtc->jet_id = MAIN_JET;
    }
    else { // Second jet
        jtc->jet_injection_properties = jirp[1];
        jtc->jet_id = COUNTER_JET;
    }

#if GEOMETRY == CARTESIAN           // 2D: x1=x, x2=z, 3D: x1=x, x2=y, x3=z
#if COMPONENTS == 2
    jtc->r_offset = sqrt(x1_prime*x1_prime + (x2_prime+current_jet_properties->injection_cone_offset)*(x2_prime+current_jet_properties->injection_cone_offset));
#elif COMPONENTS == 3
    jtc->r_offset = sqrt(x1_prime*x1_prime + x2_prime*x2_prime + (x3_prime+jtc->jet_injection_properties->injection_cone_offset)*(x3_prime+jtc->jet_injection_properties->injection_cone_offset));
#endif
    jtc->theta = acos((SELECT(x3_prime, x2_prime+jtc->jet_injection_properties->injection_cone_offset, x3_prime+jtc->jet_injection_properties->injection_cone_offset))/jtc->r_offset);
    /* jtc->phi = atan2(SELECT(0.0, 0.0, x2_prime), SELECT(0.0, 0.0, x1_prime)); */
#elif GEOMETRY == POLAR             // x1=r (or rho), x2=phi (not theta...)
    jtc->theta = x2_prime;
#elif GEOMETRY == SPHERICAL         // x1=r, x2=theta, x3=phi
    jtc->theta = x2_prime;
    jtc->phi = SELECT(0, 0, x3_prime);
#endif
}

void Jet_WithinNozzle(JetTransformedCoordinates* jtc, int* within_r, int* within_theta) {
    /* Initialise return variables */
    *within_r = 0;
    *within_theta = 0;

    /* Check we are within theta */
    if (jtc->jet_injection_properties->collimated) {
        *within_r = fabs(jtc->r * cos(jtc->theta)) <= jtc->jet_injection_properties->injection_region_radius;
        *within_theta = fabs(jtc->r * sin(jtc->theta)) <= jtc->jet_injection_properties->initial_width;
    } else {
        /* -- Check that we are within the desired r range -- */
#if JET_INJECTION_CAP == YES
        *within_r = (jtc->r <= jtc->jet_injection_properties->injection_region_radius);
#elif JET_INJECTION_CAP == NO
        *within_r = (fabs(jtc->r * cos(jtc->theta)) <= jtc->jet_injection_properties->injection_cone_height);
#endif
        if (jtc->jet_id == MAIN_JET) {
            *within_theta = (jtc->theta <= jtc->jet_injection_properties->half_opening_angle && jtc->theta >= -jtc->jet_injection_properties->half_opening_angle);
        } else if (jtc->jet_id == COUNTER_JET) {
            *within_theta = (jtc->theta <= (jtc->jet_injection_properties->half_opening_angle + CONST_PI) && jtc->theta >= (CONST_PI - jtc->jet_injection_properties->half_opening_angle));
        }
    }
}

void Jet_FillQuantitiesArray(JetTransformedCoordinates* jtc, JetIntermittencyProperties* jip, double vj[]) {
    double radial_velocity;
    /* Calculate the desired jet radial velocity */
    radial_velocity = jip->intermittency_velocity_factor * jtc->jet_injection_properties->speed;

    /* Get velocity */
    Jet_InjectionZoneVelocity(jtc->jet_injection_properties, jtc, radial_velocity, vj);
    /* Get density */
    vj[RHO] = Jet_InjectionZoneDensity(jtc->jet_injection_properties, jtc);
    /* Get pressure */
    Jet_InjectionZonePressure(jtc->jet_injection_properties, vj);
}

void Jet_NozzleToGrid(JetInjectionRegionProperties* jirp, const double r, const double theta, const double phi,
                      double* x1, double* x2, double* x3) {
#if GEOMETRY == CARTESIAN
    double tx1, tx3;
    sph2cart(r + fabs(jirp->injection_cone_offset), theta, phi, x1, x2, x3);
    *x1 -= jirp->offset[0];
    *x2 -= jirp->offset[1];
    *x3 -= jirp->offset[2] + jirp->injection_cone_offset;

#if JET_INJECTION_ROTATED == YES
    tx1 = *x1;
    tx3 = *x3;
    *x1 = tx1*cos(-jirp->injection_rotation_angle) - tx3*sin(-jirp->injection_rotation_angle);
    *x3 = tx1*sin(-jirp->injection_rotation_angle) + tx3*cos(-jirp->injection_rotation_angle);
#endif
#elif GEOMETRY == SPHERICAL
    *x1 = r - jirp->offset[0];
    *x2 = theta;
    *x3 = SELECT(0, 0, phi);
#endif
}

inline void Jet_InjectionZoneVelocity(JetInjectionRegionProperties* jirp, JetTransformedCoordinates* jtc, double radial_velocity, double vj[]) {
    double tvx1, tvx3;
#if GEOMETRY == CARTESIAN           // x1=x, x2=y
    /* Check whether we are injecting a collimated or uncollimated jet */
    if (jtc->jet_injection_properties->collimated) {
        EXPAND( vj[VX1] = 0.0;,
                vj[VX2] = SELECT(0, radial_velocity * (jtc->jet_id == MAIN_JET ? 1 : -1), 0);,
                vj[VX3] = radial_velocity * (jtc->jet_id == MAIN_JET ? 1 : -1); )
    } else {
    EXPAND( vj[VX1] = radial_velocity*sin(jtc->theta)*cos(jtc->phi);,
            vj[VX2] = SELECT(0, radial_velocity*cos(jtc->theta), radial_velocity*sin(jtc->theta)*sin(jtc->phi));,
            vj[VX3] = radial_velocity*cos(jtc->theta); )
    }

#if JET_INJECTION_ROTATED == YES
    tvx1 = vj[VX1]*cos(-jtc->jet_injection_properties->injection_rotation_angle) - vj[VX3]*sin(-jtc->jet_injection_properties->injection_rotation_angle);
    tvx3 = vj[VX1]*sin(-jtc->jet_injection_properties->injection_rotation_angle) + vj[VX3]*cos(-jtc->jet_injection_properties->injection_rotation_angle);
    vj[VX1] = tvx1;
    vj[VX3] = tvx3;
#endif

#elif GEOMETRY == POLAR             // x1=r (or rho), x2=phi (not theta...)
    vj[VX1] = radial_velocity;
    vj[VX2] = 0;
#elif GEOMETRY == SPHERICAL         // x1=r, x2=theta
    vj[VX1] = radial_velocity;
    vj[VX2] = 0;
#endif

#if RECONSTRUCT_4VEL == YES
    EXPAND( vj[VX1] *= jtc->jet_injection_properties->lorentz;,
            vj[VX2] *= jtc->jet_injection_properties->lorentz;,
            vj[VX3] *= jtc->jet_injection_properties->lorentz;);
#endif
}

inline double Jet_InjectionZoneDensity(JetInjectionRegionProperties* jirp, JetTransformedCoordinates* jtc) {
    /* The jet density falls off as 1/r^2 along the jet injection cone.
     * Depending on how the density is defined, we calculate the actual density
     * within the injection cone using a bunch of methods. */
#if JET_DENSITY_CALCULATION == DENSITY_LENGTH_SCALE
    /* Using length scales, the jet density is set to the environment density at
     * the length scale L1b. The offset radius is used within the cone, to
     * properly account for the cone offset & initial widht */
    return jirp->density * pow(jirp->L1b / jtc->r_offset, 2);
#elif JET_DENSITY_CALCULATION == DENSITY_JET_POWER
    /* Using the density calculated from jet power, we have three options.
     * If the jet is collimated, then the density is equal everywhere in the injection
     * zone. Otherwise the density varies depending on whether a conical or spherical injection
     * zone is used. */
    if (jirp->collimated) return jirp->density;
#if JET_INJECTION_ZONE == CONICAL_INJECTION_ZONE
    /* The conical case is similar to the length scale calculation, making sure to use
     * the correct offset radii. The density is limited to be 1e4 x the jet density at the 
     * head of the injection zone. */
    return jirp->density * fmin(pow((jirp->injection_region_radius + fabs(jirp->injection_cone_offset)) / jtc->r_offset, 2), 1e4);
#elif JET_INJECTION_ZONE == SPHERICAL_INJECTION_ZONE
    /* The spherical injection case has the density falling off as (1+(r/r0)^2)^-1.
     * In this case rho0 = 2rhoj, to ensure that at the injection radius r0, the
     * density is equal to the desired jet density. */
    return 2 * jirp->density / (1 + pow(jtc->r / jirp->injection_region_radius, 2));
#endif
#endif
};

inline void Jet_InjectionZonePressure(JetInjectionRegionProperties *jirp, double vj[]) {
    /* The jet pressure falls off with radius assuming adiabatic expansion.
     * In this case P0 = 2^gamma Pj, to ensure that at the injection radius r0, the
     * pressure is equal to the desired jet pressure*/
    double jet_gm = 5./3.;
    vj[PRS] = pow(2, jet_gm) * jirp->pressure * pow(vj[RHO] / (2*jirp->density), jet_gm);
}

inline int Jet_WithinCollimatedInjectionRegion(JetInjectionRegionProperties* jirp, double x1, double x2, double x3, jet_identifier jet_id) {
    if (jirp->collimated) {
        double c_r = cyl_radius(x1, x2, x3);
        double c_z = cyl_z(x1, x2, x3);
        if (jet_id == COUNTER_JET) c_z *= -1.0;
        return (c_r < jirp->initial_width) && (c_z < jirp->injection_region_radius);
    }
    return 0;
}
