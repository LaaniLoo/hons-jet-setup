#include "pluto.h"
#include "pluto_usr.h"
#include "units.h"
#include "coordinates.h"
#include "jet.h"
#include "environment.h"
#include "galaxy.h"
#include "tools_usr.h"

/* Set up our density profile specific code */
#if DENSITY_PROFILE == ENV_TABULATED
/* Tabulated density profile has its own environment properties structure */
TabulatedEnvironmentProperties ep;
TabulatedEnvironmentProperties secondary_ep;
#elif DENSITY_PROFILE == ENV_EXTERNAL
/* The external environment profile also has its own environment properties structure */
ExternalEnvironmentProperties ep;
#else
/* Just your basic, garden variety, King or Makino environment */
AnalyticEnvironmentProperties ep;
#endif

JetIntermittencyProperties jip;
JetInjectionRegionProperties jirp1, jirp2;
JetInjectionRegionProperties* jirp_array[2];

int done_once = 0;
int smoothed_sim_region = 0;

double nonrel_gamma = 5. / 3.;

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  /* We have some setup to do initially */
  if (done_once == 0)
  {
      /* Seed our PRNG */
      RandomSeed(SeedGenerator(), 0);
      /* Let's calculate our normalisations */
      SetNormalisation();

      /* ..and print them to the log */
      PrintNormalisations();

      // set gamma variables
#if EOS == IDEAL
      g_gamma = 5.0/3.0;
#endif

      /* -- THIS IS WHERE WE SET UP OUR ENVIRONMENTS -- */
#if DENSITY_PROFILE == ENV_EXTERNAL
      setup_external_environment(&ep);
      set_external_env_paths(&ep, "spherical-r200-n3.dbl", "spherical-r200-n3.grid.out");

      /* strcpy(ep.external_env_data, "spherical-r200-n3.dbl"); */
      /* strcpy(ep.external_env_out, "spherical-r200-n3.grid.out"); */

      /* PrintEnvironmentProperties(&ep); */
#elif DENSITY_PROFILE == ENV_MAKINO || DENSITY_PROFILE == ENV_KING || DENSITY_PROFILE == ENV_CONST
      setup_analytic_environment(&ep);

      /* PrintEnvironmentProperties(&ep); */
#elif DENSITY_PROFILE == ENV_TABULATED
      set_tabulated_env_path(&ep, "./data/PLUTO_data_best_fGas=0.010_Mh=4.00E+13_500m_ccc+gal_Mgal=1.00E+12galPot=true_pR=1.45E-17.csv");

      /* strcpy(ep.file_name, "./data/PLUTO_data_best_fGas=0.010_Mh=4.00E+13_500m_ccc+gal_Mgal=1.00E+12galPot=true_pR=1.45E-17.csv"); */
      /* strcpy(secondary_ep.file_name, "extra-galaxy.csv"); */
      secondary_ep.offset[0] = 100;
#endif

      PrintEnvironmentProperties(&ep);

      set_default_environment(&ep);

      Jet_LoadIntermittencyProperties(&jip);
      Jet_PrintIntermittencyProperties(&jip);

      /* Update the jet intermittency properties */
      Jet_UpdateIntermittency(&jip);


      Jet_LoadInjectionRegionProperties(&jirp1, &ep, MAIN_JET);
      Jet_PrintInjectionProperties(&jirp1, MAIN_JET);

      Jet_LoadInjectionRegionProperties(&jirp2, &ep, COUNTER_JET);
      Jet_PrintInjectionProperties(&jirp2, COUNTER_JET);

      jirp_array[0] = &jirp1;
      jirp_array[1] = &jirp2;

      /* Update the jet injection region position */
      Jet_UpdateInjectionRegionPosition(jirp_array);

      Particles_SetJetProperties(&jirp1, &jirp2, &jip);

      /* Set up some cooling variables */
#if COOLING
      g_minCoolingTemp = 1e4;
#endif

      /* Print out debugging information for supernovae feedback */

      double r_eff = get_r_eff_elliptical(g_inputParam[STELLAR_MASS]);
      double alpha = get_alpha(GALAXY_AGE);
      double T0 = get_T0();
      double cV = get_cV();
      double stellar_density = get_stellar_density_hernquist(g_inputParam[STELLAR_MASS], 0.1, 0, 0);
      double energy_inj_rate = get_energy_injection_rate(stellar_density);
      double density_loss_rate = get_mass_injection_rate(stellar_density);

      print("Supernovae feedback properties (at r=0.1): \n\n");
      print("Effective radius      = %16e kpc\n", r_eff);
      print("Alpha                 = %16e s^-1 [%16e code units]\n", alpha, alpha * vn.time);
      print("T0                    = %16e\n", T0);
      print("cV                    = %16e\n", cV);
      print("Stellar density       = %16e g cm^-3 [%16e code units]\n", stellar_density, stellar_density / vn.density);
      print("Energy injection rate = %16e erg cm^-3 s^-1 [%16e code units]\n", energy_inj_rate, energy_inj_rate * pow(vn.length, 3) / vn.power);
      print("Density loss rate     = %16e g cm^-3 s^-1 [%16e code units]\n\n", density_loss_rate, density_loss_rate / vn.density * vn.time);

      double nfw_accel = nfw_gravitational_acceleration(0.1, 0, 0);
      double hern_accel = hernquist_gravitational_acceleration(g_inputParam[STELLAR_MASS], 0.1, 0, 0);

      print("Gravitational acceleration properties (at r=0.1): \n\n");
      print("NFW accel  = %16e\n", nfw_accel);
      print("hern_accel = %16e\n", hern_accel);

      nfw_accel = nfw_gravitational_acceleration(10, 0, 0);
      hern_accel = hernquist_gravitational_acceleration(g_inputParam[STELLAR_MASS], 10, 0, 0);

      print("Gravitational acceleration properties (at r=10): \n\n");
      print("NFW accel  = %16e\n", nfw_accel);
      print("hern_accel = %16e\n", hern_accel);

      /* We've done this now */
      done_once = 1;
  }

#if FALSE
  double extra_gal_r = calculate_offset_radius(x1, x2, x3, &secondary_ep);
  if (extra_gal_r < 50) {
      v[RHO] = fmax(density_at_position(x1, x2, x3, &secondary_ep), density_at_position(x1, x2, x3, &ep));
      v[PRS] = fmax(pressure_at_position(x1, x2, x3, &secondary_ep), pressure_at_position(x1, x2, x3, &ep));
  }
#endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
  int i, j, k;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
  double v[3];

  /* Let's set up our computational domain now */
  TOT_LOOP(k,j,i){

    /* Initialise density */
    d->Vc[RHO][k][j][i] = density_at_position(x1[i], x2[j], x3[k], &ep);

    /* Initialise pressure */
    d->Vc[PRS][k][j][i] = pressure_at_position(x1[i], x2[j], x3[k], &ep);

    /* Initialise our velocities */
    velocity_at_position(x1[i], x2[j], x3[k], v, &ep);
    EXPAND( d->Vc[VX1][k][j][i] = v[0];,
            d->Vc[VX2][k][j][i] = v[1];,
            d->Vc[VX3][k][j][i] = v[2]; )

    /* Initialise our tracer values */
    int ii = 0;
    NTRACER_LOOP(ii) d->Vc[ii][k][j][i] = 0.0;
  }
}


/* **************************************************************** */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 **************************************************************** */
{ }
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 *********************************************************************** */
{
    // Variable definitions
    double *x1, *x2, *x3;
    int   i, j, k, nv;
    double vj[NVAR], environment_values[NVAR];
    double v[3];
    double R, Cr, smoothing_factor;
    int smoothing_power = 6;


    // Get spatial arrays
    x1 = grid->xgc[IDIR];
    x2 = grid->xgc[JDIR];
    x3 = grid->xgc[KDIR];

    /* Apply boundary conditions if desired in spherical coordinates */
#if GEOMETRY == SPHERICAL
    if (side == X1_BEG) {
        X1_BEG_LOOP(k, j, i) {

            /* We can have either reflective or outflow boundary conditions.
             * More cases can easily be added if desired - see pluto/Src/boundary.c
             * for the inbuilt implementations. */
#if SPH_RADIUS_LOWER_BOUNDARY == REFLECTIVE
            VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IBEG-i-1];
            d->Vc[VX1][k][j][i] *= -1.0;
#elif SPH_RADIUS_LOWER_BOUNDARY == OUTFLOW
            VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG];
#endif
        }
    }
#endif


#if GEOMETRY == CARTESIAN && DENSITY_PROFILE == ENV_EXTERNAL
    if (side != 0) {
        BOX_LOOP(box, k, j, i) {
            d->Vc[RHO][k][j][i] = density_at_position(x1[i], x2[j], x3[k], &ep);
            d->Vc[PRS][k][j][i] = pressure_at_position(x1[i], x2[j], x3[k], &ep);
            velocity_at_position(x1[i], x2[j], x3[k], v, &ep);
            EXPAND( d->Vc[VX1][k][j][i] = v[0];,
                    d->Vc[VX2][k][j][i] = v[1];,
                    d->Vc[VX3][k][j][i] = v[2]; )
        }
    }
#endif

    /* -- Internal boundary -- */
    if (side == 0 || side == X1_BEG)
    {
        int ii = 0;
        /* Update the jet intermittency properties */
        Jet_UpdateIntermittency(&jip);

        /* Update the jet injection region position */
        Jet_UpdateInjectionRegionPosition(jirp_array);

        NTRACER_LOOP(ii) vj[ii] = 0.0; /* Set all tracers to 0 */
        vj[TRC+jip.current_episode] = 1.0; /* Tracer value is 1.0 for current jet episode */

        /* Loop over entire computational domain, if the jet is active! */
        TOT_LOOP(k,j,i)
        {
            JetTransformedCoordinates jtc;
            int within_r = 0, within_theta = 0;

            /* Transform variables */
            Jet_TransformToNozzleCoordinates(x1[i], x2[j], x3[k], jirp_array, &jtc);
            Jet_WithinNozzle(&jtc, &within_r, &within_theta);
            Jet_FillQuantitiesArray(&jtc, &jip, vj);

            /* -- Inject the jet if needed -- */
            if (jip.is_active == 1)
            {
#if JET_INJECTION_ZONE == CONICAL_INJECTION_ZONE
                if (within_theta == 1 && within_r == 1)
                {
                    /* -- Set the jet values -- */
                    d->Vc[RHO][k][j][i] = vj[RHO];

                    EXPAND( d->Vc[VX1][k][j][i] = vj[VX1];,
                            d->Vc[VX2][k][j][i] = vj[VX2];,
                            d->Vc[VX3][k][j][i] = vj[VX3]; )
                    d->Vc[PRS][k][j][i] = vj[PRS];

                    // Loop over tracer values
                    NTRACER_LOOP(ii) d->Vc[ii][k][j][i] = vj[ii];

                    /* -- Flag the cells we touch so that they're not updated by PLUTO -- */
                    d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
                }
                else if (smoothed_sim_region == 0 && within_theta == 0 && within_r == 1)
                {
                    /* Calculate current cylindrical radius for point */
                    R = jtc.r_offset*sin(jtc.theta);
                    /* Calculate jet equivalent cylindrical radius */
                    Cr = jtc.r_offset*cos(jtc.theta)*tan(jtc.jet_injection_properties->half_opening_angle);
                    /* Calculate smoothing factor */
                    smoothing_factor = get_smoothing(R, Cr, smoothing_power);

                    /* Get environment values */
                    environment_values[RHO] = d->Vc[RHO][k][j][i];
                    environment_values[PRS] = d->Vc[PRS][k][j][i];

                    /* HACK!! Override r_offset in jtc structure with the equivalent r_offset for the injection cone,
                     * asuming the height isn't changing */
                    jtc.r_offset = jtc.r_offset * cos(jtc.theta) / cos(jtc.jet_injection_properties->half_opening_angle);

                    /* Get jet values */
                    vj[RHO] = Jet_InjectionZoneDensity(jtc.jet_injection_properties, &jtc);
                    /* Get pressure */
                    Jet_InjectionZonePressure(jtc.jet_injection_properties, vj);

                    /* Smooth density */
                    d->Vc[RHO][k][j][i] = environment_values[RHO] + (vj[RHO] - environment_values[RHO])*smoothing_factor;
                    /* Smooth pressure */
                    d->Vc[PRS][k][j][i] = environment_values[PRS] + (vj[PRS] - environment_values[PRS])*smoothing_factor;
                }
#elif JET_INJECTION_ZONE == SPHERICAL_INJECTION_ZONE
                if (within_r) {
                    d->Vc[RHO][k][j][i] = vj[RHO];
                    d->Vc[PRS][k][j][i] = vj[PRS];
                    d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;

                    if (within_theta) {
                        NTRACER_LOOP(ii) d->Vc[ii][k][j][i] = vj[ii];
                        EXPAND( d->Vc[VX1][k][j][i] = vj[VX1];,
                                d->Vc[VX2][k][j][i] = vj[VX2];,
                                d->Vc[VX3][k][j][i] = vj[VX3]; )
                    } else {
                        NTRACER_LOOP(ii) d->Vc[ii][k][j][i] = 0.0;
                        EXPAND( d->Vc[VX1][k][j][i] = 0.0;,
                                d->Vc[VX2][k][j][i] = 0.0;,
                                d->Vc[VX3][k][j][i] = 0.0; )
                    }
                }
#endif
            }      
            if (jtc.r_from_origin < 20) {
#if SUPERNOVAE_ENERGY_INPUT == YES || STELLAR_MASS_LOSS == YES
                double stellar_density = get_stellar_density_hernquist(g_inputParam[STELLAR_MASS], x1i, x2j, x3k);
#endif
#if SUPERNOVAE_ENERGY_INPUT == YES
                double energy_inj_rate = get_energy_injection_rate(stellar_density) * pow(vn.length, 3) / vn.power;
                d->Vc[PRS][k][j][i] += energy_inj_rate * g_dt * (nonrel_gamma - 1);
#endif

#if STELLAR_MASS_LOSS == YES
                double density_loss_rate = get_mass_injection_rate(stellar_density) / vn.density * vn.time;
                d->Vc[RHO][k][j][i] += density_loss_rate * g_dt;
#endif
            }
        }
        if (jip.is_active)
            smoothed_sim_region = 1; // Set to 1 the first time we run this TOT_LOOP, so that we only smooth once
    }
}

#if (BODY_FORCE & VECTOR)
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 *********************************************************************** */
{
  g[IDIR] = g[JDIR] = g[KDIR] = 0;

  grav_accel_at_position(x1, x2, x3, g, &ep);

    /**************** 
     * Extra galaxy stuff (potential etc)!
    double extra_gal_r = calculate_offset_radius(x1, x2, x3, &secondary_ep);
    if (extra_gal_r < 50) {
        grav_accel_at_position(x1, x2, x3, g, &secondary_ep);
    }
    ****************/
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the graviational potential as function of the coordinates.
 *
 *********************************************************************** */
{
  return potential_at_position(x1, x2, x3, &ep);
}
#endif
