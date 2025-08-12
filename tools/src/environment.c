#include "environment.h"
#include "definitions_usr.h"

void* default_environment = NULL;

inline void setup_analytic_environment(AnalyticEnvironmentProperties* ep) {
    ep->rho_0 = g_inputParam[ENV_RHO_0] / vn.density;
    ep->r_scaling = g_inputParam[ENV_R_SCALING];
    ep->central_temperature = g_inputParam[ENV_TEMP];
#if DENSITY_PROFILE == ENV_MAKINO
    ep->delta_nfw = g_inputParam[ENV_DELTA_NFW];
#elif DENSITY_PROFILE == ENV_KING
    ep->b_exponent = g_inputParam[ENV_B_EXPONENT];
#endif
    ep->central_sound_speed = sqrt((nonrel_gamma * ep->central_temperature) / (MU_NORM * KELVIN));

      ep->offset[0] = g_inputParam[ENV_X1O];
      ep->offset[1] = g_inputParam[ENV_X2O];
      ep->offset[2] = g_inputParam[ENV_X3O];


    print("> Loaded analytic environment properties");
}

inline void setup_external_environment(ExternalEnvironmentProperties* ep) {

    /* Initiliase input data IDs to -1 */
    ep->did = ep->pid = ep->gid_x = ep->gid_y = ep->gid_z = ep->vid_x = ep->vid_y = ep->vid_z = -1;

    /* Set up default variable numbering */
    ep->density_pos = 0;
    ep->pressure_pos = 1;
    ep->gx_pos = 2;
    ep->gy_pos = 3;
    ep->gz_pos = 4;
    ep->vx_pos = 5;
    ep->vy_pos = 6;
    ep->vz_pos = 7;

    /* Set up environment offset
     * NOTE: Currently UNUSED for external environments */
    ep->offset[0] = g_inputParam[ENV_X1O];
    ep->offset[1] = g_inputParam[ENV_X2O];
    ep->offset[2] = g_inputParam[ENV_X3O];

    /* Set environment velocity factor */
    ep->environment_velocity_factor = EXT_ENV_VELOCITY_FACTOR;

    print("> Loaded external environment properties");
}

void PrintEnvironmentProperties(void* bep) {
#if DENSITY_PROFILE == ENV_MAKINO || DENSITY_PROFILE == ENV_KING
    AnalyticEnvironmentProperties* ep = (AnalyticEnvironmentProperties*) bep;
    print("\n\nEnvironment properties:\n\n");
    print("rho_0               = %16e [%16e g/cm^3] \n", ep->rho_0, ep->rho_0 * vn.density);
    print("r_scaling           = %16e [%16e kpc] \n", ep->r_scaling, ep->r_scaling * vn.length * CM2KPC);
    print("central_temperature = %16e [%16e K] \n", ep->central_temperature, ep->central_temperature * vn.temperature);
    print("delta_nfw           = %16e \n", ep->delta_nfw);
    print("b_exponent          = %16e \n", ep->b_exponent);
    print("central_pressure      = %16e [%16e barye] \n", pressure_at_position(0, 0, 0, bep), pressure_at_position(0, 0, 0, bep) * vn.pressure);
    print("central_sound_speed = %16e [%16e km/s] \n", ep->central_sound_speed, ep->central_sound_speed * vn.velocity * 1e-5);
    print("X1 offset           = %16e [%16e kpc] \n", ep->offset[0], ep->offset[0] * vn.length * CM2KPC);
    print("X2 offset           = %16e [%16e kpc] \n", ep->offset[1], ep->offset[1] * vn.length * CM2KPC);
    print("X3 offset           = %16e [%16e kpc] \n", ep->offset[2], ep->offset[2] * vn.length * CM2KPC);
    print("\n\n");
#elif DENSITY_PROFILE == ENV_EXTERNAL
    ExternalEnvironmentProperties* ep = (ExternalEnvironmentProperties*) bep;

    print("\n\nLoading external environment:\n\n");
    
    /* Force loading of density, pressure, and velocity */
    density_at_position(0, 0, 0, bep);
    pressure_at_position(0, 0, 0, bep);
    double v_tmp[3];
    velocity_at_position(0, 0, 0, v_tmp, bep);
    grav_accel_at_position(0, 0, 0, v_tmp, bep);

    print("\n\nEnvironment properties:\n\n");
    print("data file            = %s\n", ep->external_env_data);
    print("out file             = %s\n", ep->external_env_out);
    print("central_density      = %16e [%16e g/cm^3] \n", density_at_position(0, 0, 0, bep), density_at_position(0, 0, 0, bep) * vn.density);
    print("central_pressure      = %16e [%16e barye] \n", pressure_at_position(0, 0, 0, bep), pressure_at_position(0, 0, 0, bep) * vn.pressure);
    print("central_temperature  = %16e [%16e K] \n", temperature_at_position(0, 0, 0, bep), temperature_at_position(0, 0, 0, bep) * vn.temperature);
    print("central_sound_speed  = %16e [%16e km/s] \n", get_sound_speed(temperature_at_position(0, 0, 0, bep)), get_sound_speed(temperature_at_position(0, 0, 0, bep)) * vn.velocity * 1e-5);
    print("X1 offset            = %16e [%16e kpc] \n", ep->offset[0], ep->offset[0] * vn.length * CM2KPC);
    print("X2 offset            = %16e [%16e kpc] \n", ep->offset[1], ep->offset[1] * vn.length * CM2KPC);
    print("X3 offset            = %16e [%16e kpc] \n", ep->offset[2], ep->offset[2] * vn.length * CM2KPC);
    print("\n\n");
#elif DENSITY_PROFILE == ENV_TABULATED
    TabulatedEnvironmentProperties* ep = (TabulatedEnvironmentProperties*) bep;
    print("\n\nEnvironment properties:\n\n");
    print("data file            = %s\n", ep->file_name);
    print("central_density      = %16e [%16e g/cm^3] \n", density_at_position(0, 0, 0, bep), density_at_position(0, 0, 0, bep) * vn.density);
    print("central_pressure      = %16e [%16e barye] \n", pressure_at_position(0, 0, 0, bep), pressure_at_position(0, 0, 0, bep) * vn.pressure);
    print("central_temperature  = %16e [%16e K] \n", temperature_at_position(0, 0, 0, bep), temperature_at_position(0, 0, 0, bep) * vn.temperature);
    print("central_sound_speed  = %16e [%16e km/s] \n", get_sound_speed(temperature_at_position(0, 0, 0, bep)), get_sound_speed(temperature_at_position(0, 0, 0, bep)) * vn.velocity * 1e-5);
    print("X1 offset            = %16e [%16e kpc] \n", ep->offset[0], ep->offset[0] * vn.length * CM2KPC);
    print("X2 offset            = %16e [%16e kpc] \n", ep->offset[1], ep->offset[1] * vn.length * CM2KPC);
    print("X3 offset            = %16e [%16e kpc] \n", ep->offset[2], ep->offset[2] * vn.length * CM2KPC);
    print("\n\n");
#endif
}

void set_default_environment(void* bep)
{
    default_environment = bep;
}

inline double density_at_position(const double x1, const double x2, const double x3, void* bep)
{
    /* First of all we set our bep pointer if it is null */
    if (bep == NULL)
        bep = default_environment;

    /* We have several different environment cases to handle, depending on what the DENSITY_PROFILE preprocessor
     * variable is set to.  First we check for the external environment case */
#if DENSITY_PROFILE == ENV_EXTERNAL

    /* Get our ExternalEnvironmentProperties structure */
    ExternalEnvironmentProperties* eep = (ExternalEnvironmentProperties*) bep;

  if (eep->did == -1) {
      /* If our external input data is not loaded, we load it now */
      print("[InputData] density_at_position\n"); // debugging
      eep->did = InputDataOpen(eep->external_env_data, eep->external_env_out," ",eep->density_pos);
  }

  /* Interpolate the density and return it */
  double density = InputDataInterpolate(eep->did,x1,x2,x3)/vn.density;
  return density;

  /* Now we check for the tabulated case */
#elif DENSITY_PROFILE == ENV_TABULATED
    HotHaloPrimitives(halo_primitives, x1, x2, x3, (TabulatedEnvironmentProperties*)bep);
    return halo_primitives[RHO];
#endif

    /* Next we check for either the makino or king profiles, which need the current spherical radius, and then go on to calculate the density for these environments */
#if DENSITY_PROFILE == ENV_MAKINO || DENSITY_PROFILE == ENV_KING
    AnalyticEnvironmentProperties* ep = (AnalyticEnvironmentProperties*) bep;
    double r = calculate_offset_radius(x1, x2, x3, bep);
#if DENSITY_PROFILE == ENV_MAKINO
    /* Makino Profile (Makino+1998) */
    if (r < 0.001) return ep->rho_0 * exp(-ep->delta_nfw);
    return ep->rho_0 * exp(-ep->delta_nfw)*pow((1+(r/ep->r_scaling)), ep->delta_nfw/(r/ep->r_scaling));
#elif DENSITY_PROFILE == ENV_KING
    /* King Profile (from Krause 2005) */
    return ep->rho_0 * pow((1+(r/ep->r_scaling)*(r/ep->r_scaling)),-1.5*ep->b_exponent); 
#endif
#endif
}

inline double temperature_at_position(const double x1, const double x2, const double x3, void* bep) {

    /* First of all we set our bep pointer if it is null */
    if (bep == NULL)
        bep = default_environment;

    /* Return the environment temperature, depending on which environment is being used */
#if DENSITY_PROFILE == ENV_EXTERNAL || DENSITY_PROFILE == ENV_TABULATED
    /* If environment is tabulated or external, we calculate temperature from density and pressure */
    return get_temperature(density_at_position(x1, x2, x3, bep), pressure_at_position(x1, x2, x3, bep));
#else
    /* Otherwise we have an isothermal environment, and just return the central temperature */
    return ((AnalyticEnvironmentProperties*) bep)->central_temperature;
#endif
}

inline double pressure_at_position(const double x1, const double x2, const double x3, void* bep) {

    /* First of all we set our bep pointer if it is null */
    if (bep == NULL)
        bep = default_environment;

    /* We have several different environment cases to handle, depending on what the DENSITY_PROFILE preprocessor
     * variable is set to.  First we check for the external environment case */
#if DENSITY_PROFILE == ENV_EXTERNAL

    /* Get our ExternalEnvironmentProperties structure */
    ExternalEnvironmentProperties* eep = (ExternalEnvironmentProperties*) bep;

  if (eep->pid == -1) {
      /* If our external input data is not loaded, we load it now */
      print("[InputData] pressure_at_position\n"); // debugging
      eep->pid = InputDataOpen(eep->external_env_data, eep->external_env_out," ",eep->pressure_pos);
  }

  /* Interpolate the pressure and return it */
  double pressure = InputDataInterpolate(eep->pid,x1,x2,x3)/vn.pressure;
  return pressure;
#elif DENSITY_PROFILE == ENV_TABULATED
    /* Now we check for the tabulated case */
    HotHaloPrimitives(halo_primitives, x1, x2, x3, (TabulatedEnvironmentProperties*)bep);
    return halo_primitives[PRS];
#else
    /* Otherwise we must have either the King or Makino profiles - these are by nature isothermal, so we calculate
     * the pressure from the current density and (central) environment temperature */
    return get_pressure(density_at_position(x1, x2, x3, (AnalyticEnvironmentProperties*)bep), temperature_at_position(x1, x2, x3, (AnalyticEnvironmentProperties*)bep));
#endif
}

inline void grav_accel_at_position(const double x1, const double x2, const double x3, double *g, void* bep) {
    /* First of all we set our bep pointer if it is null */
    if (bep == NULL)
        bep = default_environment;

#if DENSITY_PROFILE == ENV_MAKINO || DENSITY_PROFILE == ENV_KING
    AnalyticEnvironmentProperties* ep = (AnalyticEnvironmentProperties*) bep;
    double *offsets = ep->offset;
#endif

    /* The gravitational acceleration also depends on which environment we are using */
#if DENSITY_PROFILE == ENV_EXTERNAL
    /* First we check for the external case */

    /* Get our ExternalEnvironmentProperties structure */
    ExternalEnvironmentProperties* eep = (ExternalEnvironmentProperties*) bep;

    /* If our external input data is not loaded, we load it now */
  if (eep->gid_x == -1) {
    print("[InputData] BodyForceVector_x\n"); // debugging
    eep->gid_x = InputDataOpen(eep->external_env_data, eep->external_env_out," ",eep->gx_pos);
  }
  if (eep->gid_y == -1) {
    print("[InputData] BodyForceVector_y\n"); // debugging
    eep->gid_y = InputDataOpen(eep->external_env_data, eep->external_env_out," ",eep->gy_pos);
  }
  if (eep->gid_z == -1) {
    print("[InputData] BodyForceVector_z\n"); // debugging
    eep->gid_z = InputDataOpen(eep->external_env_data, eep->external_env_out," ",eep->gz_pos);
  }

  g[IDIR] = InputDataInterpolate(eep->gid_x, x1, x2, x3)/vn.acceleration;
  g[JDIR] = InputDataInterpolate(eep->gid_y, x1, x2, x3)/vn.acceleration;
  g[KDIR] = InputDataInterpolate(eep->gid_z, x1, x2, x3)/vn.acceleration;
  return;
#else

    /* Now we handle the radial acceleration profiles. */
    double radial_accel = 0.0;
    double r_prime = calculate_offset_radius(x1, x2, x3, bep);

    /* Let's calculate our radial acceleration and environment offsets */
#if DENSITY_PROFILE == ENV_TABULATED
  /* Next we check for the tabulated case */
    radial_accel = GravitationalAcceleration(x1, x2, x3, (TabulatedEnvironmentProperties*)bep);
    double *offsets = ((TabulatedEnvironmentProperties*)bep)->offset;
#elif DENSITY_PROFILE == ENV_MAKINO
    /* Handle the Makino case */
    double phi_0 = -ep->delta_nfw * (pow(ep->central_sound_speed, 2) / nonrel_gamma);
    radial_accel = -phi_0 * ep->r_scaling * (r_prime - (ep->r_scaling + r_prime)*log(1 + (r_prime/ep->r_scaling))) / (r_prime * r_prime * (ep->r_scaling + r_prime));
#elif DENSITY_PROFILE == ENV_KING
    /* Finally we handle the King profile case */
    radial_accel = - 3.0 * ep->b_exponent * (pow(ep->central_sound_speed, 2) / nonrel_gamma) * r_prime / (pow(ep->r_scaling, 2) + pow(r_prime,2));
#endif

    /* Now let's transform coordinate systems back so that they are centered on the grid origin */
#if COMPONENTS == 2 && GEOMETRY == SPHERICAL
    /* Our gravitational potential is a vector with the orign at the centre of the environment.
     * If our environment is offset, then this does not correspond to the origin of our simulation grid,
     * so we need to translate this so that it is a vector with the origin at the centre of our grid.
     * To do this, we note that our basis vector for r_prime can be written as
     *
     * r_prime_hat = (1/r_prime) * (r - a cos(theta)) r_hat + (1/r_prime) * (a sin(theta)) theta_hat
     *
     * where r_prime is our radius with respect to environment centre, r is radius with respect to grid origin,
     * a is the environment offset, and theta is our angle with respect to the grid origin.
     */
    g[IDIR] += radial_accel * ((x1 - offsets[0] * cos(x2)) / r_prime);
    g[JDIR] += radial_accel * (offsets[0] * sin(x2) / r_prime);
#elif GEOMETRY == CARTESIAN
    g[IDIR] += radial_accel * (x1 - offsets[0]) / r_prime;
    g[JDIR] += radial_accel * (x2 - offsets[1]) / r_prime;
    g[KDIR] += radial_accel * (x3 - offsets[2]) / r_prime;
#endif
#endif
}

inline double potential_at_position(const double x1, const double x2, const double x3, void* bep)
{
    /* The gravitatational potential also depends on which environment we are using.
     * Currently only the Makino and King profiles support using the gravitational potential. */
#if DENSITY_PROFILE == ENV_MAKINO || DENSITY_PROFILE == ENV_KING
    AnalyticEnvironmentProperties* ep = (AnalyticEnvironmentProperties*) bep;
    double r = calculate_offset_radius(x1, x2, x3, bep);
#endif

    /* Calculate the gravitational potential for the Makino profile */
#if DENSITY_PROFILE == ENV_MAKINO
    /* -- Makino profile coefficient -- */
    double phi_0 = -ep->delta_nfw * (pow(ep->central_sound_speed, 2) / nonrel_gamma);

    /* -- Body force potential is found from the Makino gas density profile -- */
    return phi_0 * log(1 + (r / ep->r_scaling)) / (r / ep->r_scaling);
    /* if (r <= 0.0) return 0; */

    /* Calculate the gravitational potential for the King profile */
#elif DENSITY_PROFILE == ENV_KING
    /* King body force potential, see Krause 2005 */
    double phi_0 = 1.5*ep->b_exponent*(pow(ep->central_sound_speed, 2) / nonrel_gamma);
    return phi_0*log(1+(r/ep->r_scaling)*(r/ep->r_scaling));
#elif DENSITY_PROFILE == ENV_TABULATED
    print("Gravitational potential is not supported with tabulated profiles\n");
    QUIT_PLUTO(1);
#elif DENSITY_PROFILE == ENV_EXTERNAL
    print("Gravitational potential is not supported with external profiles\n");
    QUIT_PLUTO(1);
#endif
}

inline void velocity_at_position(const double x1, const double x2, const double x3, double *v, void* bep) {
    /* The initial velocity for a given position is only relevant if we have an external environment */

    /* First of all we set our bep pointer if it is null */
    if (bep == NULL)
        bep = default_environment;
#if (DENSITY_PROFILE == ENV_EXTERNAL) && (EXT_ENV_VELOCITY == YES)
    /* Get our ExternalEnvironmentProperties structure */
    ExternalEnvironmentProperties* eep = (ExternalEnvironmentProperties*) bep;

    /* If our external input data is not loaded, we load it now */
  if (eep->vid_x == -1) {
    print("[InputData] VelocityVector_x\n"); // debugging
    eep->vid_x = InputDataOpen(eep->external_env_data, eep->external_env_out," ",eep->vx_pos);
  }
  if (eep->vid_y == -1) {
    print("[InputData] VelocityVector_y\n"); // debugging
    eep->vid_y = InputDataOpen(eep->external_env_data, eep->external_env_out," ",eep->vy_pos);
  }
  if (eep->vid_z == -1) {
    print("[InputData] VelocityVector_z\n"); // debugging
    eep->vid_z = InputDataOpen(eep->external_env_data, eep->external_env_out," ",eep->vz_pos);
  }
   v[IDIR] = (InputDataInterpolate(eep->vid_x, x1, x2, x3)/vn.velocity)*eep->environment_velocity_factor;
   v[JDIR] = (InputDataInterpolate(eep->vid_y, x1, x2, x3)/vn.velocity)*eep->environment_velocity_factor;
   v[KDIR] = (InputDataInterpolate(eep->vid_z, x1, x2, x3)/vn.velocity)*eep->environment_velocity_factor;
  return;
#else
    v[IDIR] = 0.0;
    v[JDIR] = 0.0;
    v[KDIR] = 0.0;
#endif
}

/* NOTE: The following 3 helper functions assume an IDEAL equation of state */

/* Helper function to caculate pressure from density and temperature */
inline double get_pressure(double rho, double T)
{
    return (rho * T) / (MU_NORM * KELVIN);
}

/* Helper function to calculate temperature from density and pressure */
inline double get_temperature(double rho, double P)
{
    return (P * MU_NORM * KELVIN) / (rho);
}

/* Helper function to calculate density from pressure and temperature */
inline double get_density(double P, double T)
{
    return (P * MU_NORM * KELVIN) / (T);
}

/* Helper function to calculate the sound speed for a given temperature */
inline double get_sound_speed(double T)
{
    return sqrt((nonrel_gamma * T) / (MU_NORM * KELVIN));
}

/* Gravitational force terms */
inline double nfw_gravitational_acceleration(const double x1, const double x2, const double x3) {
    /* -- Makino profile coefficient -- */
    double c_0_klypin_planck_relaxed = 7.75;
    double gamma_klypin_planck_relaxed = 0.100;
    double m_0_klypin_planck_relaxed = 4.5e5;

    double virial_mass = 1e13;

    double delta_vir = 200;

    double hubble = 0.7;

    double concentration = (c_0_klypin_planck_relaxed * (pow(virial_mass / (1e12 * (pow(hubble,(-1)))), -gamma_klypin_planck_relaxed)) * (1 + pow(virial_mass / (m_0_klypin_planck_relaxed * (1e12 * (pow(hubble,(-1))))),0.4)));

    double char_density = (delta_vir / 3.0) * (pow(concentration,3.0) / (log(1.0 + concentration) - (concentration / (1.0 + concentration))));

    double redshift = 0;
    double H_nought = 67.3e-19 / 3.086;
    double Om = 0.316;
    double Ode = 1 - Om;
    double H = H_nought * sqrt((Om * pow(1+redshift, 3) + Ode));
    double critical_density = ((3.0 * H) / (8.0 * CONST_PI * CONST_G)) * H;

    double virial_radius = cbrt(virial_mass / (delta_vir * critical_density * (4.0/3.0) * CONST_PI)) * cbrt(CONST_Msun);

    double rs = virial_radius / concentration;

    double phi_0 = - 4 * CONST_PI * CONST_G * char_density * critical_density * pow(rs, 2);

    double radius = sph_radius(x1, x2, x3);
    if (radius < 0.001)
        radius = 0.001;
    radius *= vn.length;

    double dim_r = radius / rs;

    return -(phi_0 / rs) * pow(dim_r, -2) * ((dim_r/(1+dim_r)) - log(1 + dim_r)) / vn.acceleration;
}

inline double hernquist_gravitational_acceleration(const double stellar_mass, const double x1, const double x2, const double x3) {

    double r_eff = get_r_eff_elliptical(stellar_mass); /* This is in kpc */
    double a_scale=(r_eff * vn.length)/1.815; /* Convert r_eff to cm, use this to calculate the scale radius */
    double radius = sph_radius(x1, x2, x3) * vn.length;

    return -((CONST_G * stellar_mass * CONST_Msun) / pow(radius + a_scale, 2)) / vn.acceleration;
}

void set_external_env_paths(ExternalEnvironmentProperties* ep, const char* default_data, const char* default_grid) {
    FILE * fp;

    if ((fp = fopen(g_cmdLine.env_file, "r")) != NULL) {

        /* If the environment.ini file exists, load the ExternalData & ExternalGrid keys from that */
        fclose(fp);

        ParamFileRead(g_cmdLine.env_file);
        sprintf(ep->external_env_data, "%s", ParamFileGet("ExternalData", 1));
        sprintf(ep->external_env_out, "%s", ParamFileGet("ExternalGrid", 1));
    }
    else {
        /* Else, use the supplied default values */
        strcpy(ep->external_env_data, default_data);
        strcpy(ep->external_env_out, default_grid);
    }
}

void set_tabulated_env_path(TabulatedEnvironmentProperties* ep, const char* default_data) {
    FILE * fp;

    if ((fp = fopen(g_cmdLine.env_file, "r")) != NULL) {

        /* If the environment.ini file exists, load the TabulatedData key from that */
        fclose(fp);

        ParamFileRead(g_cmdLine.env_file);
        sprintf(ep->file_name, "%s", ParamFileGet("TabulatedData", 1));
    }
    else {
        /* Else, use the supplied default values */
        strcpy(ep->file_name, default_data);
    }
}
