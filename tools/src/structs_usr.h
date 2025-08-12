#ifndef structs_usr_h
#define structs_usr_h

/** Identify the jet and counter jet
 */
typedef enum {MAIN_JET, COUNTER_JET} jet_identifier;

/** Jet injection region properties.
 * This structure holds the injection region properties of a jet
 */
typedef struct {
    double injection_region_radius; ///< The radius of the injection region
    double initial_width; ///< The desired initial with of the injection region
    double injection_cone_offset; ///< The calculated offset of the injection cone
    double half_opening_angle; ///< The half opening angle of the injection region
    double injection_cone_height; ///< The height of the injection region, if no cap is used
    double injection_rotation_angle; ///< The rotation angle of the injection region
    double solid_angle;
    double L1; ///< Inner length scale
    double L1b; ///< Radius where jet density and environment density are equal
    double density; ///< Initial density of jet in the injection region
    double pressure; ///< The pressure of the jet in the injection region
    double enthalpy; ///< The specific enthalpy of the jet in the injection region
    double lorentz; ///< The lorentz factor of the jet in the injection region
    double chi; ///< Ratio of rest-mass energy to enthalpy
    double density_contrast; ///< Ratio of proper jet density to environment density
    double enthalpy_contrast; ///< Ratio of jet enthalpy to environment enthalpy
    double energy_density_constrast; ///< Ratio of jet energy density to environment energy density
    double area; ///< Cross-sectional area for flux calculation
    double power; ///< The power of the jet
    double speed; ///< Speed of the jet
    int collimated; ///< Whether the jet is initially collimated or not  (MAY NOT WORK, NEED TO TEST)
    double initial_offset[3]; ///< The initial jet offset
    double inlet_speed[3]; ///< The jet inlet speed
    double offset[3]; ///< The current jet offset, recalculated each timestep to account for moving injection region
    double pressure_factor; ///< How much the jet is overpressured compared to the environment at the inlet
} JetInjectionRegionProperties;

/** Jet intermittency properties.
 * This structure holds the simulation properties for jet intermittency
 */
typedef struct {
    double delta_t; ///< The ramp-up/down timescale for the jet intermittency
    int episodes; ///< The number of jet episdoes
    double start_time; //< When this jet starts
    double end_time; ///< When this jet ends
    double total_active_time; ///< The total active time of the jet
    double active_time; ///< How long the jet is active for per outburst
    double quiet_time; ///< How long the jet is quiescent for per outburst
    double duty_cycle_time; ///< The total time to complete one duty cycle for the jet

    int is_active; ///< Is the jet currently active?
    int current_episode; ///< Current episode number
    int have_set_timestep; ///< Whether the timestep has been reset
    double intermittency_velocity_factor; ///< Velocity factor used for ramping jet up & down
} JetIntermittencyProperties;

typedef struct {
    double r;
    double r_from_origin;
    double r_offset;
    double theta;
    double phi;
    jet_identifier jet_id;
    JetInjectionRegionProperties* jet_injection_properties;
} JetTransformedCoordinates;

typedef struct {
    double rho_0; ///< Central density
    double r_scaling; ///< King profile core radii or Makino scale radius
    double b_exponent; ///< King profile beta exponent
    double delta_nfw; ///< NFW delta parameter
    double central_temperature; ///< Central temperature (isothermal)
    double central_sound_speed; ///< Central sound speed
    double offset[3]; ///< Environment center offset
} AnalyticEnvironmentProperties;

typedef struct {
    /* Variables to hold the input data ids */
    int did;
    int pid;
    int vid_x;
    int vid_y;
    int vid_z;
    int gid_x;
    int gid_y;
    int gid_z;

    /* Variables to hold the input data variable locations */
    int density_pos;
    int pressure_pos;
    int vx_pos;
    int vy_pos;
    int vz_pos;
    int gx_pos;
    int gy_pos;
    int gz_pos;

    char external_env_data[300]; ///< External data file name
    char external_env_out[300]; ///< External grid.out file name

    double offset[3]; ///< Environment center offset

    double environment_velocity_factor; ///< Environment velocities are multiplied by this factor
} ExternalEnvironmentProperties;

typedef struct {
    double *hot_rad; ///< Tabulated radius array
    double *hot_rho; ///< Tabulated density array
    double *hot_prs; ///< Tabulated pressure array
    double *hot_T; ///< Tabulated temperature array
    double *hot_diffPotMass; ///< Tabulated dphi_dr (using -GM(<r)/r^2)
    double *hot_diffPot; ///< Tabulated dphi_dr (using (1/rho)*(dp/dr))
    int hot_ndata; ///< Number of entries in our tabulated dataset

    char file_name[300]; ///< Tabulated environment file name
    double offset[3]; ///< Environment center offset
} TabulatedEnvironmentProperties;

#endif
