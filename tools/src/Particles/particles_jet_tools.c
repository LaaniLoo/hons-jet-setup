#include "pluto.h"
#include "units.h"
#include "tools_usr.h"

JetInjectionRegionProperties* pjirp[2];
JetIntermittencyProperties* pjip;

void Particles_SetJetProperties(JetInjectionRegionProperties* new_jirp1, JetInjectionRegionProperties* new_jirp2, JetIntermittencyProperties* new_jip) {
    pjirp[0] = new_jirp1;
    pjirp[1] = new_jirp2;
    pjip = new_jip;
}

int Particles_JetEnabled() {
    return pjip->is_active;
}

int Particles_WithinJetInjectionZone(double x1, double x2, double x3) {
    JetTransformedCoordinates jtc;
    int within_r = 0, within_theta = 0;

    Jet_TransformToNozzleCoordinates(x1, x2, x3, pjirp, &jtc);
    Jet_WithinNozzle(&jtc, &within_r, &within_theta);

    return (within_r && within_theta);
}

void Particles_InjectJetParticlesBatches(Data *data, Grid *grid) {
    int i, j, k, status, np, sb;
    int n0, n1;
    int np_inject_zone = RuntimeGet()->Nparticles_cell;
    double t0, t1, t_freq;
    double *x1  = grid->xgc[IDIR];
    double *x2  = grid->xgc[JDIR];
    double *x3  = grid->xgc[KDIR];

    double *dx1 = grid->dx[IDIR];
    double *dx2 = grid->dx[JDIR];
    double *dx3 = grid->dx[KDIR];

    double *xbeg = grid->xbeg;
    double *xend = grid->xend;

    Particle p;
    int pcount = 0;

    t0 = g_time;
    t1 = g_time + g_dt;
    t_freq = g_inputParam[PART_INJECT_FREQUENCY] / (SEC2MYR * vn.time);

    n0 = (int) (t0 / t_freq);
    n1 = (int) (t1 / t_freq);

    if (n0 != n1 || g_stepNumber == 0) {
        DOM_LOOP(k, j, i) {
            if (Particles_WithinJetInjectionZone(x1[i], x2[j], x3[k])) {
                for (np = 0; np < np_inject_zone; np++) {

                    /* Initialise our particles */
                    p.coord[0] = x1[i];
                    p.coord[1] = x2[j];
                    p.coord[2] = x3[k];
#if PARTICLES_RANDOM_MODE == PARTICLES_RANDOM_GAUSSIAN
                    EXPAND(p.coord[0] = x1[i] + GaussianRandomNumber(0.0, dx1[i] * PARTICLES_RANDOM_GAUSSIAN_FACTOR);,
                           p.coord[1] = x2[j] + GaussianRandomNumber(0.0, dx2[j] * PARTICLES_RANDOM_GAUSSIAN_FACTOR);,
                           p.coord[2] = x3[k] + GaussianRandomNumber(0.0, dx3[k] * PARTICLES_RANDOM_GAUSSIAN_FACTOR););

                    /* Clamp coords to grid */
                    p.coord[0] = clamp(p.coord[0], x1[grid->lbeg[IDIR]], x1[grid->lend[IDIR]]);
                    p.coord[1] = clamp(p.coord[1], x2[grid->lbeg[JDIR]], x2[grid->lend[JDIR]]);
                    p.coord[2] = clamp(p.coord[2], x3[grid->lbeg[KDIR]], x3[grid->lend[KDIR]]);
#endif
                    p.color = 0.0;
                    p.speed[IDIR] = 0.0;
                    p.speed[JDIR] = 0.0;
                    p.speed[KDIR] = 0.0;

#if (PARTICLES_TYPE == LAGRANGIAN) && (PARTICLES_LP_RAISE == YES)
                    p.density = 0.0;
                    p.pressure = 0.0;
                    for (sb=0; sb < PARTICLES_LP_SHK_BINS; sb++)
                        p.last_shock_time[sb] = 0.0;

                    for (sb=0; sb < NTRACER; sb++)
                        p.tracer[sb] = 0.0;
#endif

                    status = Particles_Insert(&p, data, PARTICLES_CREATE, grid);

                    if (status != TRUE) {
                        Particles_Display(&p);
                        print("! Particles_Inject(): error inserting particle\n");
                        QUIT_PLUTO(1);
                    }
                    pcount++;
                }
            }
        }
    }
}

void Particles_InjectJetParticlesStream(Data *data, Grid *grid) {
    int i, j, k, status, sb;
    int ii;
    int nd;
    double t0, t1;
    int n0, n1;
    int particles_to_inject;
    int skip_particle = 0;
    double r[2], theta[2], phi[2];
    Particle p;
    double gx1, gx2, gx3;
    double gx[3];

    int np_inject_zone = RuntimeGet()->Nparticles_cell;
    double delta_t = g_inputParam[PART_INJECT_FREQUENCY] / (SEC2MYR * vn.time);
    t0 = g_time;
    t1 = g_time + g_dt;

    n0 = (int) (t0 / delta_t * np_inject_zone * 0.5);
    n1 = (int) (t1 / delta_t * np_inject_zone * 0.5);

    particles_to_inject = g_stepNumber == 0 ? np_inject_zone : (n1 - n0);

    if (n0 != n1 || g_stepNumber == 0) {
        while (particles_to_inject > 0) {
            /* Calculate random r, theta, phi */
            if (prank == 0) {
                for (ii = 0; ii < 2; ii++) {
#if JET_INJECTION_ZONE == CONICAL_INJECTION_ZONE
                    /* For conical injection zone, radius can be anything up to injection */
                    r[ii] = RandomNumber(0, pjirp[ii]->injection_region_radius);
#elif JET_INJECTION_ZONE == SPHERICAL_INJECTION_ZONE
                    /* For spherical injection zone, need to prevent particles getting stuck in the 'deadzone' around
                     * r = 0. This depends on the resolution & opening angle, so to be safe we let the radius vary between
                     * 0.9r0 and r0 */
                    r[ii] = RandomNumber(pjirp[ii]->injection_region_radius*0.90, pjirp[ii]->injection_region_radius);
#endif
                    theta[ii] = RandomNumber(ii * (CONST_PI - pjirp[ii]->half_opening_angle), pjirp[ii]->half_opening_angle + ii * (CONST_PI - pjirp[ii]->half_opening_angle));
                    phi[ii] = RandomNumber(0, 2.*CONST_PI);
                }
            }

#ifdef PARALLEL
            MPI_Bcast(r, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(theta, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(phi, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

            for (ii = 0; ii < 2; ii++) {
                Jet_NozzleToGrid(pjirp[ii], r[ii], theta[ii], phi[ii], &gx1, &gx2, &gx3);
                gx[0] = gx1;
                gx[1] = gx2;
                gx[2] = gx3;

                skip_particle = 0;
                DIM_LOOP(nd) {  
                  if (   gx[nd] <  grid->xbeg[nd]
                      || gx[nd] >= grid->xend[nd]) {
                      skip_particle = 1;
                  }
                }

                if (skip_particle)
                    continue;

                /* Initialise our particles */
                p.coord[0] = gx1;
                p.coord[1] = gx2;
                p.coord[2] = gx3;
                p.color = 0.0;
                p.speed[IDIR] = 0.0;
                p.speed[JDIR] = 0.0;
                p.speed[KDIR] = 0.0;

#if (PARTICLES_TYPE == LAGRANGIAN) && (PARTICLES_LP_RAISE == YES)
                p.density = 0.0;
                p.pressure = 0.0;
                for (sb=0; sb < PARTICLES_LP_SHK_BINS; sb++)
                    p.last_shock_time[sb] = 0.0;

                for (sb=0; sb < NTRACER; sb++)
                    p.tracer[sb] = 0.0;
#endif

                status = Particles_Insert(&p, data, PARTICLES_CREATE, grid);

                if (status != TRUE) {
                    Particles_Display(&p);
                    print("! Particles_Inject(): error inserting particle\n");
                    QUIT_PLUTO(1);
                }
            }
            particles_to_inject -= 2;
        }
    }
}
