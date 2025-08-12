#include "read_external_pressure.h"

void ReadHotAndGravTable(TabulatedEnvironmentProperties *tep) {
    /*
     * This routine reads the data from a hot phase data file. Gas is assumed to be in hydro equilibrium.
     *
     * The data should be a six-column file in cgs units.
     * Columns are: radius, density, pressure, temperature, -GM(<r)/r^2, (1/rho)*(dp/dr). The last two are alternative formulations of dPhi/dr in HE.
     *
     * The values of these are filled into the global arrays hot_rad, hot_rho, hot_prs, and hot_dphi_dr, which
     * are used throughout the code.
     *
     * Returns 0
     *
     * */

    FILE *f;

    double buf;
    int i;

    /* Open file */
    if ((f = fopen(tep->file_name, "r")) == NULL) {
        print("Error: ReadHotData: Unable to open file\n");
        print("%s", tep->file_name);
        QUIT_PLUTO(1);
    }

    /* Scan file first to get number of lines*/
    tep->hot_ndata = 0;
    while (fscanf(f, "%le %le %le %le %le %le", &buf, &buf, &buf, &buf, &buf, &buf) != EOF) {
        tep->hot_ndata++;
    }
    
    /* Allocate memory for profile arrays */
    tep->hot_rad = ARRAY_1D(tep->hot_ndata, double);
    tep->hot_rho = ARRAY_1D(tep->hot_ndata, double);
    tep->hot_prs = ARRAY_1D(tep->hot_ndata, double);
    tep->hot_T = ARRAY_1D(tep->hot_ndata, double);
    tep->hot_diffPotMass = ARRAY_1D(tep->hot_ndata, double);
    tep->hot_diffPot = ARRAY_1D(tep->hot_ndata, double);

    /* Read data */
    fseek(f, 0, SEEK_SET);
    for (i = 0; i < tep->hot_ndata; ++i) {
        fscanf(f, "%le ", &tep->hot_rad[i]);
        fscanf(f, "%le ", &tep->hot_rho[i]);
        fscanf(f, "%le ", &tep->hot_prs[i]);
        fscanf(f, "%le ", &tep->hot_T[i]);
        fscanf(f, "%le ", &tep->hot_diffPotMass[i]);
        fscanf(f, "%le ", &tep->hot_diffPot[i]);
    }

    /* Clean up */
    fclose(f);

    /* Convert variables into code units */
    for (i = 0; i < tep->hot_ndata; ++i) {
        tep->hot_rad[i] /= vn.length;
        tep->hot_rho[i] /= vn.density;
        tep->hot_prs[i] /= vn.pressure;
        tep->hot_diffPot[i] /= vn.acceleration;

	//	hot_diffPot[i] /= ((UNIT_VELOCITY*UNIT_VELOCITY)/UNIT_LENGTH); /* these are units of p / (rho*r) */
	// something VERY STRANGE is going on with calculation unit of grav acceleration: c^2/kpc in cgs should be ~0.3, but get 3e5. Hence scale down grav pot too much, and it is never strong enough to stop outflows.
	/* hot_diffPot[i] /= 0.2916; /1* do this *by hand*!!! *1/ */
    }
    
}

double GravitationalAcceleration(const double x1, const double x2, const double x3, TabulatedEnvironmentProperties *tep) {
    if (tep->hot_rad == NULL) {
        ReadHotAndGravTable(tep);
    }
    return InterpolationWrapper(tep->hot_rad, tep->hot_diffPot, tep->hot_ndata, calculate_offset_radius(x1, x2, x3, tep));
}


/* ************************************************************** */
void HotHaloPrimitives(double *halo,
                       const double x1, const double x2, const double x3, TabulatedEnvironmentProperties *tep) {
/*
 * Return array of primitives containing Halo quantities
 *
 * double    halo         array of halo primitives
 * double    x1, x2, x3   first, second, third coordinate
 *
 **************************************************************** */

    double vel, scrh;
    double inv_unit_G;
    double frac;
    double y0, y1, y2, y3, r1, r2;
    int iv, il, ic;
    static int once01 = 0;

    if (tep->hot_rad == NULL) {
        ReadHotAndGravTable(tep);
    }

    /* The density from the table is in units of g cm^-3 */
    /* double r = sph_radius(x1 + tep->offset[0], x2 + tep->offset[1], x3 + tep->offset[2]); */

    double r = calculate_offset_radius(x1, x2, x3, tep);
    
    /* r = D_EXPAND(x1*x1, + x2*x2, + x3*x3); */
    /* r = sqrt(r); */
    
    halo[RHO] = InterpolationWrapper(tep->hot_rad, tep->hot_rho, tep->hot_ndata, r);
    halo[PRS] = InterpolationWrapper(tep->hot_rad, tep->hot_prs, tep->hot_ndata, r);
    
    D_EXPAND(halo[VX1] = 0;,
             halo[VX2] = 0;,
             halo[VX3] = 0; )

    /* A message for having initialized halo with potential*/
    if (!once01) {
        print("> Initializing hot halo distribution from table");
        once01 = 1;
    }

    return;

}

int hunter2(const double *arr, const int narr, const double val) {

    /*
     * This function returns index, il, of an array, arr, such that
     * arr[il] < val < arr[il + 1]
     */

    int il = 0;
    int ih = narr;

    while (il != (ih - 1)){
        int im = (il + ih) / 2;
        if   (val <= arr[im])
            ih = im;
        else
            il = im;
    }
    return il;
}


int hunter(const double *arr, const int narr, const double val) {

    /*
     * This function returns index, il, of an array, arr, such that
     * arr[il] < val < arr[il + 1]
     */

    int il, ir, shift;
    double arrb, arre;
    double vl, vr;
    double ldelta, rdelta, delta, dprod;
    double clamped_val;

    clamped_val = val;


    /* Beginning and end values */
    arrb = arr[0];
    arre = arr[narr - 1];

    /* Bounds check */
    if (clamped_val < arrb) {
        print("Warning: interpolation.c: hunter: interpolation out of (lower) bounds.\n");
        print("Clamping to lowest value\n");
        print("val :  %e\n", clamped_val);
        print("arrb:  %e\n", arrb);
        /* clamped_val = arrb; */
        return 0;
    }
    else if (clamped_val > arre) {
        print("Error!: interpolation.c: hunter: interpolation out of (upper) bounds.\n");
        print("val  %e\n", clamped_val);
        print("arre %e\n", arre);
        QUIT_PLUTO(1);
    }

    /* Initial linear guess, left and right indices */
    il = (int) (narr - 1) * (clamped_val - arrb) / (arre - arrb);
    ir = il + 1;

    /* Left and right values from arrays */
    vl = arr[il];
    vr = arr[ir];

    /* Partial and total differences in r */
    ldelta = clamped_val - vl;
    rdelta = vr - clamped_val;
    delta = vr - vl;

    /* If positive, then we're in the right interval */
    dprod = ldelta * rdelta;

    while (dprod < 0) {
        /* Calculate single cell shift direction */
        shift = sgn(ldelta);

        /* New indices */
        il += shift;
        ir = il + 1;

        /* Left and right values from arrays */
        vl = arr[il];
        vr = arr[ir];

        /* Partial and total differences in r */
        ldelta = clamped_val - vl;
        rdelta = vr - clamped_val;
        delta = vr - vl;

        /* If positive, then we're in the right interval */
        dprod = ldelta * rdelta;
    }

    return il;

}

double LinearInterpolate(double y1, double y2, double fc) {
    /* Linear interpolator */

    return (y1 * (1 - fc) + y2 * fc);
}


int sgn(double val) {
    /* The signum function */

    if (val > 0) return 1;

    else if (val < 0) return -1;

    else return 0;
}



double InterpolationWrapper(const double arg_arr[], const double val_arr[], const int narr, const double arg) {

/* Wrapper around hunter and CubicCatmullRomInterpolate */

    int il;
    double arg1, arg2, frac;
    double val0, val1, val2, val3;

    /* Find cell left index to interpolate at */
    il = hunter(arg_arr, narr, arg);

    /* Linear fractional location of arg in cell */
    arg1 = arg_arr[il];
    arg2 = arg_arr[il + 1];
    frac = (arg - arg1) / (arg2 - arg1);

    /* Interpolation */
    val0 = val_arr[il - 1];
    val1 = val_arr[il];
    val2 = val_arr[il + 1];
    val3 = val_arr[il + 2];

    return LinearInterpolate(val1, val2, frac);

}
