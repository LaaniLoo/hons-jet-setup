#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
    return;
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
    Image *image;

    /* Enable our PNG outputs */
    SetOutputVar("rho", PNG_OUTPUT, YES);
    SetOutputVar("prs", PNG_OUTPUT, YES);

    image = GetImage("rho");
    image->slice_plane = X13_PLANE; /* Slice in the X-Z plane */
    image->slice_coord = 0.0; /* Midplane (Y=0) slice */
    // image->max = image->min = 0.0; /* Enable autoscaling for the variable */
    image->max = pow(10,-2); // This is the actual data value, not the logged value
    image->min = pow(10,-7);
    image->logscale = 1; /* Log scaling */
    image->colormap = "red"; /* Red colormap */

    image = GetImage("prs");
    image->slice_plane = X13_PLANE; /* Slice in the X-Z plane */
    image->slice_coord = 0.0; /* Midplane (Y=0) slice */
    // image->max = image->min = 0.0; /* Enable autoscaling for the variable */
    image->max = pow(10,-6);
    image->min = pow(10,-8);
    image->logscale = 1; /* Log scaling */
    image->colormap = "blue"; /* Blue colormap */

#if DIMENSIONS == 3
    SetOutputVar("vx3", PNG_OUTPUT, YES);

    image = GetImage("vx3");
    image->slice_plane = X13_PLANE; /* Slice in the X-Z plane */
    image->slice_coord = 0.0; /* Midplane (Y=0) slice */
    // image->max = image->min = 0.0; /* Enable autoscaling for the variable */
    image->max = 1.0;
    image->min = -1.0;
    image->logscale = 0; /* Normal scaling */
    image->colormap = "br"; /* Blue-Red colormap */
#endif
}
