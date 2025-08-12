/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Restart the simulation using particles.nnnn.dbl file.
  
 To speed up restart, particles are read in blocks of \c nchunk
 elements.

 \authors B. Vaidya (bvaidya@unito.it)\n
          A. Mignone (mignone@ph.unito.it)\n
 
  \date   March 31, 2018
 */
/* ///////////////////////////////////////////////////////////////////// */
#include <pluto.h>

/* ********************************************************************* */
void Particles_Restart(Runtime *ini, Data *d, Grid *grid)
/*!
*  Routine for restarting from particles files.
*  
*
* \param [in]   d          Pointer to the PLUTO data structure.
* \param [in]  nrestart    number of the restart data file 
*                          given from command line
* \param [in]  grid        Pointer to the PLUTO grid structure.
************************************************************************  */
{
  int j, k, dir, np, count, nv, n;
  int nfields;     /* The number of structure fields to be read */
  int nelem;       /* The total number of field elements per particle */
  long int np_tot; /* Total number of particles to be read */
  long int nchunk; /* Number of particles in a chunk of data */
  long int off;
  char filename[512], line[512];
  Particle p;
  double u2, gamma;
  static double **buf;
  clock_t t_rbeg, t_rend;
  FILE *fp;
  Output *particle_output = 0;

  double particle_restart_time;
  double time_delta = 0.0;
  
  t_rbeg = clock();

  p_nparticles = 0; /* Initialize it to zero before reading particles */

/* ----------------------------------------------------------
   0.0 Initialise particle shock capturing and data structures
   ---------------------------------------------------------- */ 

  /* Initialse shock structures */
  Particles_LP_InitShock();

  /* Initialise MPI particle data type */
#ifdef PARALLEL
  Particles_StructDatatype();
#endif  

/* -----------------------------------------------------
   0.5 Build file name and read header file 
   ----------------------------------------------------- */

  /* Set p_nrestart to invalid value, to catch errors if we don't have particles to restart from */
  p_nrestart = -1;

  /* Find our particle output structure */
  for (n = 0; n < MAX_OUTPUT_TYPES; n++) {
      if (ini->output[n].type == PARTICLES_DBL_OUTPUT) {
          particle_output = &(ini->output[n]);
          break;
      }
  }

  /* Check that we found a particle_output */
  if (!particle_output) {
      print("> ERROR finding dbl particle output in runtime structure.\n");
      print("> Restarting with particles only supports .dbl outputs.\n");
      print("> Continuing without restarting particles!!!.\n");
      return;
  }

  p_nrestart = particle_output->nfile;

  if (p_nrestart < 0) {
      print("> No particle outputs to restart from.\n");
      print("> Continuing without restarting particles!!!.\n");
      return;
  }
     
  sprintf (filename,"%s/particles.%04d.dbl", particle_output->dir, p_nrestart);

  if (prank == 0){
    char A0[256], A1[256], A2[256];
    fp = fopen(filename,"rb");
    if (fp == NULL){
      print ("! Particles_Restart(): file %s does not exist.\n",filename);
      QUIT_PLUTO(1);
    }else{
      while(fgets(line,256,fp) && line[0] == '#'){
        sscanf(line, "%s %s %s",A0, A1, A2);
        if (strcmp(A1, "nparticles")   == 0) np_tot      = atoi(A2);
        if (strcmp(A1, "nfields")      == 0) nfields     = atoi(A2);
        if (strcmp(A1, "nelem")        == 0) nelem       = atoi(A2);
        if (strcmp(A1, "idCounter")    == 0) p_idCounter = atoi(A2);
        if (strcmp(A1, "dt_particles") == 0) d->Dts->invDt_particles = 0.99/atof(A2);
        if (strcmp(A1, "time")         == 0) particle_restart_time = atof(A2);
        off = ftell(fp);
      }
    }
    fclose(fp);

    /* Check that our particle restart time is close enough to our actual time */
    /* We only need to check this for thie first process */
    time_delta = g_time - particle_restart_time;
    print("> Time difference between particle restart and sim time is %12.12e\n", time_delta);

    if (fabs(time_delta) > PARTICLES_RESTART_THRESHOLD) {
        print("> ERROR! Difference between particle output and simulation time is too large.");
        print("> Continuing without restarting particles!!!.\n");
        return;
    }
  }  
 
#ifdef PARALLEL
/* --------------------------------------------------------
   1a. Broadcast the file offset and No. of particles 
       read from file to all processor and set view based
       on these values for parallel reading.
       The variable p_nparticles should not be set here
       but during the call to Particles_Insert().
   -------------------------------------------------------- */

  MPI_Datatype ParticleFields;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&off,         1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&np_tot,      1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nfields,     1,      MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nelem,       1,      MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&p_idCounter, 1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&(d->Dts->invDt_particles), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
  MPI_Type_contiguous(nelem, MPI_DOUBLE, &ParticleFields);
  MPI_Type_commit(&ParticleFields);

  MPI_File fres;
  MPI_Status stats;
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fres); 
  MPI_File_set_view(fres, off, ParticleFields, ParticleFields, "native", MPI_INFO_NULL);
#else
  fp = fopen(filename, "rb");
  fseek(fp, off, SEEK_SET);
#endif

  print("> Restarting from Particle File : %s, idCounter = %d\n",
          filename, p_idCounter);

/* -------------------------------------------------------
   1b. Each processor reads the whole file in chunks of
       nchunk particles.
       Particles are insertted to the linked list only
       if they lie on the computational domain.
   ------------------------------------------------------- */

  nchunk = np_tot/g_nprocs;
#ifdef PARALLEL
  int pw = 1;   /* Round it to closest power of 2 */
  while (pw < nchunk) pw *= 2;
  nchunk = 4*pw;  /* Increasing the size may cause memory overflow,
                   * but gives much better performances for large number
                   * of processors.                                    */
#endif
  
  print("  building linked list of %d particles [nchunk = %d]\n",np_tot, nchunk);
  if (buf == NULL) buf = ARRAY_2D(nchunk, nelem, double); 

  for (count = 1; count <= np_tot; count += nchunk){
    np = nchunk;  /* Number of particles to be read */
    if ( (count + nchunk) > np_tot) np = np_tot - count + 1;
    #ifdef PARALLEL
    MPI_File_read_all(fres, buf[0], np, ParticleFields, MPI_STATUS_IGNORE);
    #else
    fread(buf[0], sizeof(double), np*nelem, fp);
    #endif
      
  /* ------------------------------------------------------
     1c. Map array elements to structure members.
         Fields are assigned in the same order they're
         written (see the field name order in
         Particles_SetOuput).
         Since restart is done from .dbl files, we recommend
         the fields and their order is not modified.
     ------------------------------------------------------ */

    for (j = 0; j < np; j++){
      k = 0;
      #if PARTICLES_TYPE == COSMIC_RAYS
      p.id          = buf[j][k++];
      p.coord[IDIR] = buf[j][k++];
      p.coord[JDIR] = buf[j][k++];
      p.coord[KDIR] = buf[j][k++];
      p.speed[IDIR] = buf[j][k++];
      p.speed[JDIR] = buf[j][k++];
      p.speed[KDIR] = buf[j][k++];
      p.rho         = buf[j][k++];      
      p.tinj        = buf[j][k++];
      p.color       = buf[j][k++];
      #if PARTICLES_CR_WRITE_4VEL == YES
      u2 =   p.speed[IDIR]*p.speed[IDIR] 
           + p.speed[JDIR]*p.speed[JDIR] 
           + p.speed[KDIR]*p.speed[KDIR]; 
      gamma = sqrt(1.0 + u2/(PARTICLES_CR_C*PARTICLES_CR_C));
      p.speed[IDIR] /= gamma;
      p.speed[JDIR] /= gamma;
      p.speed[KDIR] /= gamma;
      #endif
      #endif
     
      #if PARTICLES_TYPE == LAGRANGIAN
      p.id          = buf[j][k++];
      p.coord[IDIR] = buf[j][k++];
      p.coord[JDIR] = buf[j][k++];
      p.coord[KDIR] = buf[j][k++];
      p.speed[IDIR] = buf[j][k++];
      p.speed[JDIR] = buf[j][k++];
      p.speed[KDIR] = buf[j][k++];
      p.tinj        = buf[j][k++];
      p.color       = buf[j][k++];
      #if PARTICLES_LP_RAISE == YES
      for(nv = 0; nv < NTRACER; nv++) {
        p.tracer[nv] = buf[j][k++];
      }

      for(nv = 0; nv < PARTICLES_LP_SHK_BINS; nv++) {
        p.last_shock_time[nv] = buf[j][k++];
      }
      p.density     = buf[j][k++];
      p.pressure    = buf[j][k++];
      #elif PARTICLES_LP_SPECTRA == YES
      p.density     = buf[j][k++];
      p.nmicro      = buf[j][k++];
      p.cmp_ratio   = buf[j][k++];
      p.shkflag     = buf[j][k++];
      p.shk_gradp   = buf[j][k++];
      p.ca          = buf[j][k++];
      p.cr          = buf[j][k++];
      for (nv = 0; nv < NFLX; nv++) p.shk_vL[nv] = buf[j][k++];
      for (nv = 0; nv < NFLX; nv++) p.shk_vL[nv] = buf[j][k++];
      for(nv = 0; nv < PARTICLES_LP_NEBINS; nv++) {
        p.eng[nv] = buf[j][k++];
      }
      for(nv = 0; nv < PARTICLES_LP_NEBINS; nv++) {
        p.chi[nv] = buf[j][k++];
      }
      #endif
      #endif
 
      Particles_Insert (&p, d, PARTICLES_RESTART, grid);
    } /* End loop on particles */  
  } /* End Loop on chunks */

#ifdef PARALLEL
  MPI_File_close(&fres);
#else
  fclose(fp);
#endif

  t_rend = clock();

 #if PARTICLES_LP_SPECTRA == YES || PARTICLES_LP_RAISE == YES
  particleNode *CurNode = d->PHead;
  Particle *pl;
  PARTICLES_LOOP(CurNode, d->PHead){
    pl = &(CurNode->p);
    Particles_LP_FixValue(pl, d, grid);
    //Particles_Display(&p);
  }
  #endif

  int tot_part = Particles_Number(d->PHead);
  print("> Restarted with %d local particles\n", tot_part);

  print("  restart took: %f (s)\n",((double)(t_rend - t_rbeg)/CLOCKS_PER_SEC));

  FreeArray2D((void **) buf);
}
