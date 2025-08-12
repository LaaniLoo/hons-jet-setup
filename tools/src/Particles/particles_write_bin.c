/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Writer for particles binary data in .dbl, .flt and 
        ASCII in .tab (only for serial version) format .
 
 \authors A. Mignone (mignone@ph.unito.it)\n
          B. Vaidya (bvaidya@unito.it)\n
 
 \date    April 08, 2018
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Particles_WriteBinary(particleNode *PHeadRef, double dt_particles,
                           Output *output, char *filename)
/*! 
 *  Write particle data in single or double precision binary format.
 *  The binary file structure is:
 *
 *  <Header (ASCII) section>
 *   .
 *   .
 *   .
 *  {field1, field2, ...}_p
 *   .
 *   .
 *   . 
 *
 * All fields are mapped from the particle structure into an array
 * and thus converted to either single or double precision, depending
 * on the function that is being called.
 * Fields can be scalars or array: for each particle nelem indicates
 * the number of elements written for each particle.
 * The field order does not have to match the structure member order but it
 * must be consistent with the sequence declared in Particles_SetOutput()
 * function.
 * Individual fields may be excluded by calling SetOutputVar() 
 * from ChangeOutputVar ().
 *
 *  \param [in]  PheadRef      Pointer to the Head Node of Particle List.
 *  \param [in]  dt_particles  Particle time step
 *  \param [in]  output        Pointer to output structure
 *  \param [in]  filename      File name of particle data: particles.nnnn.flt,
 *                             where nnnn is the data file number.
 ************************************************************************* */
{
  char     fheader[1024];
  size_t   size;
  int      dir, nv, k, nfields, nelem;
  int     *dump_var = output->dump_var;
  int      nvar      = output->nvar;
  long int i, nparticles_glob, proc_npart[g_nprocs+1];
  int     *proc_count, *proc_displacement;
  void    *arr;
  float   *farr;
  double  *darr, u[3], gamma;
  double  *glob_darr;
  float   *glob_farr;
  FILE    *file_handle;
  particleNode * CurNode;
  
/* --------------------------------------------------------
   0. Allocate memory for required fields
   -------------------------------------------------------- */

  nfields = 0; /* Count how many fields are written to disk */
  nelem   = 0; /* The number of elements to be written (some fields may
                  be arrays). */
  for (nv = 0; nv < nvar; nv++) {
    if (dump_var[nv]) {
      nfields++;
      nelem += output->field_dim[nv];
    }
  }
  darr = ARRAY_1D(nelem*p_nparticles, double);
  farr = ARRAY_1D(nelem*p_nparticles, float);
  
/* --------------------------------------------------------
   1. Loop over particles and map struct. members to array
   -------------------------------------------------------- */
    
  i  = 0;  /* Array index */

//print ("\n\n   >> nfields = %d, nelem = %d\n",nfields, nelem);
//print ("    >> PARTICLES_LP_NEBINS, NFLX = %d, %d\n",PARTICLES_LP_NEBINS,NFLX);
  PARTICLES_LOOP(CurNode, PHeadRef){

  /* -- 1a. Compute three- or four-velocity -- */

    for (dir = 0; dir < 3; dir++) u[dir] = CurNode->p.speed[dir];
    #if PARTICLES_TYPE == COSMIC_RAYS && PARTICLES_CR_WRITE_4VEL == YES
    gamma =   CurNode->p.speed[IDIR]*CurNode->p.speed[IDIR]
            + CurNode->p.speed[JDIR]*CurNode->p.speed[JDIR]
            + CurNode->p.speed[KDIR]*CurNode->p.speed[KDIR];
    gamma = 1.0/sqrt(1.0 - gamma/(PARTICLES_CR_C*PARTICLES_CR_C));
    for (dir = 0; dir < 3; dir++) u[dir] *= gamma;  /* 4-vel */
    #endif

  /* ------------------------------------------------------
     1b. Map structure members to array.
         Important: field order should match the order 
         given in Particles_SetOutput().
         Here nv scan *all* fields (nv <= nfields)
     ------------------------------------------------------ */

    nv = 0;
    #if PARTICLES_TYPE == COSMIC_RAYS
    if (dump_var[nv++]) darr[i++] = CurNode->p.id;
    if (dump_var[nv++]) darr[i++] = CurNode->p.coord[IDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.coord[JDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.coord[KDIR];
    if (dump_var[nv++]) darr[i++] = u[IDIR];
    if (dump_var[nv++]) darr[i++] = u[JDIR];
    if (dump_var[nv++]) darr[i++] = u[KDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.rho;
    if (dump_var[nv++]) darr[i++] = CurNode->p.tinj;
    if (dump_var[nv++]) darr[i++] = CurNode->p.color;
    #endif

    #if PARTICLES_TYPE == LAGRANGIAN
    if (dump_var[nv++]) darr[i++] = CurNode->p.id;
    if (dump_var[nv++]) darr[i++] = CurNode->p.coord[IDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.coord[JDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.coord[KDIR];
    if (dump_var[nv++]) darr[i++] = u[IDIR];
    if (dump_var[nv++]) darr[i++] = u[JDIR];
    if (dump_var[nv++]) darr[i++] = u[KDIR];
    if (dump_var[nv++]) darr[i++] = CurNode->p.tinj;
    if (dump_var[nv++]) darr[i++] = CurNode->p.color;
    #if PARTICLES_LP_RAISE == YES
    if (dump_var[nv++]){
      for(k = 0; k < NTRACER; k++) {
        darr[i++] = CurNode->p.tracer[k];
      }
    }

    if (dump_var[nv++]){
      for(k = 0; k < PARTICLES_LP_SHK_BINS; k++) {
        darr[i++] = CurNode->p.last_shock_time[k];
      }
    }

    if (dump_var[nv++]) darr[i++] = CurNode->p.density; 
    if (dump_var[nv++]) darr[i++] = CurNode->p.pressure; 
    #elif PARTICLES_LP_SPECTRA == YES
    if (dump_var[nv++]) darr[i++] = CurNode->p.density; 
    if (dump_var[nv++]) darr[i++] = CurNode->p.nmicro;
    if (dump_var[nv++]) darr[i++] = CurNode->p.cmp_ratio;
    if (dump_var[nv++]) darr[i++] = CurNode->p.shkflag;
    if (dump_var[nv++]) darr[i++] = CurNode->p.shk_gradp;
    if (dump_var[nv++]) darr[i++] = CurNode->p.ca;
    if (dump_var[nv++]) darr[i++] = CurNode->p.cr;

    if (dump_var[nv++]) {
      for (k = 0; k < NFLX; k++) darr[i++] = CurNode->p.shk_vL[k];
    }

    if (dump_var[nv++]) {
      for (k = 0; k < NFLX; k++) darr[i++] = CurNode->p.shk_vR[k];
    }

    if (dump_var[nv++]){
      for(k = 0; k < PARTICLES_LP_NEBINS; k++) {
        darr[i++] = CurNode->p.eng[k];
      }
    }
    
    if (dump_var[nv++]){
      for(k = 0; k < PARTICLES_LP_NEBINS; k++){
        darr[i++] = CurNode->p.chi[k];
      }
    }

    #endif
    #endif
    
  }
  
/* --------------------------------------------------------
   2. Compute the total number of particles and gather
      each proc. particle number into the proc_npart[1:]
      array. This is used to calculate displacement and count arrays.
      Also gather particle data to write using a single process.
   -------------------------------------------------------- */
    
#ifdef PARALLEL
  MPI_Allreduce(&p_nparticles, &nparticles_glob, 1,
                MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allgather(&p_nparticles, 1, MPI_LONG, &(proc_npart[1]), 1,
                MPI_LONG, MPI_COMM_WORLD);


  /* Set up root process for writing data */
  if (prank == 0) {
      /* Set up count and displacement arrays */
      proc_count = ARRAY_1D(g_nprocs, int);
      proc_displacement = ARRAY_1D(g_nprocs, int);

      /* Set up counts array for Gatherv operation */
      for (i=0; i<g_nprocs; i++) proc_count[i] = proc_npart[i+1]*nelem;

      /* Set up displacements array for Gatherv operation */
      proc_displacement[0] = 0;
      for (i=1; i<g_nprocs; i++) proc_displacement[i] = proc_displacement[i-1] + proc_count[i-1];

      /* Create global darr, that can hold ALL particles */
      glob_darr = ARRAY_1D(nelem*nparticles_glob, double);

      /* Create global farr if it's going to be needed later */
      if (output->type == PARTICLES_FLT_OUTPUT) glob_farr = ARRAY_1D(nelem*nparticles_glob, float);
  }

  /* Gather all particle data to root process.
   * send buffer is the local darr, send count is the local p_nparticles multiplied by nelem.
   * datatype is double. Recieved buffer (on root) is glob_darr, recieved counts is proc_count,
   * displacements is proc_displacement. Received datatype is also double, root is 0, communicator is world.
   */
  MPI_Gatherv(darr, nelem*p_nparticles, MPI_DOUBLE, glob_darr, proc_count, proc_displacement, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
  nparticles_glob = p_nparticles;
  glob_darr = darr;
  glob_farr = farr;
#endif

#ifdef PARALLEL
  /* Only the root process participates in the remainder for parallel operation */
  if (prank == 0)
#endif
  {

      /* Open our file handle */
      file_handle = FileOpenSingleProcess(filename, "wb");

    /* --------------------------------------------------------
       3. Write file header section.
       -------------------------------------------------------- */
      
      sprintf(fheader,"# PLUTO %s binary particle data file\n", PLUTO_VERSION);
      sprintf(fheader+strlen(fheader),"# dimensions     %d\n",DIMENSIONS);
      sprintf(fheader+strlen(fheader),"# nflux          %d\n",NFLX);
      sprintf(fheader+strlen(fheader),"# dt_particles   %12.6e\n",dt_particles);  
      if (IsLittleEndian())
          sprintf(fheader+strlen(fheader),"# endianity      little\n");
      else
          sprintf(fheader+strlen(fheader),"# endianity      big\n");

      sprintf(fheader+strlen(fheader),"# nparticles     %ld\n",nparticles_glob);  
      sprintf(fheader+strlen(fheader),"# idCounter      %ld\n",p_idCounter);
      sprintf(fheader+strlen(fheader),"# particletype   %d\n",PARTICLES_TYPE);
      if (output->type == PARTICLES_FLT_OUTPUT){  
        sprintf(fheader+strlen(fheader),"# precision      float\n");    
      }else if (output->type == PARTICLES_DBL_OUTPUT){  
        sprintf(fheader+strlen(fheader),"# precision      double\n");    
      }
      
      sprintf(fheader+strlen(fheader),"# time           %12.12e\n",g_time);  
      sprintf(fheader+strlen(fheader),"# stepNumber     %ld\n",g_stepNumber);
      sprintf(fheader+strlen(fheader),"# nfields        %d\n",nfields);
      sprintf(fheader+strlen(fheader),"# nelem          %d\n",nelem);
      sprintf(fheader+strlen(fheader),"# field_names    ");
      for (i = 0; i < nvar; i++){
        if (dump_var[i]) {
          sprintf(fheader+strlen(fheader),output->var_name[i]);
          sprintf(fheader+strlen(fheader),"  ");
        }  
      }
      sprintf(fheader+strlen(fheader),"\n");

      sprintf(fheader+strlen(fheader),"# field_dim      ");
      for (i = 0; i < nvar; i++){
        if (dump_var[i]) {
          sprintf(fheader+strlen(fheader),"%d",output->field_dim[i]);
          sprintf(fheader+strlen(fheader),"  ");
        }  
      }
      sprintf(fheader+strlen(fheader),"\n");

      sprintf(fheader+strlen(fheader),"# shk_thresh     ");
      for (i = 0; i < PARTICLES_LP_SHK_BINS; i++){
        sprintf(fheader+strlen(fheader),"%f",pow(10, ln_shk_min + dln_shk * i));
        sprintf(fheader+strlen(fheader),"  ");
      }
      sprintf(fheader+strlen(fheader),"\n");

      /* Only root writes, so we use single process function */
      FileWriteHeaderSingleProcess(fheader, file_handle);
        
    /* --------------------------------------------------------
       4. Write data 
       -------------------------------------------------------- */
     
      if (output->type == PARTICLES_DBL_OUTPUT) size = sizeof(double);
      if (output->type == PARTICLES_FLT_OUTPUT) size = sizeof(float);

      if (output->type == PARTICLES_FLT_OUTPUT){
        for (i = 0; i < nelem*nparticles_glob; i++) glob_farr[i] = (float)glob_darr[i];
        arr = (void *) glob_farr;
      }else if (output->type == PARTICLES_DBL_OUTPUT){
        arr = (void *) glob_darr;
      }

      /* Only root writes, so we use single process function */
      FileWriteArraySingleProcess(arr, nelem*nparticles_glob, size, file_handle);
      
      /* Close file handle */
      FileCloseSingleProcess(file_handle);
  }

  FreeArray1D(darr);   
  FreeArray1D(farr);   

#ifdef PARALLEL
  /* Free global arrays if we're root */
  if (prank == 0) {
      FreeArray1D(proc_count);
      FreeArray1D(proc_displacement);
      FreeArray1D(glob_darr);
      if (output->type == PARTICLES_FLT_OUTPUT) FreeArray1D(glob_farr);
  }
#endif
}

/* ********************************************************************* */
void Particles_WriteTab(particleNode* PHeadRef, char filename[128])
/*
 * Write particle coordinates, ids and speeds in Tab ASCII format 
 * only for *serial* version. 
 *
 *  \param [in]  PheadRef    Pointer to the Head Node of Particle List.
 *  \param [in]  filename    File name of particle data: particles.nnnn.tab,
 *                           where nnnn is the data file number.
 *********************************************************************** */
{
#ifdef PARALLEL
  print("! WARNING: Particle Data in Tabulated ASCII format is only written \
            with serial version");
#else
  int i;
  FILE *stream;
  particleNode* CurNode;
  
  stream = fopen(filename, "w");
  fprintf(stream, "# Nparticles: %ld\n", p_nparticles);
  fprintf(stream, "# Step:       %ld\n", g_stepNumber);
  fprintf(stream, "# Time:       %f\n",  g_time);
  
  
  CurNode = PHeadRef;

  while(CurNode != NULL) {
    fprintf(stream, "%d", CurNode->p.id);
    
    for (i = 0; i < 3; ++i) {
      fprintf(stream, "  %lf", CurNode->p.coord[i]);
    }
    for (i = 0; i < 3; ++i) {
      fprintf(stream, "  %lf", CurNode->p.speed[i]);
    }
    fprintf(stream, "\n");

    CurNode = CurNode->next;
  }
  fclose(stream);
#endif
}

