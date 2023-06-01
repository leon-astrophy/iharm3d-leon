//******************************************************************************
//*                                                                            *
//* FLUXES.C                                                                   *
//*                                                                            *
//* CALCULATES FLUID FLUXES                                                    *
//*                                                                            *
//******************************************************************************

//include header files
#include "decs.h"

//define functions
void lr_to_flux(struct GridGeom *G, struct FluidState *Sl,
struct FluidState *Sr, int dir, int loc, GridPrim *flux, GridVector *ctop);
double ndt_min(GridVector *ctop);

//******************************************************************************

//find time step
double ndt_min(GridVector *ctop) {

  // count time
  timer_start(TIMER_CMAX);

  // initialize time step
  double ndt_min = 1e20;

#if DEBUG
  int min_x1, min_x2, min_x3;
#endif

  // dt according to cfl conditions
#pragma omp parallel for collapse(3) reduction(min:ndt_min)
  ZLOOP {
    double ndt_zone = 0;
    for (int mu = 1; mu < NDIM; mu++) {
      ndt_zone += 1/(cour*dx[mu]/(*ctop)[mu][k][j][i]);
    }
    ndt_zone = 1/ndt_zone;

    if(ndt_zone < ndt_min) {
      ndt_min = ndt_zone;
#if DEBUG
      min_x1 = i; min_x2 = j; min_x3 = k;
#endif
    }
  }

#if DEBUG
  fprintf(stderr, "Timestep set by %d %d %d\n",min_x1,min_x2,min_x3);
#endif

  // count time
  timer_stop(TIMER_CMAX);

  // output dt
  return ndt_min;
}

//******************************************************************************

// calculate flux vector
double get_flux(struct GridGeom *G, struct FluidState *S, struct FluidFlux *F)
{
  //declare
  static struct FluidState *Sl, *Sr;
  static GridVector *ctop;
  double cmax[NDIM], ndts[NDIM];

  //allocate variables
  memset(cmax, 0, NDIM*sizeof(double));
  memset(ndts, 0, NDIM*sizeof(double));

  //allocate variables
  static int firstc = 1;
  if (firstc) {
    Sl  = calloc(1,sizeof(struct FluidState));
    Sr  = calloc(1,sizeof(struct FluidState));
    ctop = calloc(1,sizeof(GridVector));

    firstc = 0;
  }

  //////////////////////////
  //FLAG("First get_flux");
  //////////////////////////

  // reconstruct X-direction
  reconstruct(S, Sl->P, Sr->P, 1);

  //////////////////////////
  //FLAG("After Reconstruct");
  //////////////////////////

  // compute interface fluxes using riemann solvers, X-direction
  lr_to_flux(G, Sl, Sr, 1, FACE1, &(F->X1), ctop);

  // reconstruct Y-direction
  reconstruct(S, Sl->P, Sr->P, 2);

  // compute interface fluxes using riemann solvers, Y-direction
  lr_to_flux(G, Sl, Sr, 2, FACE2, &(F->X2), ctop);

  // reconstruct Z-direction
  reconstruct(S, Sl->P, Sr->P, 3);

  // compute interface fluxes using riemann solvers, Z-direction
  lr_to_flux(G, Sl, Sr, 3, FACE3, &(F->X3), ctop);

  // TODO don't call for static timestep
  return ndt_min(ctop);
}

//******************************************************************************

// Note that the sense of L/R flips from zone to interface during function call
void lr_to_flux(struct GridGeom *G, struct FluidState *Sr,
  struct FluidState *Sl, int dir, int loc, GridPrim *flux, GridVector *ctop)
{
  // count time
  timer_start(TIMER_LR_TO_F);

  // declare
  static GridPrim *fluxL, *fluxR;
  static GridDouble *cmaxL, *cmaxR, *cminL, *cminR, *cmax, *cmin;

  // allocate arrays
  static int firstc = 1;
  if (firstc) {
    fluxL = calloc(1,sizeof(GridPrim));
    fluxR = calloc(1,sizeof(GridPrim));
    cmaxL = calloc(1,sizeof(GridDouble));
    cmaxR = calloc(1,sizeof(GridDouble));
    cminL = calloc(1,sizeof(GridDouble));
    cminR = calloc(1,sizeof(GridDouble));
    cmax = calloc(1,sizeof(GridDouble));
    cmin = calloc(1,sizeof(GridDouble));

    firstc = 0;
  }

  ////////////////////
  //FLAG("First LR");
  ////////////////////

  // Properly offset left face
  // These are un-macro'd to bundle OpenMP thread tasks rather than memory accesses
  PLOOP {
    if (dir == 1) {
#pragma omp parallel for collapse(2)
      ZSLOOP_REVERSE(-1, N3, -1, N2, -1, N1)
        Sl->P[ip][k][j][i] = Sl->P[ip][k][j][i - 1];
    } else if (dir == 2) {
#pragma omp parallel for collapse(2)
      for (int k = (N3) + NG; k >= (-1) + NG; k--) {
        for (int i = (N1) + NG; i >= (-1) + NG; i--) {
          for (int j = (N2) + NG; j >= (-1) + NG; j--)
            Sl->P[ip][k][j][i] = Sl->P[ip][k][j - 1][i];
        }
      }
    } else if (dir == 3) {
#pragma omp parallel for collapse(2)
      for (int j = (N2) + NG; j >= (-1) + NG; j--) {
        for (int i = (N1) + NG; i >= (-1) + NG; i--) {
          for (int k = (N3) + NG; k >= (-1) + NG; k--)
            Sl->P[ip][k][j][i] = Sl->P[ip][k - 1][j][i];
        }
      }
    }
  }
  ////////////////////////////
  //FLAG("Left Face Offset");
  ////////////////////////////

  //count time
  timer_start(TIMER_LR_STATE);

  // Calculate ucon, ucov, bcon, bcov
  get_state_vec(G, Sl, loc, -1, N3, -1, N2, -1, N1);
  get_state_vec(G, Sr, loc, -1, N3, -1, N2, -1, N1);

  //count time
  timer_stop(TIMER_LR_STATE);

  ////////////////////////////
  //FLAG("Left Face Offset");
  ////////////////////////////

  //count time
  timer_start(TIMER_LR_PTOF);

  // Calculate conservative variables, left state
  prim_to_flux_vec(G, Sl, 0,   loc, -1, N3, -1, N2, -1, N1, Sl->U);
  // Calculate fluxes, left state
  prim_to_flux_vec(G, Sl, dir, loc, -1, N3, -1, N2, -1, N1, *fluxL);

  // Calculate conservative variables, left state
  prim_to_flux_vec(G, Sr, 0,   loc, -1, N3, -1, N2, -1, N1, Sr->U);
  // Calculate fluxes, left state
  prim_to_flux_vec(G, Sr, dir, loc, -1, N3, -1, N2, -1, N1, *fluxR);

  //count time
  timer_stop(TIMER_LR_PTOF);

  ////////////////////////////
  //FLAG("State, flux");
  ////////////////////////////

  //count time
  timer_start(TIMER_LR_VCHAR);

  ////////////////////////////////////////////////////////
  // TODO vectorizing these loops fails for some reason
  ////////////////////////////////////////////////////////
  // get magnetosonic speed
#pragma omp parallel
  {
#pragma omp for collapse(2) nowait
    ZSLOOP(-1, N3, -1, N2, -1, N1){
        mhd_vchar(G, Sl, i, j, k, loc, dir, *cmaxL, *cminL);
    }
#pragma omp for collapse(2)
    ZSLOOP(-1, N3, -1, N2, -1, N1) {
        mhd_vchar(G, Sr, i, j, k, loc, dir, *cmaxR, *cminR);
    }
  }

  //count time
  timer_stop(TIMER_LR_VCHAR);

  //count time
  timer_start(TIMER_LR_CMAX);

  //find maximum signal speed across inferfaces
#pragma omp parallel for collapse(2)
  ZSLOOP(-1, N3, -1, N2, -1, N1) {
    (*cmax)[k][j][i] = fabs(MY_MAX(MY_MAX(0., (*cmaxL)[k][j][i]), (*cmaxR)[k][j][i]));
    (*cmin)[k][j][i] = fabs(MY_MAX(MY_MAX(0., -(*cminL)[k][j][i]), -(*cminR)[k][j][i]));
    (*ctop)[dir][k][j][i] = MY_MAX((*cmax)[k][j][i], (*cmin)[k][j][i]);
    if (isnan(1./(*ctop)[dir][k][j][i])) {
      printf("ctop is 0 or NaN at zone: %i %i %i (%i) ", i,j,k,dir);
#if METRIC == MKS
      double X[NDIM];
      double r, th;
      coord(i, j, k, CENT, X);
      bl_coord(X, &r, &th);
      printf("(r,th,phi = %f %f %f)\n", r, th, X[3]);
#endif
      printf("\n");
      exit(-1);
    }
  }

  //count time
  timer_stop(TIMER_LR_CMAX);

  //count time
  timer_start(TIMER_LR_FLUX);

  //interface fluxes, lax-friedrichs method
#pragma omp parallel for simd collapse(3)
  PLOOP {
    ZSLOOP(-1, N3, -1, N2, -1, N1) {
      (*flux)[ip][k][j][i] = 0.5*((*fluxL)[ip][k][j][i] + (*fluxR)[ip][k][j][i] -
               (*ctop)[dir][k][j][i]*(Sr->U[ip][k][j][i] - Sl->U[ip][k][j][i]));
    }
  }
  //count time
  timer_stop(TIMER_LR_FLUX);

  //count time
  timer_stop(TIMER_LR_TO_F);
}

//******************************************************************************

// flux-ct scheme for evolving magnetic field
void flux_ct(struct FluidFlux *F)
{
  //count time
  timer_start(TIMER_FLUX_CT);

  //declare
  static struct FluidEMF *emf;

  //allocate arrays
  static int firstc = 1;
  if (firstc) {
    emf = calloc(1,sizeof(struct FluidEMF));
    firstc = 0;
  }

  //first, calculate emf from fluxes
#pragma omp parallel
  {
    // This and the following are /not/ just ZLOOPs
#pragma omp for simd collapse(2)
    ZSLOOP(0, N3, 0, N2, 0, N1) {
      emf->X3[k][j][i] =  0.25*(F->X1[B2][k][j][i] + F->X1[B2][k][j-1][i]
                              - F->X2[B1][k][j][i] - F->X2[B1][k][j][i-1]);
      emf->X2[k][j][i] = -0.25*(F->X1[B3][k][j][i] + F->X1[B3][k-1][j][i]
                              - F->X3[B1][k][j][i] - F->X3[B1][k][j][i-1]);
      emf->X1[k][j][i] =  0.25*(F->X2[B3][k][j][i] + F->X2[B3][k-1][j][i]
                              - F->X3[B2][k][j][i] - F->X3[B2][k][j-1][i]);
    }

    // Then, Rewrite EMFs as fluxes, after Toth
#pragma omp for simd collapse(2) nowait
    ZSLOOP(0, N3 - 1, 0, N2 - 1, 0, N1) {
      F->X1[B1][k][j][i] =  0.;
      F->X1[B2][k][j][i] =  0.5*(emf->X3[k][j][i] + emf->X3[k][j+1][i]);
      F->X1[B3][k][j][i] = -0.5*(emf->X2[k][j][i] + emf->X2[k+1][j][i]);
    }
#pragma omp for simd collapse(2) nowait
    ZSLOOP(0, N3 - 1, 0, N2, 0, N1 - 1) {
      F->X2[B1][k][j][i] = -0.5*(emf->X3[k][j][i] + emf->X3[k][j][i+1]);
      F->X2[B2][k][j][i] =  0.;
      F->X2[B3][k][j][i] =  0.5*(emf->X1[k][j][i] + emf->X1[k+1][j][i]);
    }
#pragma omp for simd collapse(2)
    ZSLOOP(0, N3, 0, N2 - 1, 0, N1 - 1) {
      F->X3[B1][k][j][i] =  0.5*(emf->X2[k][j][i] + emf->X2[k][j][i+1]);
      F->X3[B2][k][j][i] = -0.5*(emf->X1[k][j][i] + emf->X1[k][j+1][i]);
      F->X3[B3][k][j][i] =  0.;
    }
  } // omp parallel

  //count time
  timer_stop(TIMER_FLUX_CT);
}
