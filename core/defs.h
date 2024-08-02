//******************************************************************************
//*                                                                            *
//* DEFS.H                                                                     *
//*                                                                            *
//* GLOBAL VARIABLE DEFINITIONS                                                *
//*                                                                            *
//******************************************************************************

#pragma once

/*---------------------------------------------------------------*/

// Zone flags. TODO move these to the heap
GridInt pflag;
GridInt fail_save;
GridInt fflag;

/*---------------------------------------------------------------*/
/* arrays */

// backup variables
#if DEBUG
struct FluidFlux preserve_F;
GridPrim preserve_dU;
#endif

// Leon's patch, extra variables for cooling //
#if COOLING
GridDouble omg_gr; // angular velocity 
GridDouble t_gr; // target temperature
#endif

/*---------------------------------------------------------------*/
// Parameters

// physical
double a;
double gam;
double Rhor;
double tp_over_te;

// geometry
double Rin, Rout, hslope;
double poly_norm, poly_xt, poly_alpha, mks_smooth;
double cour;
double dV, dx[NDIM], startx[NDIM];
double x1Min, x1Max, x2Min, x2Max, x3Min, x3Max;
double dt, dt_light;
double t, tf;
double rcurr, hcurr;
int istart, istop, jstart, jstop;
int nstep;
int is_restart;

// fluid dumps
double DTd;
double DTf;
double DTl;
int DTr;
int DTp;
int dump_cnt;
double tdump, tlog;

// derived logged output
double mdot, edot, ldot;
double mdot_eh, edot_eh, ldot_eh;
int icurr, jcurr, kcurr;

//electronic variables
#if ELECTRONS
double game, gamp;
double fel0;
double tptemin, tptemax;
#endif

// number of threads
int nthreads;

// MPI stuff 
int global_start[3];
int global_stop[3];

/*---------------------------------------------------------------*/

// Leon's patch, unit coversion //
double R_isco;

// Leon's patch, unit conversion //
double M_unit, mbh; 
double Mbh, L_unit, T_unit, RHO_unit, U_unit, B_unit;

/*---------------------------------------------------------------*/

