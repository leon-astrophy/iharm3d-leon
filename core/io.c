//******************************************************************************
//*                                                                            *
//* IO.C                                                                       *
//*                                                                            *
//* HDF5 OUTPUT AND RESTART                                                    *
//*                                                                            *
//******************************************************************************

//include main headers
#include "decs.h"

//include hdf5 headers
#include "hdf5_utils.h"

//include other headers
#include <sys/stat.h>
#include <ctype.h>

/////////////////////////////////////
// TODO move bsq to this scope? (no)
/////////////////////////////////////

// delcare
#if DEBUG
#define OUT_H5_TYPE H5T_IEEE_F64LE
#define OUT_TYPE double
#else
#define OUT_H5_TYPE H5T_IEEE_F32LE
#define OUT_TYPE float
#endif

// delcare
#define HDF_STR_LEN 20

//***************************************************************************************

//output dump
void dump(struct GridGeom *G, struct FluidState *S)
{
  dump_backend(G, S, IO_REGULAR);
}

//***************************************************************************************

//write dump
void dump_backend(struct GridGeom *G, struct FluidState *S, int type)
{
  //count time
  timer_start(TIMER_IO);

  //declare
  static GridDouble *data;
  char fname[80];

  //choose what to output
  #if ELECTRONS
  #if ALLMODELS
  const char varNames[NVAR][HDF_STR_LEN] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3",
                            "KTOT", "KEL0", "KEL1", "KEL2", "KEL3"};
  #else
  const char varNames[NVAR][HDF_STR_LEN] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3",
                            "KTOT", "KEL0"};
  #endif
  #else
  const char varNames[NVAR][HDF_STR_LEN] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3"}; //Reserve some extra
  #endif

  //allocate arrays
  static int firstc = 1;
  if(firstc) {
    data = calloc(1,sizeof(GridDouble));
    firstc = 0;
  }

  //Don't re-dump the grid after a restart
  if (dump_cnt == 0) dump_grid(G);

  //dump file type
  if (type == IO_REGULAR) {
    sprintf(fname, "dumps/dump_%08d.h5", dump_cnt);
  } else if (type == IO_ABORT) {
    sprintf(fname, "dumps/dump_abort.h5");    
  }

  //mpi stuff
  if(mpi_io_proc()) fprintf(stdout, "DUMP %s\n", fname);
  
  //open hdf5 file
  hdf5_create(fname);

  //Can use H5T_VARIABLE for any-length strings. But not compatible with parallel IO
  hid_t string_type = hdf5_make_str_type(HDF_STR_LEN);

  // Write header
  hdf5_set_directory("/");
  hdf5_make_directory("header");
  hdf5_set_directory("/header/");
  
  //Adding problem-specific data (problem type, initial data)
  hdf5_make_directory("problem");
  hdf5_set_directory("/header/problem/");
  save_problem_data(string_type);
  
  //header files
  hdf5_set_directory("/header/");
  hdf5_write_single_val(VERSION, "version", string_type);
  int has_electrons = ELECTRONS;
  hdf5_write_single_val(&has_electrons, "has_electrons", H5T_STD_I32LE);

  //write metric type
#if METRIC == MINKOWSKI
  hdf5_write_single_val("MINKOWSKI", "metric", string_type);
#elif METRIC == MKS
#if DEREFINE_POLES
  hdf5_write_single_val("MMKS", "metric", string_type);
#else
  hdf5_write_single_val("MKS", "metric", string_type);
#endif //DEREFINE
#endif //MKS

  ////////////////////////////////////////////
  // TODO match below instead of hard-coding
  ////////////////////////////////////////////
  // grid files
  char gridfile[HDF_STR_LEN] = "grid.h5"; 
  hdf5_write_single_val(&gridfile, "gridfile", string_type);

  // write reconstruction scheme
#if RECONSTRUCTION == LINEAR
  hdf5_write_single_val("LINEAR", "reconstruction", string_type);
#elif RECONSTRUCTION == PPM
  hdf5_write_single_val("PPM", "reconstruction", string_type);
#elif RECONSTRUCTION == WENO
  hdf5_write_single_val("WENO", "reconstruction", string_type);
#elif RECONSTRUCTION == MP5
  hdf5_write_single_val("MP5", "reconstruction", string_type);
#endif

  // write n1, n2, n3
  int n1 = N1TOT, n2 = N2TOT, n3 = N3TOT;
  hdf5_write_single_val(&n1, "n1", H5T_STD_I32LE);
  hdf5_write_single_val(&n2, "n2", H5T_STD_I32LE);
  hdf5_write_single_val(&n3, "n3", H5T_STD_I32LE);

  // write number of primitive variables
  int n_prims = NVAR;
  hdf5_write_single_val(&n_prims, "n_prim", H5T_STD_I32LE);

  // In case we do passive variables
  int n_prims_passive = 0;
  hdf5_write_single_val(&n_prims_passive, "n_prims_passive", H5T_STD_I32LE);
  hdf5_write_str_list(varNames, "prim_names", HDF_STR_LEN, n_prims);

  //gas constant
  hdf5_write_single_val(&gam, "gam", H5T_IEEE_F64LE);

  //electronic variables
#if ELECTRONS
  hdf5_write_single_val(&game, "gam_e", H5T_IEEE_F64LE);
  hdf5_write_single_val(&gamp, "gam_p", H5T_IEEE_F64LE);
  hdf5_write_single_val(&tptemin, "tptemin", H5T_IEEE_F64LE);
  hdf5_write_single_val(&tptemax, "tptemax", H5T_IEEE_F64LE);
  hdf5_write_single_val(&fel0, "fel0", H5T_IEEE_F64LE);
#endif
  hdf5_write_single_val(&cour, "cour", H5T_IEEE_F64LE);
  hdf5_write_single_val(&tf, "tf", H5T_IEEE_F64LE);
  hdf5_add_units("tf", "code");

  //Geometry
  hdf5_make_directory("geom");
  hdf5_set_directory("/header/geom/");

  // grid variables
  hdf5_write_single_val(&(startx[1]), "startx1", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(startx[2]), "startx2", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(startx[3]), "startx3", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(dx[1]), "dx1", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(dx[2]), "dx2", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(dx[3]), "dx3", H5T_IEEE_F64LE);

  // number of dimensions
  int n_dim = NDIM;
  hdf5_write_single_val(&n_dim, "n_dim", H5T_STD_I32LE);

  //ker schild coordinates variables
#if METRIC == MKS
#if DEREFINE_POLES
  hdf5_make_directory("mmks");
  hdf5_set_directory("/header/geom/mmks/");
  hdf5_write_single_val(&poly_xt, "poly_xt", H5T_IEEE_F64LE);
  hdf5_write_single_val(&poly_alpha, "poly_alpha", H5T_IEEE_F64LE);
  hdf5_write_single_val(&mks_smooth, "mks_smooth", H5T_IEEE_F64LE);
#else
  hdf5_make_directory("mks");
  hdf5_set_directory("/header/geom/mks/");
#endif
  hdf5_write_single_val(&Rin, "r_in", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Rout, "r_out", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Rhor, "r_eh", H5T_IEEE_F64LE);

  ////////////////////////////////////////////////////////////////
  // I don't need this in code but it's in the spec to output it
  ////////////////////////////////////////////////////////////////
  double z1 = 1 + pow(1 - a*a,1./3.)*(pow(1+a,1./3.) + pow(1-a,1./3.));
  double z2 = sqrt(3*a*a + z1*z1);
  double Risco = 3 + z2 - sqrt((3-z1)*(3 + z1 + 2*z2));
  hdf5_write_single_val(&Risco, "r_isco", H5T_IEEE_F64LE);
  hdf5_write_single_val(&hslope, "hslope", H5T_IEEE_F64LE);
  hdf5_write_single_val(&a, "a", H5T_IEEE_F64LE);
#endif

  // hdf5 file directory
  hdf5_set_directory("/");

  //////////////////////////////////
  // TODO do partial/full dumps
  //////////////////////////////////

  // what are these?
  int is_full_dump = 1; 
  hdf5_write_single_val(&is_full_dump, "is_full_dump", H5T_STD_I32LE);
  hdf5_write_single_val(&t, "t", H5T_IEEE_F64LE);
  hdf5_add_units("t", "code");
  hdf5_write_single_val(&dt, "dt", H5T_IEEE_F64LE);
  hdf5_add_units("dt", "code");
  hdf5_write_single_val(&nstep, "n_step", H5T_STD_I32LE);
  hdf5_write_single_val(&dump_cnt, "n_dump", H5T_STD_I32LE);
  hdf5_write_single_val(&DTd, "dump_cadence", H5T_IEEE_F64LE);
  hdf5_write_single_val(&DTf, "full_dump_cadence", H5T_IEEE_F64LE);

  // Write primitive variables
  pack_write_vector(S->P, NVAR, "prims", OUT_H5_TYPE);
  hdf5_add_units("prims", "code");

  // Write jcon (not recoverable from prims)
  pack_write_vector(S->jcon, NDIM, "jcon", OUT_H5_TYPE);
  hdf5_add_units("jcon", "code");

  // lorentz factors
  ZLOOP (*data)[k][j][i] = mhd_gamma_calc(G, S, i, j, k, CENT);
  pack_write_scalar((*data), "gamma", OUT_H5_TYPE);

  // Space for any extra items
  // Currently debug/diagnostic output, on full dumps only
  hdf5_make_directory("extras");
  hdf5_set_directory("/extras/");

  // Preserve git commit or tag of the run -- see Makefile
#ifdef GIT_VERSION
  hdf5_write_single_val(QUOTE(GIT_VERSION), "git_version", string_type);
#endif

  // write error stuff
  if (is_full_dump) {
    ///////////////////////////////////////////////
    // TODO need sync here for consistent output?
    ///////////////////////////////////////////////
    ZLOOP (*data)[k][j][i] = flux_ct_divb(G, S, i, j, k);
    pack_write_scalar((*data), "divB", OUT_H5_TYPE);

    pack_write_int(fail_save, "fail");
    ZLOOP fail_save[k][j][i] = 0;

    pack_write_int(fflag, "fixup");
  }

  // write conservative varaiables, fluxes, and source terms
#if DEBUG
    pack_write_vector(S->U, NVAR, "U", OUT_H5_TYPE);
    pack_write_vector(preserve_F.X1, NVAR, "X1", OUT_H5_TYPE);
    pack_write_vector(preserve_F.X2, NVAR, "X2", OUT_H5_TYPE);
    pack_write_vector(preserve_F.X3, NVAR, "X3", OUT_H5_TYPE);
    pack_write_vector(preserve_dU, NVAR, "dU", OUT_H5_TYPE);
#endif

  ///////////////////////////////////////////////////////////////////
  // Extra physical variables.  These are just recalculated in post
  // ZLOOP (*data)[k][j][i] = bsq_calc(G, S, i, j, k);
  // pack_write_scalar((*data), "bsq", OUT_H5_TYPE);
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  // Could write extra vectors here, but they take up too much space
  // pack_write_vector(S->bcon, NDIM, "bcon", OUT_H5_TYPE);
  // hdf5_add_units("bcon", "code");
  ///////////////////////////////////////////////////////////////////
  
  //close files
  hdf5_close();

  //count time
  timer_stop(TIMER_IO);
}

//***************************************************************************************

//output grid files
#define NGRIDVARS 9
void dump_grid(struct GridGeom *G)
{
  //allocate variables
  GridDouble *x[NGRIDVARS];
  for (int d = 0; d < NGRIDVARS; d++) x[d] = calloc(1,sizeof(GridDouble));

  //grid file headers
  const char *coordNames[] = {"X", "Y", "Z", "r", "th", "phi", "X1", "X2", "X3"};

  //grid file names
  char *fname = "dumps/grid.h5";
  
  //mpi stuff
  if(mpi_io_proc()) fprintf(stdout, "GRID %s\n", fname);

  //create hdf5
  hdf5_create(fname);

  //hdf5 directory
  hdf5_set_directory("/");

  // Batch fill grid var buffers since lots of them are the same
  ZLOOP {
	double xp[4];
    coord(i, j, k, CENT, xp);
    #if METRIC == MINKOWSKI
    (*x[0])[k][j][i] = xp[1];
    (*x[1])[k][j][i] = xp[2];
    (*x[2])[k][j][i] = xp[3];
    (*x[3])[k][j][i] = 0;
    (*x[4])[k][j][i] = 0;
    (*x[5])[k][j][i] = 0;
    #elif METRIC == MKS
    double r, th;
    bl_coord(xp, &r, &th);
    (*x[0])[k][j][i] = r*cos(xp[3])*sin(th);
    (*x[1])[k][j][i] = r*sin(xp[3])*sin(th);
    (*x[2])[k][j][i] = r*cos(th);
    (*x[3])[k][j][i] = r;
    (*x[4])[k][j][i] = th;
    (*x[5])[k][j][i] = xp[3];
    #endif

    (*x[6])[k][j][i] = xp[1];
    (*x[7])[k][j][i] = xp[2];
    (*x[8])[k][j][i] = xp[3];
  }

  // write to hdf5 
  for (int d = 0; d < NGRIDVARS; d++){
    pack_write_scalar((*x[d]), coordNames[d], OUT_H5_TYPE);
  }

  // determinant and lapse function
  pack_write_axiscalar(G->gdet[CENT], "gdet", OUT_H5_TYPE);
  pack_write_axiscalar(G->lapse[CENT], "lapse", OUT_H5_TYPE);

  // covariant and contravariant metric tensors
  pack_write_Gtensor(G->gcon[CENT], "gcon", OUT_H5_TYPE);
  pack_write_Gtensor(G->gcov[CENT], "gcov", OUT_H5_TYPE);

  // free arrays
  for (int d = 0; d < NGRIDVARS; d++) free(x[d]);

  //close hdf5
  hdf5_close();
}
