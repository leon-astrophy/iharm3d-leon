//******************************************************************************
//*                                                                            *
//* MAIN.C                                                                     *
//*                                                                            *
//* ESTABLISH MPI COMMUNICATION, LOOP OVER TIME, COMPLETE                      *
//*                                                                            *
//******************************************************************************

// header files
#include "decs.h"
#include "defs.h"

// more header files
#include <time.h>
#include <sys/stat.h>

// positron //
#include "positrons.h"

// main body of the program
int main(int argc, char *argv[])
{

  //mpi stuff
  mpi_initialization(argc, argv);

  //print out license
  if (mpi_io_proc()) {
    fprintf(stdout, "\n          ************************************************************\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *                          IHARM3D                         *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *          Gammie, McKinney & Toth ApJ 509:444, 2003       *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *  B S Prather                                             *\n");
    fprintf(stdout, "          *  G N Wong                                                *\n");
    fprintf(stdout, "          *  B R Ryan                                                *\n");
    fprintf(stdout, "          *  J C Dolence                                             *\n");
    fprintf(stdout, "          *  C F Gammie                                              *\n");
    fprintf(stdout, "          *  S M Ressler                                             *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *                          SYNTAX                          *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *    -p /path/to/param.dat                                 *\n");
    fprintf(stdout, "          *    -o /path/to/output/dir                                *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          ************************************************************\n\n");
  }

  // Read command line arguments, parameter files
  char pfname[STRLEN] = "param.dat";
  char outputdir[STRLEN] = ".";
  for (int n = 0; n < argc; n++) {
    // Check for argv[n] of the form '-*'
    if (*argv[n] == '-' && *(argv[n]+1) != '\0' && *(argv[n]+2) == '\0' &&
        n < argc-1) {
      if (*(argv[n]+1) == 'o') { // Set output directory path
        strcpy(outputdir, argv[++n]);
      }
      if (*(argv[n]+1) == 'p') { // Set parameter file path
        strcpy(pfname, argv[++n]);
      }
    }
  }

  // Read parameter file before we move away from invocation dir
  set_core_params();
  set_problem_params();
  read_params(pfname);

  // Remove 'abort' file if it exists
  char abort_fname[256] = "abort";
  remove(abort_fname);

  // Chdir to output directory and make output folders
  if( chdir(outputdir) != 0) {
    fprintf(stderr, "Output directory does not exist!\n");
    exit(2);
  }

  //mpi stuff
  if (mpi_io_proc()) {
    int is_error = mkdir("dumps/", 0777) || mkdir("restarts/", 0777);
    if (is_error == -1 && errno != EEXIST){
      fprintf(stderr, "Could not make dumps/restarts directory.  Is output dir writeable?\n");
      exit(1);
    }
  }

  //get number of threads
  #pragma omp parallel
  {
    #pragma omp master
    {
      nthreads = omp_get_num_threads();
    }
  }
  omp_set_num_threads(1);

  ///////////////////////////////////////////////////////////////////////
  // TODO centralize more allocations here with safe, aligned _mm_malloc
  ///////////////////////////////////////////////////////////////////////
  // allocate arrays
  struct GridGeom *G = calloc(1,sizeof(struct GridGeom));
  struct FluidState *S = calloc(1,sizeof(struct FluidState));

  // Leon's patch. calculate isco radius here //
  double z1 = 1 + pow(1 - a*a,1./3.)*(pow(1+a,1./3.) + pow(1-a,1./3.));
  double z2 = sqrt(3*a*a + z1*z1);
  R_isco = 3 + z2 - sqrt((3-z1)*(3 + z1 + 2*z2));

  // Perform initializations, either directly or via checkpoint
  is_restart = restart_init(G, S);
  if (!is_restart) {
    
    // Geneate initial conditions
    init(G, S);

  // Leon's patch, initialize positrons //
#if POSITRONS
  init_positrons(G,S);
#endif

    // Set globals
    nstep = 0;
    t = 0;
    dump_cnt = 0;

    // Zero the pflag array
    zero_arrays();
    if (mpi_io_proc())
      fprintf(stdout, "Initial conditions generated\n\n");
  
  }

  // Leon's patch, initialize cooling //
#if COOLING
  init_cooling(G);
#endif  

  // Leon's patch, positron related //
#if POSITRONS
  // Leon's patch, set units only if we are doing pair productions //
  set_units();
#endif

  // In case we're restarting and these changed
  tdump = t + DTd;
  tlog = t + DTl;

  // Initial diagnostics
  diag(G, S, DIAG_INIT);
  if (!is_restart) restart_write(S);

  //mpi stuff
  if (mpi_io_proc())
    fprintf(stdout, "t = %e tf = %e\n", t, tf);

//*******************************************************************************
//*
//*    MAIN LOOP
//*
//*******************************************************************************

  //mpi stuff
  if (mpi_io_proc())
    fprintf(stdout, "\nEntering main loop\n");

  //initialize
  time_init();
  int dumpThisStep = 0;

  //advance in time
  while (t < tf) {

    // Handle abort case
    if ( access(abort_fname, F_OK) != -1 ) {
      if (mpi_io_proc()) {
        fprintf(stdout, "\nFound 'abort' file. Quitting now.\n\n");
      }
      diag(G, S, DIAG_ABORT);
      restart_write_backend(S, IO_ABORT);
      mpi_finalize();
      return 0;
    }

    //set some stuff
    dumpThisStep = 0;
    timer_start(TIMER_ALL);

    // Step variables forward in time
    step(G, S);
    nstep++;

    // Don't step beyond end of run
    if (t + dt > tf) {
      dt = tf - t;
    }

    //mpi stuff
    if (mpi_io_proc()) {
      fprintf(stdout, "t = %10.5g dt = %10.5g n = %8d\n", t, dt, nstep);
    }

    // File I/O with set frequencies
    if (t < tf) {
      if (t >= tdump) {
        dumpThisStep = 1;
        diag(G, S, DIAG_DUMP);
        tdump += DTd;
      }
      if (t >= tlog) {
        diag(G, S, DIAG_LOG);
        tlog += DTl;
      }
      if (nstep % DTr == 0)
        restart_write(S);
    }

    //count time
    timer_stop(TIMER_ALL);

    //report code efficiencies
    if (nstep % DTp == 0)
      report_performance();
//t = tf;
  }
  
//*******************************************************************************
//*
//*    END MAIN LOOP
//*
//*******************************************************************************

  // output diagonstic variables 
  if (dumpThisStep == 0) diag(G, S, DIAG_FINAL);

  // mpi stuff
  mpi_finalize();

  return 0;
}
