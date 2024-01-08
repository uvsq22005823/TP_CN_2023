/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"

#define ALPHA 0
#define JAC 1
#define GS 2

// To measure performance
#define MAX_SAMPLES 500
#define n 12   // Size of matrix
#define r 500 // Number of repetitions

int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, lab, kv;
  int *ipiv;
  int info;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *SOL, *EX_SOL, *X;
  double *AB;
  double *MB;

  double temp, relres;

  double opt_alpha;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  /* Size of the problem */
  NRHS=1;
  nbpoints=n;
  la=nbpoints-2;

  /* Dirichlet Boundary conditions */
  T0=5.0;
  T1=20.0;


  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  SOL=(double *) calloc(la, sizeof(double));
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  /* Setup the Poisson 1D problem */
  /* General Band Storage */
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=0;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  // For time measures
  struct timespec t1, t2;
  double elapsed = 0.0;
  double samples[MAX_SAMPLES];
  double samples_deux[MAX_SAMPLES];


  double size_b = (double)(sizeof(double) * la * lab);;
  double size_kib = size_b / (1024.0);
  double size_mib = size_b / (1024.0 * 1024.0);

  char* title = calloc(4, sizeof(char));
  char* title_two = calloc(4, sizeof(char));

  AB = (double *) malloc(sizeof(double)*lab*la);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

  /* uncomment the following to check matrix A */
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  /********************************************/
  /* Solution (Richardson with optimal alpha) */

  /* Computation of optimum alpha */
  opt_alpha = richardson_alpha_opt(&la);
  printf("Optimal alpha for simple Richardson iteration is : %lf",opt_alpha);

  /* Solve */
  double tol=1e-3;
  int maxit=1000;
  double *resvec;
  int nbite=0;

  resvec=(double *) calloc(maxit, sizeof(double));

  /* Solve with Richardson alpha */
  if (IMPLEM == ALPHA) {

    // Performance measure
    strcpy(title, "ALP");
    double* SOL_copy = (double *) malloc(sizeof(double) * la);
    double* resvec_copy = malloc(maxit * sizeof(double));
    int nbite_copy = nbite;
    for (size_t i = 0; i < MAX_SAMPLES; ++i)
    {
      do
      {
        // NOTE: might add a small error on time measure but I don't want to measure copying time
        double total = 0;
        for (size_t j = 0; j < r; ++j)
        {
          cblas_dcopy(la, SOL, 1, SOL_copy, 1);
          cblas_dcopy(maxit, resvec, 1, resvec_copy, 1);
          nbite_copy = nbite;

          clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
          // NOTE: Add calculation of optimum alpha here ? Or maybe measure the time it takes?
          richardson_alpha(AB, RHS, SOL_copy, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec_copy, &nbite_copy);
          clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
          total += (double)(t2.tv_nsec - t1.tv_nsec);
        }

        elapsed = total / (double)r;
      }
      while (elapsed <= 0.0);
      samples[i] = elapsed;
    }
    cblas_dcopy(la, SOL_copy, 1, SOL, 1);
    cblas_dcopy(maxit, resvec_copy, 1, resvec, 1);
    nbite = nbite_copy;
    free(SOL_copy);
    free(resvec_copy);


    // richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
  }

  /* Richardson General Tridiag */

  // NOTE: I measure perf for matrix extraction but results makes it feel unrelevent

  /* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
  kv = 1;
  ku = 1;
  kl = 1;
  MB = (double *) malloc(sizeof(double)*(lab)*la);
  if (IMPLEM == JAC) {

    // Performance measure
    strcpy(title, "JAC");
    double* MB_copy = (double *) malloc(sizeof(double)*(lab)*la);
    for (size_t i = 0; i < MAX_SAMPLES; ++i)
    {
      do
      {
        // NOTE: might add a small error on time measure but I don't want to measure copying time
        double total = 0;
        for (size_t j = 0; j < r; ++j)
        {
          cblas_dcopy(lab * la, MB, 1, MB_copy, 1);

          clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
          extract_MB_jacobi_tridiag(AB, MB_copy, &lab, &la, &ku, &kl, &kv);
          clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
          total += (double)(t2.tv_nsec - t1.tv_nsec);
        }

        elapsed = total / (double)r;
      }
      while (elapsed <= 0.0);
      samples[i] = elapsed;
    }
    cblas_dcopy(lab * la, MB_copy, 1, MB, 1);
    free(MB_copy);

    // extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);


  } else if (IMPLEM == GS) {

    // Performance measure
    strcpy(title, "GS");
    double* MB_copy = (double *) malloc(sizeof(double)*(lab)*la);
    for (size_t i = 0; i < MAX_SAMPLES; ++i)
    {
      do
      {
        // NOTE: might add a small error on time measure but I don't want to measure copying time
        double total = 0;
        for (size_t j = 0; j < r; ++j)
        {
          cblas_dcopy(lab * la, MB, 1, MB_copy, 1);

          clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
          extract_MB_gauss_seidel_tridiag(AB, MB_copy, &lab, &la, &ku, &kl, &kv);
          clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
          total += (double)(t2.tv_nsec - t1.tv_nsec);
        }

        elapsed = total / (double)r;
      }
      while (elapsed <= 0.0);
      samples[i] = elapsed;
    }
    cblas_dcopy(lab * la, MB_copy, 1, MB, 1);
    free(MB_copy);

    // extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
  }
  /* Solve with General Richardson */
  if (IMPLEM == JAC || IMPLEM == GS) {
    write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "MB.dat");


    // Performance measure
    strcpy(title_two, "RMB");
    double* SOL_copy = malloc(la * sizeof(double));
    double* resvec_copy = malloc(maxit * sizeof(double));
    int nbite_copy = nbite;
    for (size_t i = 0; i < MAX_SAMPLES; ++i)
    {
      do
      {
        // NOTE: might add a small error on time measure but I don't want to measure copying time
        double total = 0;
        for (size_t j = 0; j < r; ++j)
        {
          cblas_dcopy(la, SOL, 1, SOL_copy, 1);
          cblas_dcopy(maxit, resvec, 1, resvec_copy, 1);
          nbite_copy = nbite;

          clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
          richardson_MB(AB, RHS, SOL_copy, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec_copy, &nbite_copy);
          clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
          total += (double)(t2.tv_nsec - t1.tv_nsec);
        }

        elapsed = total / (double)r;
      }
      while (elapsed <= 0.0);
      samples_deux[i] = elapsed;
    }
    nbite = nbite_copy;
    cblas_dcopy(maxit, resvec_copy, 1, resvec, 1);
    cblas_dcopy(la, SOL_copy, 1, SOL, 1);
    free(resvec_copy);
    free(SOL_copy);


    // richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
  }

  /* Write solution */
  write_vec(SOL, &la, "SOL.dat");

  /* Write convergence history */
  write_vec(resvec, &nbite, "RESVEC.dat");

  /* Relative forward error */
  relres = relative_forward_error(SOL, EX_SOL, &la);

  printf("\nThe relative forward error is relres = %e\n",relres);



    // Printing measures
  sort_double(samples, MAX_SAMPLES);
  double min  = samples[0];
  double max  = samples[MAX_SAMPLES - 1];
  double mean = mean_double(samples, MAX_SAMPLES);
  double dev  = stddev_double(samples, MAX_SAMPLES, mean);

  //Size in MiB / time in seconds
  double mbps = size_mib / (mean / 1e9);

  printf("%10s; %15s; %15s; %10s; %10s; %15s; %15s; %15s; %26s; %10s\n",
	 "title",
	 "KiB", "MiB", "n", "r",
   "min", "max", "mean", "stddev (%)", "MiB/s");
  printf("%10s; %15.3lf; %15.3lf; %10d; %10d; %15.3lf; %15.3lf; %15.3lf; %15.3lf (%6.3lf %%); %10.3lf\n",
        title,
        3 * size_kib,
        3 * size_mib,
        n,
        r,
        min,
        max,
        mean,
        dev,
        (dev * 100.0 / mean),
        mbps);

  if (IMPLEM != ALPHA)
  {
    // Printing measures for sample deux
    sort_double(samples_deux, MAX_SAMPLES);
    min  = samples_deux[0];
    max  = samples_deux[MAX_SAMPLES - 1];
    mean = mean_double(samples_deux, MAX_SAMPLES);
    dev  = stddev_double(samples_deux, MAX_SAMPLES, mean);

    //Size in MiB / time in seconds
    mbps = size_mib / (mean / 1e9);

    printf("%10s; %15.3lf; %15.3lf; %10d; %10d; %15.3lf; %15.3lf; %15.3lf; %15.3lf (%6.3lf %%); %10.3lf\n",
          title_two,
          3 * size_kib,
          3 * size_mib,
          n,
          r,
          min,
          max,
          mean,
          dev,
          (dev * 100.0 / mean),
          mbps);
  }




  free(title);
  free(title_two);
  free(resvec);
  free(RHS);
  free(SOL);
  free(EX_SOL);
  free(X);
  free(AB);
  free(MB);
  printf("\n\n--------- End -----------\n");
}
