/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"

#define TRF 0
#define TRI 1
#define SV 2

// To measure performance
#define MAX_SAMPLES 500
#define n 10   // Size of matrix
#define r 5000 // Number of repetitions


int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info = 1;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;

  double relres;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  NRHS=1;
  nbpoints=n;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
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
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  printf("Solution with LAPACK\n");
  ipiv = (int *) calloc(la, sizeof(int));

  /* LU Factorization */
  if (IMPLEM == TRF) {
    strcpy(title, "TRF");

    // Performance measure
    int* ipiv_copy = (int *) calloc(la, sizeof(int));
    double* AB_copy = (double *) malloc(sizeof(double) * lab * la);
    for (size_t i = 0; i < MAX_SAMPLES; ++i)
    {
      do
      {
        // NOTE: might add a small error on time measure but I don't want to measure copying time
        double total = 0;
        for (size_t j = 0; j < r; ++j)
        {
          // Copy AB and ipiv
          cblas_dcopy(lab * la, AB, 1, AB_copy, 1);
          my_icopy(la, ipiv, ipiv_copy);

          clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
          dgbtrf_(&la, &la, &kl, &ku, AB_copy, &lab, ipiv_copy, &info);
          clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
          total += (double)(t2.tv_nsec - t1.tv_nsec);
        }

        elapsed = total / (double)r;
      }
      while (elapsed <= 0.0);
      samples[i] = elapsed;
    }
    cblas_dcopy(lab * la, AB_copy, 1, AB, 1);
    my_icopy(la, ipiv_copy, ipiv);
    free(ipiv_copy);
    free(AB_copy);

    // dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  }

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  if (IMPLEM == TRI) {
    strcpy(title, "TRI");
    int* ipiv_copy = (int *) calloc(la, sizeof(int));
    double* AB_copy = (double *) malloc(sizeof(double) * lab * la);
    for (size_t i = 0; i < MAX_SAMPLES; ++i)
    {
      do
      {
        double total = 0;
        for (size_t j = 0; j < r; ++j)
        {
          // Copy AB and ipiv
          cblas_dcopy(lab * la, AB, 1, AB_copy, 1);
          my_icopy(la, ipiv, ipiv_copy);

          clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
          dgbtrftridiag(&la, &la, &kl, &ku, AB_copy, &lab, ipiv_copy, &info);
          clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
          total += (double)(t2.tv_nsec - t1.tv_nsec);
        }

        elapsed = total / (double)r;
      }
      while (elapsed <= 0.0);
      samples[i] = elapsed;
    }
    cblas_dcopy(lab * la, AB_copy, 1, AB, 1);
    my_icopy(la, ipiv_copy, ipiv);
    free(ipiv_copy);
    free(AB_copy);

    // dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  }

  if (IMPLEM == TRI || IMPLEM == TRF){
    strcpy(title_two, "TRS");
    /* Solution (Triangular) */
    if (info==0){

      // Performance measure
      double* RHS_copy = calloc(la, sizeof(double));
      for (size_t i = 0; i < MAX_SAMPLES; ++i)
      {
        do
        {
          double total = 0;
          for (size_t j = 0; j < r; ++j)
          {
            cblas_dcopy(la, RHS, 1, RHS_copy, 1);

            clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
            dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS_copy, &la, &info);
            clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
            total += (double)(t2.tv_nsec - t1.tv_nsec);
          }

          elapsed = total / (double)r;
        }
        while (elapsed <= 0.0);
        samples_deux[i] = elapsed;
      }
      cblas_dcopy(la, RHS_copy, 1, RHS, 1);
      // WARNING Invalid pointer when trying to free this
      // free(RHS_copy);

      // dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    }
    else
      {printf("\n INFO = %d\n",info);}
  }

  /* It can also be solved with dgbsv */
  if (IMPLEM == SV) {
    strcpy(title, "SV");
    int* ipiv_copy = (int *) calloc(la, sizeof(int));
    double* AB_copy = (double *) malloc(sizeof(double) * lab * la);
    double* RHS_copy = (double *) malloc(sizeof(double)*la);
    for (size_t i = 0; i < MAX_SAMPLES; ++i)
    {
      do
      {
        double total = 0;
        for (size_t j = 0; j < r; ++j)
        {
          cblas_dcopy(lab * la, AB, 1, AB_copy, 1);
          my_icopy(la, ipiv, ipiv_copy);
          cblas_dcopy(la, RHS, 1, RHS_copy, 1);

          clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
          dgbsv_(&la, &kl, &ku, &NRHS, AB_copy, &lab, ipiv, RHS_copy, &la, &info);
          clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
          total += (double)(t2.tv_nsec - t1.tv_nsec);
        }

        elapsed = total / (double)r;
      }
      while (elapsed <= 0.0);
      samples[i] = elapsed;
    }
    cblas_dcopy(lab * la, AB_copy, 1, AB, 1);
    my_icopy(la, ipiv_copy, ipiv);
    cblas_dcopy(la, RHS_copy, 1, RHS, 1);
    free(RHS_copy);
    free(ipiv_copy);
    free(AB_copy);

    // dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
  }

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  relres = relative_forward_error(RHS, EX_SOL, &la);

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

  if (IMPLEM != SV)
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
  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
