/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"
#include <omp.h>

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double *RHS, *RHS2, *EX_SOL, *X;  // Vecteurs
  double **AAB;
  double *AB;  // Matrice

  double temp, relres;

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  RHS2=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
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

  AB = (double *) malloc(sizeof(double)*lab*la);

  // set_GB_operator_colMajor_poisson1D_Id(AB, &lab, &la, &kv);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

  // print matrice vérifier c'est ok
  // int indice = 0;
  // for (size_t i = 0; i < la; ++i) {
  //   for (size_t j = 0; j < lab; ++j)
  //   {
  //     printf("%f ", AB[indice]);
  //     ++indice;
  //   }
  //   printf("\n");
  // }


  // write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  printf("Solution with LAPACK\n");
  /* LU Factorization */
  info=0;
  ipiv = (int *) calloc(la, sizeof(int));
  dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  // ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

  // write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");


  /*
#define LAPACK_dgbtrs_base LAPACK_GLOBAL(dgbtrs,DGBTRS)
void LAPACK_dgbtrs_base(
    char const* trans,
    lapack_int const* n, lapack_int const* kl, lapack_int const* ku, lapack_int const* nrhs,
    double const* AB, lapack_int const* ldab, lapack_int const* ipiv,
    double* B, lapack_int const* ldb,
    lapack_int* info
#ifdef LAPACK_FORTRAN_STRLEN_END
    , size_t
#endif
);
Bon ben si jamais LAPACK_FORTRAN_STRLEN_END est défini
   */


  // CBLAS DGBMV
  cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1.0, AB+1, lab, EX_SOL, 1, 0.0, RHS2, 1);
  // Bon l'appel compile et s'exécute c'est déjà ça j'ai envie de dire
  // ku = 1 => alpha = 1
  // beta = 0
  // => AB+1 * lab

  // Validation (4.3)
  double norme_RHS = cblas_dnrm2(la, RHS, 1);
  cblas_daxpy(la, -1, RHS2, 1, RHS, 1);  // DAXPY c'est Double A * X Plus Y
  double norm_mv = cblas_dnrm2(la, RHS, 1);
  norm_mv /= norme_RHS;  // Should give ~1.0

  // printf("%f\n", norm_mv);








  // ex 5.1
  double t1 = 0.0, t2 = 0.0;

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  t1 = omp_get_wtime();
  dgbtrf_(&la, &lab, &kl, &ku, AB, &lab, ipiv, &info); // Crée une facto LU dans A avec pivot partiel (bloc BLAS3)
  // nb ligne, nb colonne, nb sous-diag, nb sur-diag, leading dim, indice pivot??? (est un truc out anyway donc osef), info(littéralement une info, osef)
  // utiliser LAPACK_dgbtrs  // Solves A * X = B avec A facto LU faites par DGBTRF
  int n = 1;
  LAPACK_dgbtrs("N", &n, &kl, &ku, &lab, AB, &lab, ipiv, RHS, &lab, &info);
  t2 = omp_get_wtime();
  double elapsed_s = (double)(t2 - t1);

  printf("elpased time %f\n", elapsed_s);
  // TODO: étudier perf







  /* Solution (Triangular) */
  if (info==0){
    LAPACK_dgbtrs("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);  // Remplacé dgbtrs_ par LAPACK_dgbtrs
    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
  }else{
    printf("\n INFO = %d\n",info);
  }

  /* It can also be solved with dgbsv */
  // TODO : use dgbsv

  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  // TODO : Compute relative norm of the residual

  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
