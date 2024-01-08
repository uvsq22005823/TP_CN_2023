/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <math.h>

/*
 * Richardson
 * x0
 * r0 = b - Ax0
 * while ||rk+1|| > E (E très petit ; rk+1 -> (b - A * xk+1))
 *    xk+1 = xk + M-1 (b - A * xk)
 *
 * M-1 est, en fait, alpha (pour richardson_alpha)
*/
void eig_poisson1D(double* eigval, int *la){
}

double eigmax_poisson1D(int *la){
  return 0;
}

double eigmin_poisson1D(int *la){
  return 0;
}

// Calcule le alpha optimal
// TODO: démontrer formule dans le rapport
// On sait qu'il vaut 2 / (lambdamin+lambdamax) (avec lamda valeur propre de A)
double richardson_alpha_opt(int *la){
  double h = 1.0 / ((*la) +1.0);
  double lambda_min = 4*sin((M_PI * h)/2)*sin(( M_PI * h)/2);
  double lambda_max = 4*sin(((*la) * M_PI * h)/2)*sin(((*la) * M_PI * h)/2);
  return 2 / (lambda_min + lambda_max);
  // return 0.5;
}

// Calcule Richardson
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  /* AB => A (matrice)
   * RHS => b (vecteur)
   * X => x (vecteur)
   * alpha_rich => alpha
   * lab => number of columns
   * la => number of rows
   * ku => sur diags
   * kl => sous diags
   * tol => epsilon I guess
   * maxit => nb itérations max
   * resvec => vecteur contenant les résidus ? (vecteur)
   * nbite => nb itérations actuel ??
   */

  /*
  xk+1 = xk + M-1 (b - A * xk)
  dgbmv("N", la, lab, kl, ku, alpha, a, la, x, incx, beta, y, incy)
    => y = alpha*A*x + beta*y

  xk+1 = xk + M-1 * c
  c <= b - A * xk
    => -(A*xk) + b

    => AB -> A -> a
    => X -> xk -> x
    => RHS -> b -> y
    ==> alpha = -1
    ==> beta = 1
  ==> dgbmv("N", la, la, kl, ku, -1, AB, lab, X, 1, 1, RHS, 1)

  Calculer residu entre les deux opérations

    xk + M-1 * c    vecteur + (ici scalaire) * vecteur

    daxpy(n, da, dx, incx, dy, incy)

    X -> xk => dy
    alpha_rich -> M-1 => da
    c => dx ~> RHS
    ==> incx = 1
    ==> incy = 1
    ==> n = la
  ==> daxpy(la, alpha_rich, X, 1, RHS, 1)
  */

  const double norme_b_const = cblas_dnrm2(*la, RHS, 1);
  double residu_b = cblas_dnrm2(*la, RHS, 1) / norme_b_const;

  double* y = malloc(sizeof(double) * *la);

  while (residu_b > *tol)
  {
    cblas_dcopy(*la, RHS, 1, y, 1);

    // y = alpha*A*x + beta*y
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, y, 1);

    // calcul residu
    residu_b = cblas_dnrm2(*la, y, 1) / norme_b_const;
    resvec[*nbite] = residu_b;

    // Vecteur = Vecteur + scalaire * ( vecteur - matrice * vecteur)
    //*X = *X + *alpha_rich * (*RHS - *AB * *X);
    cblas_daxpy(*la, *alpha_rich, y, 1, X, 1);

    ++*nbite;
    if (*nbite == *maxit)
      break;
  }
  free(y);
}

// A matrice de poisson!  => AB
// D diagonale
// E diagonale inférieure
// F diagonale supérieure

// Construit M en extrayant D de A (A et M en stockage bande)
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  for (int i = 0; i < *la; i++)
    MB[i * (*lab) + 1] = AB[i * (*lab) + 1];
}

// Construit M en extrayant D et E de A (A et M en stockage bande)
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  for(int i = 0; i < *la; i++)
  {
    MB[*lab * i + 1] = AB[*lab * i + 1]; // M = D
    MB[*lab * i + 2] = AB[*lab * i + 2]; // M = M - E
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  /*
   A = D-E-F
   Jacobi: M = D
   Gauss Seidel = M = D-E
   */
  const double norme_b_const = cblas_dnrm2(*la, RHS, 1);
  double residu_b = cblas_dnrm2(*la, RHS, 1) / norme_b_const;

  double* y = malloc(sizeof(double) * *la);

  // Declarations for dgbtrs
  int* pivot = calloc(*la, sizeof(int));  // ipiv in dgbtrs ; array pivot indices
  int info = 0;  // Should always remain 0
  const int NRHS = 1;  // Number of columns of y
  const int ku_mb = 0;  // Pas de sur diagonale ici

  while (residu_b > *tol)
  {
    cblas_dcopy(*la, RHS, 1, y, 1);

    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, y, 1);

    residu_b = cblas_dnrm2(*la, y, 1) / norme_b_const;
    resvec[*nbite] = residu_b;

    LAPACK_dgbsv(la, kl, &ku_mb, &NRHS, MB, lab, pivot, y, la, &info);
    if (info != 0)
    {
      //NOTE: crash here?
      fprintf(stderr, "\nIllegal or null value in dgbsv call!\n");
      // break;
    }

    // Vecteur = Vecteur + scalaire * ( vecteur - matrice * vecteur)
    //*X = *X + *alpha_rich * (*RHS - *AB * *X);
    cblas_daxpy(*la, 1, y, 1, X, 1);

    ++*nbite;
    if (*nbite == *maxit)
      break;
  }
  free(y);
  free(pivot);
}
