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
  double lambda_min = 2-(2 * cos(M_PI * 0 / (*la + 1)));
  double lambda_max = 2-(2 * cos(M_PI * *la / (*la + 1)));
  for (size_t i = 0; i < *la; ++i)  // Peut-être pas utile ?
  {
    double tmp = 2-(2 * cos(M_PI * i / (*la + 1)));
    if (tmp < lambda_min)
      lambda_min = tmp;
    else if (tmp > lambda_max)
      lambda_max = tmp;
  }
  return 2 / (lambda_min + lambda_max);
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
  ==> dgbmv("N", la, lab, kl, ku, -1, AB, la, X, 1, 1, RHS, 1)

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

  // BUG marche pas (stagne à un moment)
  const double norme_b_const = cblas_dnrm2(*la, RHS, 1);
  double residu_b = cblas_dnrm2(*la, RHS, 1) / norme_b_const;
  resvec[*nbite] = residu_b;
  printf("\n %lf \n", residu_b);
  printf("\n tol %lf \n", *tol);
  double* y = malloc(sizeof(double) * *la);
  // cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *lab, *kl, *ku, -1, AB, *la, X, 1, 1, RHS, 1);
  // dgbmv_("N", la, lab, kl, ku, -1, AB, la, X, 1, 1, RHS, 1);
  while (residu_b > *tol)
  {
    // Vecteur = Vecteur + scalaire * ( vecteur - matrice * vecteur)
    //*X = *X + *alpha_rich * (*RHS - *AB * *X);
    cblas_dcopy(*la, RHS, 1, y, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1, AB, *lab, X, 1, 1, y, 1);
    residu_b = cblas_dnrm2(*la, y, 1) / norme_b_const;
    printf("\n %lf \n", residu_b);
    cblas_daxpy(*la, *alpha_rich, y, 1, X, 1);
    // dgbmv_("N", la, lab, kl, ku, -1, AB, la, X, 1, RHS, resvec, 1);  // (*RHS - *AB * *X)
    // daxpy_(la, alpha_rich, RHS, 1, X, 1);  // (*X = *X + *alpha * resultatprécédent)
    ++*nbite;
    if (*nbite == *maxit)
      break;
    resvec[*nbite] = residu_b;
  /*
   calcul residu
   while residu > tol
     if nbite==maxit
       break
     calcul xk
     ++nbite
     calcul residu

   */
  }
  free(y);
}

// A matrice de poisson!
// D diagonale
// E diagonale inférieure
// F diagonale supérieure

// Construit M en extrayant D de A (A et M en stockage bande)*
// TODO
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

// Construit M en extrayant D et E de A (A et M en stockage bande)
// TODO
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  /*
   A = D-E-F
   Jacobi: M = D
   Gauss Seidel = M = D-E
   */
  double norme_b = cblas_dnrm2(*la, RHS, 1);
  resvec[*nbite] = norme_b;
  double* y = malloc(sizeof(double) * *la);
  while ((cblas_dnrm2(*la, RHS, 1) / norme_b) > *tol)
  {
    // Vecteur = Vecteur + Matrice * vecteur~>( vecteur - matrice * vecteur)
    //*X = *X + *MB * (*RHS - *AB * *X);
    cblas_dcopy(*lab, RHS, 1, y, 1);
    //TODO Changer valeurs pour refléter changements faits dans l'autre fonction
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *lab, *kl, *ku, -1, AB, *lab, X, 1, 1, y, 1);
    norme_b = cblas_dnrm2(*la, y, 1);
    // X = X + MB * y
    // TODO utiliser dgbtrs pour avoir MB-1
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *lab, *lab, *kl, *ku, -1, MB, *lab, y, 1, 1, X, 1);
    ++*nbite;
    if (*nbite == *maxit)
      break;
    resvec[*nbite] = norme_b;
  }
  free(y);
}

