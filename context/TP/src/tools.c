// #include "tools.h"
#include <stddef.h>
#include <math.h>

void sort_double(double *restrict a, size_t n)
{
  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = i + 1; j < n; ++j)
    {
      if (a[i] > a[j])
      {
        double tmp = a[i];
        a[i] = a[j];
        a[j] = tmp;
      }
    }
  }
}


double mean_double(double *restrict a, size_t n)
{
  double m = 0.0;

  for (size_t i = 0; i < n; ++i)
    m += a[i];

  m /= (double)n;

  return m;
}


double stddev_double(double *restrict a, size_t n, double m)
{
  double d = 0.0;

  for (size_t i = 0; i < n; ++i)
    d += (a[i] - m) * (a[i] - m);

  d /= (double)(n - 1);

  return sqrt(d);
}


// Couldn't find a way to use the actual BLAS icopy
void my_icopy(int taille, const int* source, int* dest)
{
  for (int i = 0; i < taille; ++i)
  {
    dest[i] = source[i];
  }
}
