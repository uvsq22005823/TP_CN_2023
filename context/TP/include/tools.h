/**************************************************/
/* tools.h                                        */
/* Header for tools used in perf measures         */
/*                                                */
/**************************************************/
#ifndef _TOOLS_H_

#define _TOOLS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void sort_double(double *restrict a, size_t n);
double mean_double(double *restrict a, size_t n);
double stddev_double(double *restrict a, size_t n, double m);
void my_icopy(int taille, const int* source, int* dest);

#endif
