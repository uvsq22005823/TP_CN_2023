/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

// AB => Matrice stockée en General Band
// Creating 1D Poisson matrix in AB
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  // lab is row size and la is column size
  // i.e. lab number of columns and la number of rows
	int k_int = (int)*kv;
  int i = 0;
  int j = 0;
	int indice = 0;

  for (int k = 0; k < k_int; ++k)
    {
      AB[indice] = 0;
      ++indice;
    }

  AB[indice] = 0;
  ++indice;
  AB[indice] = -1;
  ++indice;
  AB[indice] = 2;
  ++indice;
  ++i;

  while (i != (int)*la - 1)
  {
    for (int k = 0; k < k_int; ++k)
    {
      AB[indice] = 0;
      ++indice;
    }

    AB[indice] = -1;
    ++indice;
    AB[indice] = 2;
    ++indice;
    AB[indice] = -1;
		++indice;

    ++i;
  }

  AB[indice] = 0;
  ++indice;
  AB[indice] = -1;
  ++indice;
  AB[indice] = 2;
  ++indice;
}


void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  // Comme celle au dessus mais avec des 1 à la place des 2 et 0 à la place de -1
  int k_int = (int)*kv;
  int i = 0;
  int j = 0;
	int indice = 0;

  for (int k = 0; k < k_int; ++k)
    {
      AB[indice] = 0;
      ++indice;
    }

  AB[indice] = 0;
  ++indice;
  AB[indice] = 0;
  ++indice;
  AB[indice] = 1;
  ++indice;
  ++i;

  while (i != (int)*la - 1)
  {
    for (int k = 0; k < k_int; ++k)
    {
      AB[indice] = 0;
      ++indice;
    }

    AB[indice] = 0;
    ++indice;
    AB[indice] = 1;
    ++indice;
    AB[indice] = 0;
		++indice;

    ++i;
  }

  AB[indice] = 0;
  ++indice;
  AB[indice] = 0;
  ++indice;
  AB[indice] = 1;
  ++indice;
}


// Danx Ax = b ; fonction qui permet de créer b ; ici x pas dans espace
// b de la forme (T0, 0, 0, ..., 0, T1) car g=0
// (valeur initiale de b)
// BC0 => T0 ; BC1 => T1
// la => leading dimension
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  RHS[0] = *BC0;
  for (int i=1; i < *la-1; ++i)
  {
    RHS[i] = 0;
  }
  RHS[*la-1] = *BC1;
}


void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
}

void set_grid_points_1D(double* x, int* la){
}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
	fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

int indexABCol(int i, int j, int *lab){
  return 0;
}
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}
