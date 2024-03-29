/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"


void set_GB_operator_rowMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
	int ii, jj, kk;
	for (kk = 0 ; kk < (*kv); kk++){
		for ( jj = 0 ; jj < (*la) ; jj++){
			AB[jj + kk * (*la)] = 0 ; 
		}
	}
	for (kk = 0 ; kk < (*la) ; kk++){
		AB[ (*kv)*(*la) + kk ] = -1 ;
		AB[ (*kv)*(*la) + kk + 1 * (*la)] = 2 ;
		AB[ (*kv)*(*la) + kk + 2 * (*la) ] = -1 ;			 		
	}
  AB[0] = 0.0;
  if (*kv == 1) { AB[(*la)] = 0 ;}	
  AB[ (*lab)*(*la) - 1 ] = 0.0;
}
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj = 0 ; jj < (*la) ; jj++ ){
    kk = jj * ( *lab );
    if (*kv >= 0){
      for (ii = 0 ; ii < *kv ; ii++){
	AB[kk + ii] = 0.0;
      }
    }
    AB[kk + *kv] = -1.0;
    AB[kk + *kv+1] = 2.0;
    AB[kk + *kv+2] = -1.0;
  }
  AB[0] = 0.0;
  if (*kv >= 1) { AB[*kv] = 0 ;}
  
  AB[ (*lab)*(*la) - 1 ] = 0.0;
}

void set_GB_operator_rowMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
	int ii, jj, kk;
	for (kk = 0 ; kk < (*kv); kk++){
		for ( jj = 0 ; jj < (*la) ; jj++){
			AB[jj + kk * (*la)] = 0 ; 
		}
	}
	for (kk = 0 ; kk < (*la) ; kk++){
		AB[ (*kv)*(*la) + kk ] = 0 ;
		AB[ (*kv)*(*la) + kk + 1 * (*la)] = 1 ;
		AB[ (*kv)*(*la) + kk + 2 * (*la) ] = 0 ;			 		
	}
  AB[0] = 0.0;
  if (*kv == 1) { AB[(*la)] = 0 ;}	
  AB[ (*lab)*(*la) - 1 ] = 0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj = 0 ; jj < (*la) ; jj++){
    kk = jj * (*lab);
    if (*kv >= 0){
      for (ii = 0 ; ii < *kv ; ii++){
	AB[kk + ii] = 0.0;
      }
    }
    AB[kk + *kv] = 0.0;
    AB[kk + *kv+1] = 1.0;
    AB[kk + *kv+2] = 0.0;
  }
  AB[1] = 1.0; // !!! erreur 
  AB[(*lab) * (*la) - 1] = 0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0] = *BC0;
  RHS[(*la)-1] = *BC1;
  for (jj = 1 ; jj < (*la)-1 ; jj++){
    RHS[jj] = 0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T = (*BC1) - (*BC0);
  for ( jj = 0 ; jj < (*la) ; jj++ ){
    EX_SOL[jj] = (*BC0) + X[jj] * DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h = 1.0 / (1.0 * ((*la) + 1));
  for ( jj = 0 ; jj < (*la) ; jj++ ){
    x[jj]=(jj+1) * h;
  }
}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for ( ii = 0 ; ii < (*lab) ; ii++){
      for ( jj = 0 ; jj < (*la) ; jj++){
	fprintf(file,"%lf\t",AB[ii * (*la) + jj]);
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
    for ( ii = 0 ; ii < (*la) ; ii++){
      for ( jj = 0 ; jj < (*lab) ; jj++){
	fprintf(file,"%lf\t",AB[ii * (*lab) + jj]);
      }
      fprintf(file,"\n");
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
    for ( jj = 0 ; jj < (*la) ; jj++){
      fprintf(file,"%lf\n", vec[jj] );
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
    for ( jj = 0 ; jj < (*la) ; jj++){
      fprintf(file,"%lf\t%lf\n", x[jj], vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void eig_poisson1D(double* eigval, int *la){
  int ii;
  double scal;
  for ( ii = 0 ; ii< *la; ii++){
    scal = (1.0 * ii + 1.0) * M_PI_2*(1.0 / (*la+1));
    eigval[ii] = sin(scal);
    eigval[ii] = 4 * eigval[ii] * eigval[ii];
  } 
}

double eigmax_poisson1D(int *la){
  double eigmax;
  eigmax = sin(*la *M_PI_2 * (1.0/(*la+1)));
  eigmax = 4 * eigmax * eigmax;
  return eigmax;
}

double eigmin_poisson1D(int *la){
  double eigmin;
  eigmin = sin(M_PI_2 * (1.0/(*la+1)));
  eigmin = 4 * eigmin * eigmin;
  return eigmin;
}



void LU_tri_col(double* A, double* L, double* U, int* la){
	// A est matrice bande triangulaire colonne
	// L et U reprennent alors le meme format 
	/*	
	function [L,U] = lutrig(A)
	[n, n2] = size(A); 
	L = eye(n,n);
	U = zeros(n,n)
	U(1,1) = A(1,1) ; 
	L(2,1) = A(2,1) / A(1,1) ;
	U(1,2) = A(1,2);
	for i = [2:n-1]
	    U(i,i+1)=A(i,i+1);
	    U(i,i)=A(i,i)-A(i-1,i)*L(i,i-1);
	    L(i+1,i)=A(i+1,i)/U(i,i);
	end
	U(n,n)=A(n,n)-A(n-1,n)*L(n,n-1);
	end
	*/
	
	for (int i = 0 ; i < (*la) ; i ++){
		L[1 + i*3] = 1.0 ; 
	}	
	U[1] = A[1] ; 
	U[2] = A[2] ; 
	L[3] = A[3] / A[1] ; 
	for (int i = 1 ; i < (*la) - 1 ; i++){
		U[2 + 3 * i] = A[2 + 3 * i] ; 
		U[1 + 3 * i] = A[1 + 3 * i] - A[-1 + 3 * i]*L[3 * i] ; 
		L[3 + 3 * i] = A[3 + 3 * i]/U[1 + 3 * i] ; 
	}	 
	U[1 + 3 * ((*la) - 1 )] = A[1 + 3 * ((*la) - 1 )] - A[-1 + 3 * ((*la) - 1 )]*L[3 * ((*la) - 1 )] ; 
	
}

void allocation(double** x_col, double** x_row, double** y_col, double** y_row, double** AB_col, double** AB_row, int la,int lab){
	*x_col = malloc(sizeof(double)*la);
	*x_row = malloc(sizeof(double)*la);
	*y_col = malloc(sizeof(double)*la);
	*y_row = malloc(sizeof(double)*la);
	*AB_col = malloc(sizeof(double)*lab*la);
	*AB_row = malloc(sizeof(double)*lab*la);
}
void deallocation(double** x_col, double** x_row, double** y_col, double** y_row, double** AB_col, double** AB_row){
	free(*AB_col);
	free(*AB_row);
	free(*x_col);
	free(*x_row);
	free(*y_col);
	free(*y_row);
}




double richardson_alpha_opt(int *la){
  //TODO
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit){
  //TODO
}
