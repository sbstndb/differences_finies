/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/

#include "lib_poisson1D.h"

unsigned long long rdtsc(void)
{
  unsigned long long a, d;
  
  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
  
  return (d << 32) | a;
}



int main(int argc,char *argv[]){


	int it = 100 ; 

	int ierr;
	int jj;
	int nbpoints, la;
	int ku, kl, kv, lab;
	int *ipiv;
	int info;
	int NRHS;
	double T0, T1;
	double *RHS_col, *RHS_row, *EX_SOL, *X_col, *X_row;
	double *AB_col, *AB_row;
	// kv est a 0 pour le produit matriciel , ou 1 sinon 
	double temp, relres;

	NRHS=1;
	nbpoints=80;
	la=nbpoints-2;
	T0=-5.0;
	T1=5.0;


	RHS_col=(double *) malloc(sizeof(double)*la);
	RHS_row=(double *) malloc(sizeof(double)*la);
	EX_SOL=(double *) malloc(sizeof(double)*la);
	X_col=(double *) malloc(sizeof(double)*la);
	X_row=(double *) malloc(sizeof(double)*la);
	
	set_grid_points_1D(X_col, &la);
	set_grid_points_1D(X_row, &la);
		
	set_dense_RHS_DBC_1D(RHS_col,&la,&T0,&T1);
	set_dense_RHS_DBC_1D(RHS_row,&la,&T0,&T1);
	
	set_analytical_solution_DBC_1D(EX_SOL, X_col, &la, &T0, &T1);	
	set_analytical_solution_DBC_1D(EX_SOL, X_row, &la, &T0, &T1);



	kv=1;
	ku=1;
	kl=1;
	lab=kv+kl+ku+1;

	AB_col = (double *) malloc(sizeof(double)*lab*la);
	AB_row = (double *) malloc(sizeof(double)*lab*la);

	info=0;

	/* working array for pivot used by LU Factorization */
	ipiv = (int *) calloc(la, sizeof(int));

	set_grid_points_1D(X_col, &la);
	set_grid_points_1D(X_row, &la);
	set_dense_RHS_DBC_1D(RHS_col,&la,&T0,&T1);
	set_dense_RHS_DBC_1D(RHS_row,&la,&T0,&T1);
	set_analytical_solution_DBC_1D(EX_SOL, X_col, &la, &T0, &T1);


	// for col
	printf("--------- Poisson COL ---------\n\n");
	set_GB_operator_colMajor_poisson1D(AB_col, &lab, &la, &kv);	

	double before = (double)rdtsc();
	for (int i = 0 ; i < it ; i++){	
		//write_GB_operator_colMajor_poisson1D(AB_col, &lab, &la, "AB_col.dat");
		info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, la, kl, ku, NRHS, AB_col, lab, ipiv, RHS_col, la);
	}
	double after = (double)rdtsc();
	printf(" %lf\n", (after - before)/it);	
	
	write_xy(RHS_col, X_col, &la, "SOL_col.dat");
 	printf("\n INFO DGBSV COL = %d\n",info);

	// for row
	printf("--------- Poisson ROW ---------\n\n");	
	set_GB_operator_rowMajor_poisson1D(AB_row, &lab, &la, &kv);
	//write_GB_operator_rowMajor_poisson1D(AB_row, &lab, &la, "AB_row.dat");
	
	before = (double)rdtsc();
	for (int i = 0 ; i < it ; i++){		
		info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, la, kl, ku, NRHS, AB_row, la, ipiv, RHS_row, NRHS);	 
	}
	after = (double)rdtsc();
	printf(" %lf\n", (after - before)/it);
	
	write_xy(RHS_row, X_row, &la, "SOL_row.dat");
	printf("\n INFO DGBSV ROW = %d\n",info);
}
