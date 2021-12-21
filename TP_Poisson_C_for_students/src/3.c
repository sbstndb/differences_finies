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

	int it = 1 ; 

	int ierr;
	int jj;
	int nbpoints, la;
	int ku, kl, kv, lab;
	int *ipiv;
	int info;
	int NRHS;
	double T0, T1;
	double *RHS_col, *RHS_row, *EX_SOL, *coord_col, *coord_row, *x_col, *x_row;
	double *AB_col, *AB_row;
	// kv est a 0 pour le produit matriciel , ou 1 sinon 
	double temp, relres_col, relres_row;

	NRHS=1;
	nbpoints=8;
	la=nbpoints-2;
	T0=-5.0;
	T1=5.0;


	RHS_col=(double *) malloc(sizeof(double)*la);
	RHS_row=(double *) malloc(sizeof(double)*la);
	EX_SOL=(double *) malloc(sizeof(double)*la);
	coord_col=(double *) malloc(sizeof(double)*la);
	coord_row=(double *) malloc(sizeof(double)*la);
	x_col=(double *) malloc(sizeof(double)*la);
	x_row=(double *) malloc(sizeof(double)*la);

	
	set_grid_points_1D(coord_col, &la);
	set_grid_points_1D(coord_row, &la);
		
	set_dense_RHS_DBC_1D(RHS_col,&la,&T0,&T1);
	set_dense_RHS_DBC_1D(RHS_row,&la,&T0,&T1);
	
	set_analytical_solution_DBC_1D(EX_SOL, coord_col, &la, &T0, &T1);	
	set_analytical_solution_DBC_1D(EX_SOL, coord_row, &la, &T0, &T1);



	kv=1;
	ku=1;
	kl=1;
	lab=kv+kl+ku+1;

	AB_col = (double *) malloc(sizeof(double)*lab*la);
	AB_row = (double *) malloc(sizeof(double)*lab*la);

	info=0;

	/* working array for pivot used by LU Factorization */
	ipiv = (int *) calloc(la, sizeof(int));

	set_grid_points_1D(coord_col, &la);
	set_grid_points_1D(coord_row, &la);
	set_dense_RHS_DBC_1D(RHS_col,&la,&T0,&T1);
	set_dense_RHS_DBC_1D(RHS_row,&la,&T0,&T1);
	set_dense_RHS_DBC_1D(x_col,&la,&T0,&T1);
	set_dense_RHS_DBC_1D(x_row,&la,&T0,&T1);
	set_analytical_solution_DBC_1D(EX_SOL, coord_col, &la, &T0, &T1);


	// for col
	printf("--------- Poisson COL ---------\n\n");
	set_GB_operator_colMajor_poisson1D(AB_col, &lab, &la, &kv);	

	write_GB_operator_colMajor_poisson1D(AB_col, &lab, &la, "AB_col_3.dat");
	double before = (double)rdtsc();
	for (int i = 0 ; i < it ; i++){	
		info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, la, kl, ku, NRHS, AB_col, lab, ipiv, x_col, la);
	}
	double after = (double)rdtsc();
	printf(" %lf\n", (after - before)/it);	

	write_xy(x_col, coord_col, &la, "SOL_col_3.dat");
 	printf("\n INFO DGBSV COL = %d\n",info);

	// for row
	printf("--------- Poisson ROW ---------\n\n");	
	set_GB_operator_rowMajor_poisson1D(AB_row, &lab, &la, &kv);
	write_GB_operator_rowMajor_poisson1D(AB_row, &lab, &la, "AB_row_3.dat");
	
	before = (double)rdtsc();
	for (int i = 0 ; i < it ; i++){		
		info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, la, kl, ku, NRHS, AB_row, la, ipiv, x_row, NRHS);	 
	}
	after = (double)rdtsc();
	printf(" %lf\n", (after - before)/it);
	
	write_xy(x_row, coord_row, &la, "SOL_row_3.dat");
	printf("\n INFO DGBSV ROW = %d\n",info);
	
	/* Relative residual */
	temp = cblas_ddot(la, x_col, 1, x_col, 1);
	temp = sqrt(temp);
	cblas_daxpy(la, -1.0, x_col, 1, EX_SOL, 1);
	relres_col = cblas_ddot(la, EX_SOL, 1, EX_SOL, 1);
	relres_col = sqrt(relres_col);
	relres_col = relres_col / temp;
	
		
	set_analytical_solution_DBC_1D(EX_SOL, coord_row, &la, &T0, &T1);	
	temp = cblas_ddot(la, x_row, 1, x_row, 1);
	temp = sqrt(temp);
	cblas_daxpy(la, -1.0, x_row, 1, EX_SOL, 1);
	relres_row = cblas_ddot(la, EX_SOL, 1, EX_SOL, 1);
	relres_row = sqrt(relres_row);
	relres_row = relres_row / temp;	
	set_analytical_solution_DBC_1D(EX_SOL, coord_row, &la, &T0, &T1);	
	
	// erreur max
	double tmp_row = 0.0 ; 
	double tmp_col = 0.0 ; 
	double norm_inf = 0.0 ; 
	for (int i = 0 ; i < la  ; i++){
		if (fabs(x_col[i] - EX_SOL[i]) > tmp_col){
			tmp_col = fabs(x_col[i] - EX_SOL[i]);
			
		}
		if (fabs(x_row[i] - EX_SOL[i]) > tmp_row){
			tmp_row = fabs(x_row[i] - EX_SOL[i]);
		}
		if (fabs(EX_SOL[i] * EX_SOL[i]) > norm_inf){
			norm_inf = EX_SOL[i] * EX_SOL[i];
		}
	}

	printf("\nThe relative residual error is relres = %e (col)\n",relres_col);
	printf("\nThe relative residual error is relres = %e (row)\n",relres_row);

	printf("\nThe inf relative residual error is relres = %e (col)\n",(tmp_col/norm_inf));
	printf("\nThe inf relative residual error is relres = %e (row)\n",(tmp_row/norm_inf));
	
	printf("\nThe inf norm is relres = %e (row)\n",norm_inf);


	free(RHS_col);
	free(RHS_row);
	free(EX_SOL);
	free(x_col);
	free(x_row);
	free(coord_col);
	free(coord_row);
	free(AB_col);
	free(AB_row);
	free(ipiv);

	printf("\n\n--------- End -----------\n");	
	
}
