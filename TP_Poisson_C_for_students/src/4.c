/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/

#include "lib_poisson1D.h"
#include "cblas.h"


// chronometre par compteur de tics rdtsc
unsigned long long rdtsc(void)
{
  unsigned long long a, d;
  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
  return (d << 32) | a;
}



int main(int argc,char *argv[]){

	int nbpoints, la;
	int ku, kl, kv, lab;
	
	double it = 1 ; 

	double *x_col, *x_row, *y_col, *y_row ; 
	double *AB_col, *AB_row;

	nbpoints=8;
	la=nbpoints-2;
	
	kv=0;
	ku=1;
	kl=1;
	lab=kv+kl+ku+1;


	x_col=(double *) malloc(sizeof(double)*la);
	x_row=(double *) malloc(sizeof(double)*la);
	y_col=(double *) malloc(sizeof(double)*la);
	y_row=(double *) malloc(sizeof(double)*la);
	
	
	for(int i = 0 ; i <  la; i++){
		x_col[i] = 1.0 ;  
		y_col[i] = 1.0 ; 
	}
	
	
	for(int i = 0 ; i < la ; i++){
		x_row[i] = 1.0 ;  
		y_row[i] = 1.0 ; 
	}		

	AB_col = (double *) malloc(sizeof(double)*lab*la);
	AB_row = (double *) malloc(sizeof(double)*lab*la);


	double alpha ; 
	double beta ; 
	int stride ; 
	alpha = 1.0; 
	beta = 0.0 ;
	stride = 1 ; 

	printf("--------- Poisson COL ---------\n\n");
	set_GB_operator_colMajor_poisson1D(AB_col, &lab, &la, &kv);	
	write_GB_operator_colMajor_poisson1D(AB_col, &lab, &la, "AB_col.dat");	
	
	
	double before = (double)rdtsc();
	for (int i = 0 ; i < it ; i++){
		
		cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, alpha, AB_col, lab, x_col, stride , beta, y_col,  stride);
	}

	double after = (double)rdtsc();
	printf(" %lf\n", (after - before)/it);	

	write_xy(x_col, x_col, &la, "x_col.dat");
	write_xy(y_col, y_col, &la, "y_col.dat");
	
	
	printf("--------- Poisson ROW ---------\n\n");
	set_GB_operator_rowMajor_poisson1D(AB_row, &lab, &la, &kv);	
	write_GB_operator_rowMajor_poisson1D(AB_row, &lab, &la, "AB_row.dat");	
	
	
	before = (double)rdtsc();

	for (int i = 0 ; i < it ; i++){
		cblas_dgbmv(CblasRowMajor, CblasNoTrans, la, la, kl, ku, alpha, AB_row, lab, x_row, stride , beta, y_row,  stride);

	}

	after = (double)rdtsc();
	printf(" %lf\n", (after - before)/it);	

	write_xy(x_row, x_row, &la, "x_row.dat");
	write_xy(y_row, y_row, &la, "y_row.dat");	
	
	
	
	
	
	
	
	free(AB_col);
	free(AB_row);
	free(x_col);
	free(x_row);
	free(y_col);
	free(y_row);

	
}
