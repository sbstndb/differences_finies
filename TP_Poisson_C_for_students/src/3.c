/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/

#include "lib_poisson1D.h"


int size = 5 ;
//int N = 10 ; // number of linear equations
//int KL = 2 ; // number of subdiagonal ; 
//int KU = 2 ; // number of superdiagonals ; 
//int NHRS = 1 ; // number of right hand sides ;


int NRHS=1;
double T0=-5.0;
double T1=5.0;



double *AB_col, *AB_row;

int main(int argc,char *argv[]){

	int kv=2;
	int ku=1;
	int kl=1;
	int lab = kv + ku + kl +1;

	AB_col = calloc(lab *size, lab * size * sizeof(size));
	set_GB_operator_colMajor_poisson1D(AB_col, &lab, &size, &kv);
	write_GB_operator_colMajor_poisson1D(AB_col, &lab, &size, "AB_col.dat");
	
	AB_row = calloc(lab *size, lab * size * sizeof(size));
	set_GB_operator_rowMajor_poisson1D(AB_row, &lab, &size, &kv);
	write_GB_operator_rowMajor_poisson1D(AB_row, &lab, &size, "AB_row.dat");	


}
