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

	double *AB, *L, *U;

	nbpoints=8;
	la=nbpoints-2;
	
	kv=0;
	ku=1;
	kl=1;
	lab=kv+kl+ku+1;


	AB = (double *) malloc(sizeof(double)*lab*la);
	L = (double *) malloc(sizeof(double)*lab*la);
	U = (double *) malloc(sizeof(double)*lab*la); 

	for (int i = 0 ; i < la * lab ; i++){
		AB[i] = 0 ; 
		L[i] = 0 ; 
		U[i] = 0 ;
	}


	printf("--------- Poisson COL ---------\n\n");
	set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);	

	
	double before = (double)rdtsc();
	
	for (int i = 0 ; i < it ; i++){
		LU_tri_col(AB, L, U, &la);
	}

	double after = (double)rdtsc();
	printf(" %lf\n", (after - before)/it);	
	
	write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");		
	write_GB_operator_colMajor_poisson1D(L, &lab, &la, "L.dat");	
	write_GB_operator_colMajor_poisson1D(U, &lab, &la, "U.dat");	
	free(AB);
	free(L); 
	free(U);
	
}
