/******************************************/
/* Exercice 4  poisson1D_direct.c         */
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
	double tmp_x, tmp_y;

	double *x_col, *x_row, *y_col, *y_row , *y_col_init, *y_row_init; 
	double *AB_col, *AB_row;

	nbpoints=8;
	la=nbpoints-2;
	printf("Dimension : %d\n",la);
	kv=0;
	ku=1;
	kl=1;
	lab=kv+kl+ku+1;


	//////////////////////////////////////////
	// premier cas de verification : alpha = 0
	// vecteurs x et y aleatoires
	// AB initialise aleatoire
	// On doit avoir y = beta * y_old
	//////////////////////////////////////////
	printf("\n--------- Exercice 4 - verification 1  ---------\n");	
	allocation(&x_col, &x_row,  &y_col,  &y_row,  &AB_col,  &AB_row,  la, lab) ; 
	
	y_row_init = malloc(sizeof(double)*la);
	y_col_init = malloc(sizeof(double)*la);

	for(int i = 0 ; i <  la; i++){
		tmp_x = 0; 
		tmp_y = (double)rand() / (double)RAND_MAX; 
		x_col[i] = tmp_x ; 
		y_col[i] = tmp_y ; 
		x_row[i] = tmp_x ;
		y_row[i] = tmp_y ;
		y_col_init[i] = tmp_x ;
		y_row_init[i] = tmp_y ;
		for (int j = 0 ; j < lab ; j++){
			AB_row[i * lab + j] = tmp_y; 
			AB_col[i * lab + j] = tmp_y ; 

		}		
	}
		
	write_xy(y_col, y_col, &la, "y_col_init_4_1.dat");	
	write_xy(y_col, y_row, &la, "y_row_init_4_1.dat");	
	
	double alpha ; 
	double beta ; 
	int stride ; 
	alpha = 0.0; 
	beta = 1.8765432 ;
	stride = 1 ; 

	
	cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, alpha, AB_col, lab, x_col, stride , beta, y_col,  stride);


	cblas_dgbmv(CblasRowMajor, CblasNoTrans, la, la, kl, ku, alpha, AB_row, lab, x_row, stride , beta, y_row,  stride);
	
	double res_x_col ; 
	double res_x_row ; 
	double res_x_ana ; 
	
	double res_y_col ; 
	double res_y_row ; 
	double res_y_ana ;
	
	int verification_col = 1 ;
	int verification_row = 1 ;  
		
	for (int i = 0 ; i < la ; i++){
		res_x_col = x_col[i] ; 
		res_y_col = y_col[i] ;
		res_x_row = x_row[i] ; 
		res_y_row = y_row[i] ;
		res_x_ana = 0.0 ; 
		res_y_ana = y_col_init[i] * beta ; 		
		
		if ((fabs(res_x_col - res_x_ana) < 1.0e-10) && (fabs(res_y_col - res_y_ana) < 1.0e-10)){
			verification_col = 0;			
		}	
		if ((fabs(res_x_row - res_x_ana) < 1.0e-10) && (fabs(res_y_row - res_y_ana) < 1.0e-10)){
			verification_row = 0;		
		}
	}
	
	
	if (verification_col == 1){
		printf("Resultat  COL 1 OK \n");
	}
	else{
		printf("Resultat  COL 1 FAUX \n");	
	}
	if (verification_row == 1){
		printf("Resultat  ROW 1 OK \n");
	}
	else{
		printf("Resultat  ROW 1 FAUX \n");	
	}

	

	write_GB_operator_colMajor_poisson1D(AB_col, &lab, &la, "AB_col_4_1.dat");	
	write_GB_operator_rowMajor_poisson1D(AB_row, &lab, &la, "AB_row_4_1.dat");		
	write_xy(x_col, x_col, &la, "x_col_4_1.dat");
	write_xy(y_col, y_col, &la, "y_col_4_1.dat");
	write_xy(x_row, x_row, &la, "x_row_4_1.dat");
	write_xy(y_row, y_row, &la, "y_row_4_1.dat");	
	
	deallocation(&x_col, &x_row,  &y_col,  &y_row,  &AB_col,  &AB_row) ; 

	

	//////////////////////////////////////////
	// second cas de verification : beta = 0
	// vecteurs x et y initialisees nuls sauf sur une valeur precise
	// AB initialise nul sauf sur une valeur precise
	// tel que Ax = v avec v nul sauf sur une valeur precise
	//////////////////////////////////////////
	printf("\n--------- Exercice 4 - verification 2  ---------\n");	
	allocation(&x_col, &x_row,  &y_col,  &y_row,  &AB_col,  &AB_row,  la, lab) ; 

	for(int i = 0 ; i <  la; i++){
		tmp_x = 0; 
		tmp_y = 0; 
		x_col[i] = tmp_x ; 
		y_col[i] = tmp_y ; 
		x_row[i] = tmp_x ;
		y_row[i] = tmp_y ;
		for (int j = 0 ; j < lab ; j++){
			AB_row[i * lab + j] = 0.0 ; 
			AB_col[i * lab + j] = 0.0 ; 

		}		
	}
	
	int position = (int)(la/2);

	float value = 5.54345 ; 
	AB_row[2 * la + position] = value ; 
	AB_col[3 * position] = value ; 
	x_col[position] = value ; 
	x_row[position] = value ; 
	
	write_xy(y_col, y_col, &la, "y_col_init_4_2.dat");	
	write_xy(y_col, y_row, &la, "y_row_init_4_2.dat");	
	
 
	alpha = 1.543; 
	beta = 0.0 ;


	
	cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, alpha, AB_col, lab, x_col, stride , beta, y_col,  stride);


	cblas_dgbmv(CblasRowMajor, CblasNoTrans, la, la, kl, ku, alpha, AB_row, lab, x_row, stride , beta, y_row,  stride);
	
	double res_col = y_col[position-1] ; 
	double res_row = y_row[position-1] ; 
	double res_ana = alpha * value * value ; 
	
	printf("Valeur attendue : %lf\n",res_ana) ;		 
	
	if ( fabs(res_col - res_ana) < 1.0e-10){
		printf("Resultat  COL 2 OK \n");			
	}
	else{
		printf("Resultat  COL 2 FAUX \n");
	}	
	if (fabs(res_row - res_ana) < 1.0e-10){
		printf("Resultat  ROW 2 OK \n"); 	
	}
	else{
		printf("Resultat  ROW 2 FAUX \n");	
	}
	
	
	write_GB_operator_colMajor_poisson1D(AB_col, &lab, &la, "AB_col_4_2.dat");	
	write_GB_operator_rowMajor_poisson1D(AB_row, &lab, &la, "AB_row_4_2.dat");		
	write_xy(x_col, x_col, &la, "x_col_4_2.dat");
	write_xy(y_col, y_col, &la, "y_col_4_2.dat");
	write_xy(x_row, x_row, &la, "x_row_4_2.dat");
	write_xy(y_row, y_row, &la, "y_row_4_2.dat");	
	
	deallocation(&x_col, &x_row,  &y_col,  &y_row,  &AB_col,  &AB_row) ; 
	

	//////////////////////////////////////////
	// troisieme cas de verification : beta = 0
	// vecteurs x et y aleatoires
	// AB identite
	// On doit avoir y = alpha * x
	//////////////////////////////////////////
	printf("\n--------- Exercice 4 - verification 3  ---------\n");	
	allocation(&x_col, &x_row,  &y_col,  &y_row,  &AB_col,  &AB_row,  la, lab) ; 
	
	y_row_init = malloc(sizeof(double)*la);
	y_col_init = malloc(sizeof(double)*la);

	for(int i = 0 ; i <  la; i++){
		tmp_x = (double)rand() / (double)RAND_MAX; 
		tmp_y = (double)rand() / (double)RAND_MAX; 
		x_col[i] = tmp_x ; 
		y_col[i] = tmp_y ; 
		x_row[i] = tmp_x ;
		y_row[i] = tmp_y ;
		y_col_init[i] = tmp_y ;
		y_row_init[i] = tmp_y ;		
	}

	set_GB_operator_colMajor_poisson1D_Id(AB_col, &lab, &la, &kv);
	set_GB_operator_rowMajor_poisson1D_Id(AB_row, &lab, &la, &kv);
		
	write_xy(y_col, y_col, &la, "y_col_init_4_3.dat");	
	write_xy(y_col, y_row, &la, "y_row_init_4_3.dat");	
	

	alpha = 2.43234; 
	beta = 0.0 ;
	stride = 1 ; 

	
	cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, alpha, AB_col, lab, x_col, stride , beta, y_col,  stride);


	cblas_dgbmv(CblasRowMajor, CblasNoTrans, la, la, kl, ku, alpha, AB_row, lab, x_row, stride , beta, y_row,  stride);
	

	
	verification_col = 1 ;
	verification_row = 1 ;  
		
	for (int i = 0 ; i < la ; i++){
		res_y_col = y_col[i] ;
		res_y_row = y_row[i]  ;
		res_y_ana = alpha * x_col[i]  ; 		
		
		if (fabs(res_y_col - res_y_ana) > 1.0e-10){
			verification_col = 0;		
		}	
		if (fabs(res_y_row - res_y_ana) > 1.0e-10){
			verification_row = 0;		
		}
	}
	
	
	if (verification_col == 1){
		printf("Resultat  COL 1 OK \n");
	}
	else{
		printf("Resultat  COL 1 FAUX \n");	
	}
	if (verification_row == 1){
		printf("Resultat  ROW 1 OK \n");
	}
	else{
		printf("Resultat  ROW 1 FAUX \n");	
	}

	

	write_GB_operator_colMajor_poisson1D(AB_col, &lab, &la, "AB_col_4_3.dat");	
	write_GB_operator_rowMajor_poisson1D(AB_row, &lab, &la, "AB_row_4_3.dat");		
	write_xy(x_col, x_col, &la, "x_col_4_3.dat");
	write_xy(y_col, y_col, &la, "y_col_4_3.dat");
	write_xy(x_row, x_row, &la, "x_row_4_3.dat");
	write_xy(y_row, y_row, &la, "y_row_4_3.dat");	
	
	deallocation(&x_col, &x_row,  &y_col,  &y_row,  &AB_col,  &AB_row) ; 




//////////////////////////////////////////
	// quatrieme cas de verification : beta et alpha non nuls
	// vecteurs x et y aleatoires
	// AB identite
	// On doit avoir y = beta y_old + alpha * x
	//////////////////////////////////////////
	printf("\n--------- Exercice 4 - verification 4  ---------\n");	
	allocation(&x_col, &x_row,  &y_col,  &y_row,  &AB_col,  &AB_row,  la, lab) ; 
	
	y_row_init = malloc(sizeof(double)*la);
	y_col_init = malloc(sizeof(double)*la);

	for(int i = 0 ; i <  la; i++){
		tmp_x = (double)rand() / (double)RAND_MAX; 
		tmp_y = (double)rand() / (double)RAND_MAX; 
		x_col[i] = tmp_x ; 
		y_col[i] = tmp_y ; 
		x_row[i] = tmp_x ;
		y_row[i] = tmp_y ;
		y_col_init[i] = tmp_y ;
		y_row_init[i] = tmp_y ;		
	}

	set_GB_operator_colMajor_poisson1D_Id(AB_col, &lab, &la, &kv);
	set_GB_operator_rowMajor_poisson1D_Id(AB_row, &lab, &la, &kv);
		
	write_xy(y_col, y_col, &la, "y_col_init_4_4.dat");	
	write_xy(y_col, y_row, &la, "y_row_init_4_4.dat");	
	

	alpha = 2.43234; 
	beta = 1.43 ;
	stride = 1 ; 

	
	cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, alpha, AB_col, lab, x_col, stride , beta, y_col,  stride);


	cblas_dgbmv(CblasRowMajor, CblasNoTrans, la, la, kl, ku, alpha, AB_row, lab, x_row, stride , beta, y_row,  stride);
	

	
	verification_col = 1 ;
	verification_row = 1 ;  
		
	for (int i = 0 ; i < la ; i++){
		res_y_col = y_col[i] ;
		res_y_row = y_row[i]  ;
		res_y_ana = alpha * x_col[i] + beta * y_col_init[i]  ; 		
		
		if (fabs(res_y_col - res_y_ana) > 1.0e-10){
			verification_col = 0;		
		}	
		if (fabs(res_y_row - res_y_ana) > 1.0e-10){
			verification_row = 0;		
		}
	}
	
	
	if (verification_col == 1){
		printf("Resultat  COL 1 OK \n");
	}
	else{
		printf("Resultat  COL 1 FAUX \n");	
	}
	if (verification_row == 1){
		printf("Resultat  ROW 1 OK \n");
	}
	else{
		printf("Resultat  ROW 1 FAUX \n");	
	}

	

	write_GB_operator_colMajor_poisson1D(AB_col, &lab, &la, "AB_col_4_4.dat");	
	write_GB_operator_rowMajor_poisson1D(AB_row, &lab, &la, "AB_row_4_4.dat");		
	write_xy(x_col, x_col, &la, "x_col_4_4.dat");
	write_xy(y_col, y_col, &la, "y_col_4_4.dat");
	write_xy(x_row, x_row, &la, "x_row_4_4.dat");
	write_xy(y_row, y_row, &la, "y_row_4_4.dat");	
	
	deallocation(&x_col, &x_row,  &y_col,  &y_row,  &AB_col,  &AB_row) ; 




	
}
