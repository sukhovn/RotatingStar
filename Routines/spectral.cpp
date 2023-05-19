#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <spectral.h>
#include <fftw3.h>

#define ODDM(a) ((a)%2 ? -1.0 : 1.0)

//Chebyshev spectral methods functions

double chpoint(int n, int i){//Gives out Chebyshev-Gauss-Lobatto point
	return cos(M_PI * i/(1.0 * n));
}

void charray(int n, double *vec){//Fills out an array with Chebyshev-Gauss-Lobatto points
	for(int i = 0; i <= n; i++)
		vec[i] = cos(M_PI * i/(1.0 * n));

	return;
}

double ch_x(double xst, double xend, int n, int i){//Gives out Chebyshev-Gauss-Lobatto point
	return xst + 0.5 * (xend - xst) * (cos(M_PI * (n - i)/(1.0 * n)) + 1.0);
}

void ch_x_array(double xst, double xend, int n, double *vec){//Fills out an array with Chebyshev-Gauss-Lobatto points
	for(int i = 0; i <= n; i++)
		vec[i] = xst + 0.5 * (xend - xst) * (cos(M_PI * (n - i)/(1.0 * n)) + 1.0);

	return;
}

double chpoint_st(int n, int i){
	return cos(M_PI * (1.0 * i + 0.5)/(1.0 * n));
}

void charray_st(int n, double *vec){
	for(int i = 0; i < n; i++)
		vec[i] = cos(M_PI * (1.0 * i + 0.5)/(1.0 * n));

	return;
}

double ch_x_st(double xst, double xend, int n, int i){
	return xst + 0.5 * (xend - xst) * (cos(M_PI * (1.0 * (n - i) - 0.5)/(1.0 * n)) + 1.0);
}

void ch_x_array_st(double xst, double xend, int n, double *vec){
	for(int i = 0; i < n; i++)
		vec[i] = xst + 0.5 * (xend - xst) * (cos(M_PI * (1.0 * (n - i) - 0.5)/(1.0 * n)) + 1.0);

	return;
}

int ch_x_inv_low(double xst, double xend, int n, double x){
	double tmp = (2.0 * x - xst - xend)/(xend - xst);
	if(tmp >= 1.0) tmp = 0.0;
	else if(tmp <= -1.0) tmp = 1.0;
	else tmp = acos(tmp)/M_PI;
	return floor(n*(1.0 - tmp));
}

int ch_x_inv_high(double xst, double xend, int n, double x){
	double tmp = n*(1.0 - acos((2.0 * x - xst - xend)/(xend - xst))/M_PI);
	if(tmp == (int) tmp) return tmp++;
	else return ceil(tmp);
}

void fill_cheb_Dn(int n, double **mat){
	int i, j;

	double *vec = (double *) malloc((n+1) * sizeof(double));
	charray(n, vec);

	mat[0][0] = (2.0 * n * n + 1.0)/6.0;
	for(j = 1; j < n; j++)
		mat[0][j] = 2.0 * ODDM(j)/(1.0 - vec[j]);
	mat[0][n] = 0.5 * ODDM(n);
	for(i = 1; i < n; i++){
		mat[i][0] = -0.5 * ODDM(i)/(1.0 - vec[i]);
		for(j = 1; j < n; j++)
			if(i != j)
				mat[i][j] = ODDM(i+j)/(vec[i] - vec[j]);
		mat[i][n] = 0.5 * ODDM(i+n)/(1.0 + vec[i]);		
		mat[i][i] = - 0.5 * vec[i]/(1.0 - vec[i] * vec[i]);
	}
	mat[n][0] = - 0.5 * ODDM(n);
	for(j = 1; j < n; j++){
		mat[n][j] = -2.0 * ODDM(n+j)/(1.0 + vec[j]);
	}
	mat[n][n] = -(2.0 * n * n + 1.0)/6.0;

	free(vec);
	return;
}

double cheb_Dn1_el(int n, int i, int j){ //Gives out an element of 1st derivative Chebyshev-Gauss-Lobatto matrix
	if(i == 0){
		if(j == 0) return (2.0 * n * n + 1.0)/6.0;
		else if(j == n) return 0.5 * ODDM(n);
		else return 2.0 * ODDM(j)/(1.0 - chpoint(n, j));
	}
	else if (i == n){
		if(j == 0) return - 0.5 * ODDM(n);
		else if(j == n) return -(2.0 * n * n + 1.0)/6.0;
		else return -2.0 * ODDM(n+j)/(1.0 + chpoint(n, j));
	}
	else{
		if(j == 0) return -0.5 * ODDM(i)/(1.0 - chpoint(n, i));
		else if(j == n) return 0.5 * ODDM(i+n)/(1.0 + chpoint(n, i));
		else if(j == i)	return - 0.5 * chpoint(n, j)/(1.0 - pow(chpoint(n, j), 2.0));
		else return ODDM(i+j)/(chpoint(n, i) - chpoint(n, j));
	}
}

double cheb_Dn2_el(int n, int i, int j){ //Gives out an element of 2nd derivative Chebyshev-Gauss-Lobatto matrix
	double xi, xj;
	if(i == 0){
		if(j == 0) return (1.0 * n * n * n * n - 1.0)/15.0;
		else if(j == n) return ODDM(n) * (n * n - 1.0)/3.0;
		else{
			xj = chpoint(n, j);
			return 2.0/3.0 * ODDM(j) * ((2.0 * n * n + 1.0) * (1.0 - xj) - 6.0)/pow(1.0 - xj, 2.0);}
	}
	else if (i == n){
		if(j == 0) return ODDM(n) * (n * n - 1.0)/3.0;
		else if(j == n) return (1.0 * n * n * n * n - 1.0)/15.0;
		else{
			xj = chpoint(n, j);
			return 2.0/3.0 * ODDM(j+n) * ((2.0 * n * n + 1.0) * (1.0 + xj) - 6.0)/pow(1.0 + xj, 2.0);}
	}
	else{
		xi = chpoint(n, i);
		if(j == 0) return 0.5 * ODDM(i) * (xi*xi + xi - 2.0)/(1.0 - xi*xi)/pow(xi - 1.0, 2.0);
		else if(j == n) return 0.5 * ODDM(i+n) * (xi*xi - xi - 2.0)/(1.0 - xi*xi)/pow(xi + 1.0, 2.0);
		else if(j == i)	return -((n*n - 1.0) * (1.0 - xi*xi) + 3.0)/3.0/pow(1.0 - xi*xi, 2.0);
		else {
			xj = chpoint(n, j);
			return ODDM(i+j) * (xi*xi + xi*xj - 2.0)/(1.0 - xi*xi)/pow(xi - xj, 2.0);}			
	}
}

double cheb_Dn_el(int npoint, int ord, int i, int j){//Very inefficient, only works reasonably for ords 1 and 2
	if(ord == 0){
		if(i == j) return 1.0;
		else return 0.0;
	}
	if(ord == 1){
		return cheb_Dn1_el(npoint, i, j);
	}
	if(ord == 2){
		return cheb_Dn2_el(npoint, i, j);
	}	
	if(ord == -1){//No value at all
		return 0.0;
	}
	if(ord < -1){
		printf("Error! In cheb_Dn_el order of differentiation is invalid.\n");
		return 0.0;
	}
	double tmp = 0.0;
	for(int k = 0; k <= npoint; k++)
		tmp += cheb_Dn_el(npoint, ord - 1, i, k) * cheb_Dn1_el(npoint, k, j);
	return tmp;
}

double cheb_Dn_double(bool side, bool parity, int n, int ord, int i, int j){
	if(side){
		if(parity) return cheb_Dn_el(2*n-1, ord, 2*n-i-1, 2*n-j-1) - cheb_Dn_el(2*n-1, ord, 2*n-i-1, j);
		else return cheb_Dn_el(2*n-1, ord, 2*n-i-1, 2*n-j-1) + cheb_Dn_el(2*n-1, ord, 2*n-i-1, j);
	}
	else{
		if(parity) return cheb_Dn_el(2*n-1, ord, n-i-1, n-j-1) - cheb_Dn_el(2*n-1, ord, n-i-1, n+j);
		else return cheb_Dn_el(2*n-1, ord, n-i-1, n-j-1) + cheb_Dn_el(2*n-1, ord, n-i-1, n+j);
	}
}

double cheb_Ex_el(int n, int j, double xst, double xend, double x){
	double xc = (xst + xend - 2*x)/(xend - xst);
	if(xc <= -1.0) xc = M_PI;
	else if(xc >= 1.0) xc = 0.0;
	else xc = acos(xc);
	double denom = xc - M_PI*j/(1.0*n);
	if(denom == 0) return 1.0;

	if(j == 0){
		if(xc == 0.0) return 1.0;
		else return sin(n * xc) / tan(0.5 * xc) / (2.0 * n);
	}

	double sign = -1.0;
	if(j%2) sign = 1.0;

	if(j == n){
		if(xc == M_PI) return 1.0;
		else return sign * sin(n * xc) * tan(0.5 * xc) / (2.0 * n);
	}

	return (sin((n - 0.5)*denom)/sin(0.5*denom) +
				sin((n - 0.5)*(xc + M_PI*j/(1.0*n)))/sin(0.5*(xc + M_PI*j/(1.0*n)))
				- 2.0 * sign * cos(n*xc)) / (2.0 * n);
}

void zero_mat(int nrow, int ncol, double **mat){
	for(int i = 0; i < nrow; i++){
		for(int j = 0; j < ncol; j++)
			mat[i][j] = 0.0;
	}
}

double mult_cheb_Dn(int npoint, int nfunc, int off, int ord, double xst, double xend, int i, double *vec){
	double tmp = 0.0;
	double dydx = 2.0/(xend - xst);

	for(int j = 0; j < npoint + 1; j++)
		tmp += cheb_Dn_el(npoint, ord, npoint - i, npoint - j) * vec[nfunc*j+off];

	return tmp * pow(dydx, ord);
}

void row_cheb_Dn(int npoint, int nfunc, int off, int ord, double xst, double xend, int i, double *vec){
	double dydx = 2.0/(xend - xst);

	for(int j = 0; j < nfunc*(npoint+1); j++)
		vec[j] = 0.0;
	for(int j = 0; j < npoint + 1; j++)
		vec[nfunc*j+off] = cheb_Dn_el(npoint, ord, npoint - i, npoint - j) * pow(dydx, ord);

	return;
}

void row_add_cheb_Dn(int npoint, int nfunc, int off, int ord, double xst, double xend, int i, double coef, double *vec){
	double dydx = 2.0/(xend - xst);

	for(int j = 0; j < npoint + 1; j++)
		vec[nfunc*j+off] += coef * cheb_Dn_el(npoint, ord, npoint - i, npoint - j) * pow(dydx, ord);

	return;
}

//This part relies on fft transformations

void cheb_spect(double *y, int n, double *spect){
	double coef = n;

	fftw_plan planCOS = fftw_plan_r2r_1d(n+1, spect, spect, FFTW_REDFT00, FFTW_ESTIMATE);
	
	for(int i = 0; i <= n; i++)
		spect[i] = y[i]/coef;

	fftw_execute_r2r(planCOS, spect, spect);

	fftw_destroy_plan(planCOS);

	return;
}

void cheb_spect_st(double *y, int n, double *spect){
	double coef = n;

	fftw_plan planCOS = fftw_plan_r2r_1d(n, spect, spect, FFTW_REDFT10, FFTW_ESTIMATE);

	for(int i = 0; i < n; i++)
		spect[i] = y[i]/coef;

	fftw_execute_r2r(planCOS, spect, spect);

	fftw_destroy_plan(planCOS);

	return;
}

void cheb_diff(double xst, double xend, double *specti, int n, double *specto){
	double tmp;
	double sp1 = specti[n-1];
	specto[n-1] = - 4.0*n*specti[n]/(xend-xst);
	specto[n] = 0.0;
	for(int i = n-2; i >= 0; i--){
		tmp = specti[i];
		specto[i] = specto[i+2] - 4.0*(i+1)*sp1/(xend-xst);
		sp1 = tmp;
	}
	return;
}

void cheb_int(double xst, double xend, double *specti, int n, double *specto){
	double tmp1 = specti[0], tmp2;
	for(int i = 1; i < n; i++){
		tmp2 = specti[i];
		specto[i] = -(xend-xst)*(tmp1 - specti[i+1])/(4.0*i);
		tmp1 = tmp2;
	}
	specto[0] = 0.0;
	specto[n] = -(xend-xst)*tmp1/(4.0*n);
	return;
}

double chebev(double xst, double xend, double *spect, int n, double x){
	double d = 0.0, dd = 0.0, sv, y, y2;
	int j;
	if ((x-xst)*(x-xend) > 0.0){
		printf("x is not in the range xst to xend\n");
		return NAN;
	}

	y2 = 2.0 * (y = (xst + xend - 2.0 * x)/(xend - xst)); //Change of variable.
	
	for (j = n-1; j >= 1; j--){ //Clenshaw’s recurrence.
		sv = d;
		d = y2 * d - dd + spect[j];
		dd = sv;
	}

	return y * d - dd + 0.5 * spect[0]; //Last step is different.
}

double chebev_d(double xst, double xend, double *spect, int n, double x){
	double d = 0.0, dd = 0.0, sv, y, y2;
	int j;
	if ((x-xst)*(x-xend) > 0.0){
		printf("x is not in the range xst to xend\n");
		return NAN;
	}

	y2 = 2.0 * (y = (xst + xend - 2.0 * x)/(xend - xst)); //Change of variable.
	
	for (j = n-2; j >= 1; j--){ //Clenshaw’s recurrence.
		sv = d;
		d = (y2 * (j+2) * d - (j+3) * dd)/(1.0*(j+1)) + spect[j+1];
		dd = sv;
	}

	return (4.0 * y * d - 3.0 * dd + spect[1]) * 2.0 / (xst - xend); //Last step is different.
}

double chebev_int(double xst, double xend, double *spect, int n, double x){
	double d = 0.0, dd = 0.0, sv, y, y2;
	int j;
	if ((x-xst)*(x-xend) > 0.0){
		printf("x is not in the range xst to xend\n");
		return NAN;
	}

	y2 = 2.0 * (y = (xst + xend - 2.0 * x)/(xend - xst)); //Change of variable.
	
	for (j = n-1; j >= 1; j--){ //Clenshaw’s recurrence.
		sv = d;
		d = y2 * d - dd - (spect[j-1] - spect[j+1])/(2.0*j);
		dd = sv;
	}

	return 0.5*(y * d - dd)*(xend-xst); //Last step is different.
}

double cheb_integrate(double xst, double xend, double *spect, int n){
	double sum = 0.5*spect[0];
	for(int i = 2; i <= n; i+=2) sum -= spect[i]/(1.0*i*i - 1.0);
	return (xend - xst)*sum;
}

//Fourier spectral methods functions
//Periodic function [0, 2\pi]
double four_Dn1_el(int n, int i, int j){ //Gives out an element of 1st derivative Chebyshev-Gauss-Lobatto matrix
	if(i == j) return 0.0;
	else{
		if((i+j)%2) return -0.5*tan(M_PI*(0.5*n-i+j)/(1.0*n));
		else return 0.5*tan(M_PI*(0.5*n-i+j)/(1.0*n));
	}
}

double four_Dn2_el(int n, int i, int j){ //Gives out an element of 1st derivative Chebyshev-Gauss-Lobatto matrix
	if(i == j) return -(2.0 + n*n)/12.0;
	else{
		if((i+j)%2) return 0.5/pow(sin(M_PI*(i-j)/(1.0*n)), 2.0);
		else return -0.5/pow(sin(M_PI*(i-j)/(1.0*n)), 2.0);
	}
}

double four_Dn_el(int npoint, int ord, int i, int j){
	if(ord == 0){
		if(i == j) return 1.0;
		else return 0.0;
	}
	if(ord == 1){
		return four_Dn1_el(npoint, i, j);
	}
	if(ord == 2){
		return four_Dn2_el(npoint, i, j);
	}	
	if(ord == -1){//No value at all
		return 0.0;
	}
	if(ord < -1){
		printf("Error! In four_Dn_el order of differentiation is invalid.\n");
		return 0.0;
	}
	double tmp = 0.0;
	for(int k = 0; k < npoint; k++)
		tmp += four_Dn_el(npoint, ord - 1, i, k) * four_Dn1_el(npoint, k, j);
	return tmp;
}

//Periodic function [0, \pi]

double cos_Dn_el(int npoint, int ord, int i, int j){
	if(j == 0) return four_Dn_el(2*npoint, ord, i, 0);
	else if(j == npoint) return four_Dn_el(2*npoint, ord, i, npoint);
	else return four_Dn_el(2*npoint, ord, i, j) + four_Dn_el(2*npoint, ord, i, 2*npoint-j);
}

double sin_Dn_el(int npoint, int ord, int i, int j){
	if(j == 0 || j == npoint) return 0.0;
	else return four_Dn_el(2*npoint, ord, i, j) - four_Dn_el(2*npoint, ord, i, 2*npoint-j);
}

double four_DN_Dn_el(int npoint, int ord, int i, int j){
	if(j == 0) return 0.0;
	else if(j == npoint) return four_Dn_el(4*npoint, ord, i, npoint) - four_Dn_el(4*npoint, ord, i, 3*npoint);
	else return four_Dn_el(4*npoint, ord, i, j) + four_Dn_el(4*npoint, ord, i, 2*npoint - j) - four_Dn_el(4*npoint, ord, i, 2*npoint + j) - four_Dn_el(4*npoint, ord, i, 4*npoint - j);
}

double four_ND_Dn_el(int npoint, int ord, int i, int j){
	if(j == 0) return four_Dn_el(4*npoint, ord, i, 0) - four_Dn_el(4*npoint, ord, i, 2*npoint);
	else if(j == npoint) return 0.0;
	else return four_Dn_el(4*npoint, ord, i, j) - four_Dn_el(4*npoint, ord, i, 2*npoint - j) - four_Dn_el(4*npoint, ord, i, 2*npoint + j) + four_Dn_el(4*npoint, ord, i, 4*npoint - j);
}