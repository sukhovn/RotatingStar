#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <interpolation.h>

//This routine closely follows Numerical Recipies, it is a realization of Neville interpolation algorithm
//The algorithm takes two arrays with n points xa and ya and uses them to find
//an interpolated value y at the point x, giving an error estimate dy
int polint(double *xa, double *ya, int n, double x, double *y, double *dy){
	int i, m, ns = 0;
	double den, dif, dift, ho, hp, w;

	double *c = (double *) malloc(n * sizeof(double));
	double *d = (double *) malloc(n * sizeof(double));

	dif = fabs(x - xa[0]);
	for(i = 0; i < n; i++){
		if((dift = fabs(x - xa[i])) < dif){//Finding the cl
			ns = i;
			dif = dift;
		}
		c[i] = d[i] = ya[i]; 
	}

	*y = ya[ns--];
	for(m = 1; m < n; m++){
		for(i = 0; i < n-m; i++){
			ho = xa[i] - x;
			hp = xa[i+m] - x;
			w = c[i+1] - d[i];
			if((den = ho - hp) == 0.0) {
				printf("Error in routine polint\n");
				return -1;
			}
			den = w/den;
			d[i] = hp * den;
			c[i] = ho * den;
		}
		*y += (*dy = (2 * ns < (n-m-2) ? c[ns+1] : d[ns--]));
	}

	free(c);
	free(d);
	return 0;
}

//This routine closely follows Numerical Recipies, it is a realization of cubic spline interpolation algorithm
//It takes a set of n function values y on a set of points x together with boundary first derivatives yp1 and ypn
//If either yp1 or ypn is greateer then 0.99e30 it discards it and assumes the second derivative at the boundary to be zero
//The routine then finds quibic spline second derivatives at the points x 
void spline(double *x, double *y, int n, double yp1, double ypn, double *y2){
	int i;
	double p, qn, sig, un;

	double *u = (double *) malloc((n-1) * sizeof(double));
	//Setting up y2[0] or the equation
	if (yp1 > 0.99e30)
		y2[0] = u[0] = 0;
	else{
		y2[0] = -0.5;
		u[0] = 3.0/(x[1] - x[0]) * ((y[1] - y[0])/(x[1] - x[0]) - yp1);
	}

	for(i = 1; i < n-1; i++){
		sig = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
		p = sig * y2[i-1] + 2.0;
		y2[i] = (sig - 1.0)/p;
		u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1]);
		u[i] = (6.0 * u[i]/(x[i+1] - x[i-1]) - sig * u[i-1])/p;
	}

	if(ypn > 0.99e30)
		qn = un = 0.0;
	else{
		qn = 0.5;
		un = 3.0/(x[n-1] - x[n-2]) * (ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
	}
	y2[n-1] = (un - qn * u[n-2])/(qn * y2[n-2] + 1.0);
	for(i = n-2; i >= 0; i--)
		y2[i] = y2[i] * y2[i+1] + u[i];

	free(u);
	return;
}

//The same function but finds spline second derivatives for multiple functions
void spline_multiple(double *x, double *y, int n, int nfunc, double *yp1, double *ypn, double *y2){
	int i;
	double p, qn, sig, un;

	double *u = (double *) malloc(nfunc*(n-1) * sizeof(double));

	//Setting up y2[0] or the equation
	for(int fn = 0; fn < nfunc; fn++){
		if (yp1[fn] > 0.99e30)
			y2[fn] = u[fn] = 0.0;
		else{
			y2[fn] = -0.5;
			u[fn] = 3.0/(x[1] - x[0]) * ((y[nfunc+fn] - y[fn])/(x[1] - x[0]) - yp1[fn]);
		}
	}

	for(i = 1; i < n-1; i++){
		sig = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
		for(int fn = 0; fn < nfunc; fn++){
			int indx = nfunc*i+fn;
			p = sig * y2[indx-nfunc] + 2.0;
			y2[indx] = (sig - 1.0)/p;
			u[indx] = (y[indx+nfunc] - y[indx])/(x[i+1] - x[i]) - (y[indx] - y[indx-nfunc])/(x[i] - x[i-1]);
			u[indx] = (6.0 * u[indx]/(x[i+1] - x[i-1]) - sig * u[indx-nfunc])/p;
		}
	}

	for(int fn = 0; fn < nfunc; fn++){
		if(ypn[fn] > 0.99e30)
			qn = un = 0.0;
		else{
			qn = 0.5;
			un = 3.0/(x[n-1] - x[n-2]) * (ypn[fn] - (y[nfunc*(n-1)+fn] - y[nfunc*(n-2)+fn])/(x[n-1] - x[n-2]));
		}
		y2[nfunc*(n-1)+fn] = (un - qn * u[nfunc*(n-2)+fn])/(qn * y2[nfunc*(n-2)+fn] + 1.0);
	}

	for(i = n-2; i >= 0; i--){
		for(int fn = 0; fn < nfunc; fn++){
			int indx = nfunc*i+fn;
			y2[indx] = y2[indx] * y2[indx+nfunc] + u[indx];
		}
	}

	free(u);
	return;
}


//This routine closely follows Numerical Recipies, it performs cubic spline
//It uses ya function values at xa point together with the derivatives y2a from previous routine
//It interpolates value y at the point x

int splint(double *xa, double *ya, double *y2a, int n, double x, double *y){
	int klo, khi, k;
	double h, a, b;

	klo = 0; khi = n-1;
	while(khi - klo > 1){
		k = (khi + klo) >> 1;
		if (xa[k] > x) khi = k;
		else klo = k;
	}
	h = xa[khi] - xa[klo];
	if(h == 0.0){
		printf("Bad xa input to routine\n");
		return -1;
	}
	a = (xa[khi] - x)/h;
	b = (x - xa[klo])/h;
	*y = a * ya[klo] + b *ya[khi] + ((a*a*a - a) * y2a[klo] + (b*b*b - b) * y2a[khi]) * h*h/6.0;
	return 0;
}

//The same function but interpolates multiple functions at the same time

int splint_multiple(double *xa, double *ya, double *y2a, int n, int nfunc, double x, double *y){
	int klo, khi, k;
	double h, a, b;

	klo = 0; khi = n-1;
	while(khi - klo > 1){
		k = (khi + klo) >> 1;
		if (xa[k] > x) khi = k;
		else klo = k;
	}
	h = xa[khi] - xa[klo];
	if(h == 0.0){
		printf("Bad xa input to routine\n");
		return -1;
	}
	a = (xa[khi] - x)/h;
	b = (x - xa[klo])/h;
	for(int fn = 0; fn < nfunc; fn++)
		y[fn] = a * ya[nfunc*klo+fn] + b *ya[nfunc*khi+fn] + ((a*a*a - a) * y2a[nfunc*klo+fn] + (b*b*b - b) * y2a[nfunc*khi+fn]) * h*h/6.0;
	return 0;
}

//This routine performs cubic spline and returns the derivative
//It uses ya function values at xa point together with the derivatives y2a from previous routine
//It interpolates value of the derivative dy/dx at the point x

int dsplint(double *xa, double *ya, double *y2a, int n, double x, double *yp){
	int klo, khi, k;
	double h, a, b;

	klo = 0; khi = n-1;
	while(khi - klo > 1){
		k = (khi + klo) >> 1;
		if (xa[k] > x) khi = k;
		else klo = k;
	}
	h = xa[khi] - xa[klo];
	if(h == 0.0){
		printf("Bad xa input to routine\n");
		return -1;
	}
	a = (xa[khi] - x)/h;
	b = (x - xa[klo])/h;
	*yp = (ya[khi] - ya[klo])/h + ((3.0*b*b - 1.0) * y2a[khi] - (3.0*a*a - 1.0) * y2a[klo]) * h/6.0;
	return 0;
}

//The same function but interpolates values of first derivatives of multiple functions at the same time

int dsplint_multiple(double *xa, double *ya, double *y2a, int n, int nfunc, double x, double *yp){
	int klo, khi, k;
	double h, a, b;

	klo = 0; khi = n-1;
	while(khi - klo > 1){
		k = (khi + klo) >> 1;
		if (xa[k] > x) khi = k;
		else klo = k;
	}
	h = xa[khi] - xa[klo];
	if(h == 0.0){
		printf("Bad xa input to routine\n");
		return -1;
	}
	a = (xa[khi] - x)/h;
	b = (x - xa[klo])/h;
	for(int fn = 0; fn < nfunc; fn++)
		yp[fn] = (ya[nfunc*khi+fn] - ya[nfunc*klo+fn])/h + ((3.0*b*b - 1.0) * y2a[nfunc*khi+fn] - (3.0*a*a - 1.0) * y2a[nfunc*klo+fn]) * h/6.0;
	return 0;
}