#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

//This routine closely follows Numerical Recipies, it is a realization of Neville interpolation algorithm
//The algorithm takes two arrays with n points xa and ya and uses them to find
//an interpolated value y at the point x, giving an error estimate dy

int polint(double *xa, double *ya, int n, double x, double *y, double *dy);

//This routine closely follows Numerical Recipies, it is a realization of cubic spline interpolation algorithm
//It takes a set of n function values y on a set of points x together with boundary first derivatives yp1 and ypn
//If either yp1 or ypn is greateer then 0.99e30 it discards it and assumes the second derivative at the boundary to be zero
//The routine then finds quibic spline second derivatives at the points x 

void spline(double *x, double *y, int n, double yp1, double ypn, double *y2);

//The same function but finds spline second derivatives for multiple functions

void spline_multiple(double *x, double *y, int n, int nfunc, double *yp1, double *ypn, double *y2);

//This routine closely follows Numerical Recipies, it performs cubic spline
//It uses ya function values at xa point together with the derivatives y2a from previous routine
//It interpolates value y at the point x

int splint(double *xa, double *ya, double *y2a, int n, double x, double *y);

//The same function but interpolates multiple functions at the same time

int splint_multiple(double *xa, double *ya, double *y2a, int n, int nfunc, double x, double *y);

//This routine performs cubic spline and returns the derivative
//It uses ya function values at xa point together with the derivatives y2a from previous routine
//It interpolates value of the derivative dy/dx at the point x

int dsplint(double *xa, double *ya, double *y2a, int n, double x, double *yp);

//The same function but interpolates values of first derivatives of multiple functions at the same time

int dsplint_multiple(double *xa, double *ya, double *y2a, int n, int nfunc, double x, double *yp);

#endif