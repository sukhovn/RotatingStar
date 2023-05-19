#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#ifndef SPECTRAL_H
#define SPECTRAL_H

//Chebyshev spectral methods functions

double chpoint(int n, int i);
//Gives out Chebyshev-Gauss-Lobatto point

void charray(int n, double *vec);
//Fills out an array with Chebyshev-Gauss-Lobatto points

double ch_x(double xst, double xend, int n, int i);
//Gives out Chebyshev-Gauss-Lobatto point adjusted to [xst, xend] segment

void ch_x_array(double xst, double xend, int n, double *vec);
//Fills out an array with Chebyshev-Gauss-Lobatto points adjusted to [xst, xend] segment

int ch_x_inv_low(double xst, double xend, int n, double x);
//Inverts ch_x expression and finds the closest to x CGL point from the left

int ch_x_inv_high(double xst, double xend, int n, double x);
//Inverts ch_x expression and finds the closest to x CGL point from the right

//The same set of routines on a staggered grid
//Collocation points are now zeroes of Chebyshev T_n(x) polynomials

double chpoint_st(int n, int i);
//Gives out Chebyshev-Gauss-Lobatto point on a staggered grid

void charray_st(int n, double *vec);
//Fills out an array with Chebyshev-Gauss-Lobatto points on a staggered grid

double ch_x_st(double xst, double xend, int n, int i);
//Gives out Chebyshev-Gauss-Lobatto point adjusted to [xst, xend] segment on a staggered grid

void ch_x_array_st(double xst, double xend, int n, double *vec);
//Fills out an array with Chebyshev-Gauss-Lobatto points adjusted to [xst, xend] segment on a staggered grid

//Matrix routines

void fill_cheb_Dn(int n, double **mat);
//Makes nth order Chebyshev differentiation matrix

double cheb_Dn1_el(int n, int i, int j);
//Gives out an element of 1st derivative Chebyshev differentiation matrix

double cheb_Dn2_el(int n, int i, int j);
//Gives out an element of 2nd derivative Chebyshev differentiation matrix

double cheb_Dn_el(int npoint, int ord, int i, int j);
//Gives out an element of ord-th derivative Chebyshev differentiation matrix, explicit for ord = 1 and ord = 2, recursively for higer orders

double cheb_Dn_double(bool side, bool parity, int n, int ord, int i, int j);
//Gives out an element of ord-th derivative Chebyshev differentiation matrix for a doubled domain

double cheb_Ex_el(int n, int j, double xst, double xend, double x);
//Gives extrapolation coeffient to point x of the collocation polynomial 

void zero_mat(int nrow, int ncol, double **mat);
//Zeroes nrow x ncol matrix

double mult_cheb_Dn(int npoint, int nfunc, int off, int ord, double xst, double xend, int i, double *vec);
//Multiplies a vector *vec by i-th row of Chebyshev D-ord matrix, following the mask 
//vec[nfunc*i+off]

void row_cheb_Dn(int npoint, int nfunc, int off, int ord, double xst, double xend, int i, double *vec);
//Fills a row *vec with i-th row of Chebyshev D-ord matrix, following the mask 
//vec[nfunc*i+off]

void row_add_cheb_Dn(int npoint, int nfunc, int off, int ord, double xst, double xend, int i, double coef, double *vec);
//Adds to a *vec an i-th row of Chebyshev D-ord matrix multiplied by coef, following the mask 
//vec[nfunc*i+off]

//This part relies on fft transformations

void cheb_spect(double *y, int n, double *spect);
//Finds Chebyshev polynomial coefficients *spect given a function set on Chebyshev-Gauss-Lobatto grid

void cheb_spect_st(double *y, int n, double *spect);
//Finds Chebyshev polynomial coefficients *spect given a function set on a staggered Chebyshev-Gauss-Lobatto grid

void cheb_diff(double xst, double xend, double *specti, int n, double *specto);
//Finds Chebyshev polynomial coefficients for the derivative *spect given coefficients for the function

void cheb_int(double xst, double xend, double *specti, int n, double *specto);
//Finds Chebyshev polynomial coefficients for the antiderivative *spect given coefficients for the function

double chebev(double xst, double xend, double *spect, int n, double x);
//Evaluates a function at a arbitrary point x between xst and xend given Chebyshev polynomial expansion coefficients

double chebev_d(double xst, double xend, double *spect, int n, double x);
//Evaluates function's 1st derivative at an arbitrary point x between xst and xend given Chebyshev polynomial expansion coefficients

double chebev_int(double xst, double xend, double *spect, int n, double x);
//Evaluates function's integral at an arbitrary point x between xst and xend given Chebyshev polynomial expansion coefficients

double cheb_integrate(double xst, double xend, double *spect, int n);
//Evaluates function's integral at an arbitrary point x between xst and xend given Chebyshev polynomial expansion coefficients

//Fourier spectral methods functions

double four_Dn1_el(int n, int i, int j);
//Gives out an element of 1st derivative Fourier differentiation matrix

double four_Dn2_el(int n, int i, int j);
//Gives out an element of 2nd derivative Fourier differentiation matrix

double four_Dn_el(int npoint, int ord, int i, int j);
//Gives out an element of ord-th derivative Fourier differentiation matrix, explicit for ord = 1 and ord = 2, recursively for higer orders

double cos_Dn_el(int npoint, int ord, int i, int j);
//Gives out an element of ord-th derivative cosine differentiation matrix, explicit for ord = 1 and ord = 2, recursively for higer orders

double sin_Dn_el(int npoint, int ord, int i, int j);
//Gives out an element of ord-th derivative sine differentiation matrix, explicit for ord = 1 and ord = 2, recursively for higer orders

double four_DN_Dn_el(int npoint, int ord, int i, int j);
//Gives out an element of ord-th derivative sine differentiation matrix, explicit for ord = 1 and ord = 2, recursively for higer orders

double four_ND_Dn_el(int npoint, int ord, int i, int j);
//Gives out an element of ord-th derivative sine differentiation matrix, explicit for ord = 1 and ord = 2, recursively for higer orders

#endif