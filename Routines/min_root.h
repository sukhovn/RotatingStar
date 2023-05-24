#include <iostream>
#include <fstream>

#ifndef FMAX
#define FMAX(a, b) (a > b ? a : b)
#endif

#ifndef FMIN
#define FMIN(a, b) (a > b ? b : a)
#endif

#ifndef SIGN
#define SIGN(a, b) (b >= 0 ? a : -a)
#endif

#ifndef SHFT2
#define SHFT2(a, b, c) (a) = (b); (b) = (c);
#endif

#ifndef SHFT3
#define SHFT3(a, b, c, d) (a) = (b); (b) = (c); (c) = (d);
#endif

#ifndef MIN_ROOT_H
#define MIN_ROOT_H

//Secant solver constants
#define TOLF 10e-10
#define MAXITS 2000 //Limit for maximal number of iterations

//This algorithm implemets secant method.
//It takes xst and xend and function values there as initial guesses and solves equation f(x) = val
//It exits with zero if an equation was solved to a desired accuracy
//It exits with -1 in case of a fail or an error

template<class Callable>
int secant_solve(const Callable &func, double xl, double xr, double *x){
	double yl = func(xl);
	if(fabs(yl) < TOLF){ //Checking if x is a zero
		*x = xl;
		return 0;
	}

	double yr = func(xr);
	if(fabs(yr) < TOLF){ //Checking if x is a zero
		*x = xr;
		return 0;
	}

	double ytmp;
	for(int its = 0; its < MAXITS; its++){
		*x = (xl*yr - xr*yl)/(yr - yl);
		ytmp = func(*x);
		
		if(fabs(ytmp) < TOLF){ //Checking if x is a zero
			return 0;
		}
		if(ytmp*yl > 0.0){
			xl = *x;
			yl = ytmp;
		}
		else{
			xr = *x;
			yr = ytmp;
		}
	}

	std::cout << "No root was found, maximal number of iterations was exceeded\n";
	return -1;
}

#define CGOLD 0.3819660
#define TINY 1.0e-20

//Minimizer constants
#define TOL 1.0e-6

#define GLIMIT 100.0

#define ITMAX 100
#define ZEPS 1.0e-10

template<typename Callable>
class Minimization{
private:
	const Callable &func;

	bool bracket_set;

	double xa;
	double xb;
	double xc;

	double fa;
	double fb;
	double fc;

public:
	double xmin;
	double fmin;
	
	Minimization(const Callable &funci): func(funci), bracket_set(false) {};

	//Setting bracketing manualy
	int set_bracket(double ax, double bx, double cx){
		double dum;
		if(ax == bx || bx == cx || cx == ax){
			std::cout << "Bad bracket choice in Minimization" << std::endl;
			return -1;
		}

		xa = ax; xb = bx; xc = cx;

		if(xa > xb){dum = xa; xa = xb; xb = dum;}//Bubble sorting points
		if(xb > xc){dum = xc; xc = xb; xb = dum;}
		if(xa > xb){dum = xa; xa = xb; xb = dum;}

		fa = func(xa);
		fb = func(xb);
		fc = func(xc);

		if(fb > fa || fb > fc){
			std::cout << "Bad bracket choice in Minimization" << std::endl;
			std::cout << "(" << xa << ", " << fa << "), ";
			std::cout << "(" << xb << ", " << fb << "), ";
			std::cout << "(" << xc << ", " << fc << ")" << std::endl;
			return -1;
		}

		bracket_set = true;
		return 0;
	}
	
	//Searches the extremum via golden ratio search
	double golden(void){
		if(bracket_set == false){
			std::cout << "Bracketing is not set, minimization fail in Minimization" << std::endl;
			return 0.0;
		}

		double x0 = xa;
		double x3 = xc;
		double f1, f2, x1, x2;
		if(fabs(x3 - xb) > fabs(xb - x0)){
			x1 = xb;
			x2 = xb + CGOLD * (x3 - xb);
		}
		else{
			x2 = xb;
			x1 = xb - CGOLD * (xb - x0);		
		}
		f1 = func(x1);
		f2 = func(x2);
		while (fabs(x3 - x0) > TOL * (fabs(x1) + fabs(x2))){
			if (f2 < f1){
				SHFT3(x0, x1, x2, (1.0 - CGOLD) * x1 + CGOLD * x3);
				SHFT2(f1, f2, func(x2));
			} else{
				SHFT3(x3, x2, x1, (1.0 - CGOLD) * x2 + CGOLD * x0);
				SHFT2(f2, f1, func(x1));
			}
		}
		if(f1 < f2){
			xmin = x1;
			return f1;
		} else {
			xmin = x2;
			return f2;
		}
	}

	//Searches the extremum via Brent method
	double brent(void){
		if(bracket_set == false){
			std::cout << "Bracketing is not set, minimization fail in Minimization" << std::endl;
			return 0.0;
		}
		int iter;
		double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
		double e = 0.0;

		a = (xa < xc ? xa : xc);
		b = (xa > xc ? xa : xc);
		x = w = v = xb;
		fw = fv = fx = func(x);

		for(iter = 1; iter <= ITMAX; iter++){
			xm = 0.5 * (a + b);
			tol2 = 2.0 * (tol1 = TOL * fabs(x) + ZEPS);
			if(fabs(x - xm) <= (tol2 - 0.5 * (b-a))){
				xmin = x;
				fmin = fx;
				return xmin;
			}

			if(fabs(e) > tol1){
				r = (x - w) * (fx - fv);
				q = (x - v) * (fx - fw);
				p = (x - v) * q - (x - w) * r;
				q = 2.0 * (q - r);
				if (q > 0.0) p = -p;
				q = fabs(q);
				etemp = e;
				e = d;

				if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q*(a-x) || p >= q*(b-x)) d = CGOLD * (e = (x >= xm ? a-x : b-x));
				else {
					d = p/q;
					u = x+d;
					if (u-a < tol2 || b-u < tol2)
						d = SIGN(tol1,xm-x);
				}
			} else {
				d = CGOLD * (e = (x >= xm ? a-x : b-x));
			}
			u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
			fu = func(u);
			if (fu <= fx) {
				if (u >= x) a=x; else b=x;
				SHFT3(v, w, x, u);
				SHFT3(fv, fw, fx, fu);
			} else {
				if (u < x) a = u; else b=u;
				if (fu <= fw || w == x) {
					v=w;
					w=u;
					fv=fw;
					fw=fu;
				} else if (fu <= fv || v == x || v == w) {
					v=u;
					fv=fu;
				}
			}
		}
		std::cout << "Too many iterations in brent\n";
		xmin = x;
		fmin = fx;
		return xmin;
	}
};

#endif