#include <cmath>
#include <limits>
#include <container.h>

#ifndef EQUATIONS_H
#define EQUATIONS_H

class Radial_Compactification{
private:
	double xst;
	double cc;

	std::vector<bool> domains_compactified;
public:
	void adjust_compactification(const Domain_Partitioning &part, bool print_selection = false){
		int ndomains = 0;
		xst = part.partition_start(0);
		double xend = part.partition_end(0);
		domains_compactified.clear();
		for(const Bulk_Domain *dom = part.first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			if(dom->domain_start(0) > xst) xst = dom->domain_start(0);
			++ndomains;
			domains_compactified.push_back(false);
		}
		if(print_selection){
			std::cout << "Selected compactification region: [" << xst << ',' << xend << ']' << std::endl;
		}
		cc = 1.0/(xend - xst);
		for(const Bulk_Domain *dom = part.first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			if(dom->domain_end(0) > xst + COMP_TOL)
				domains_compactified[dom->bulk_domain_number] = true;
		}
	}

	Radial_Compactification(const Domain_Partitioning &part, bool print_selection = false){
		adjust_compactification(part, print_selection);
	}

	double xc(int dm, double *x) const{
		if(dm == 0){
			double xshift = x[0] - xst;
			if(xshift <= 0.0) return x[0];
			if(xshift*cc >= 1.0) return std::numeric_limits<double>::infinity();
			return xst + xshift/(1 - xshift*cc);
		}
		else
			return x[dm];
	}

	double dxc_factor(double x) const{
		double xshift = x - xst;
		if(xshift <= 0.0) return 1.0;
		if(xshift*cc >= 1.0) return std::numeric_limits<double>::infinity();
		return 1.0/std::pow(1 - xshift*cc, 2.0);
	}

	double dxc(int dnum, int dm, double *x) const{
		if(!domains_compactified[dnum])
			return 1.0;
		if(dm == 0){
			double xshift = x[0] - xst;
			if(xshift*cc >= 1.0) return std::numeric_limits<double>::infinity();
			return 1.0/std::pow(1 - xshift*cc, 2.0);
		}
		else
			return 1.0;
	}

	double ddxc(int dnum, int dm, double *x) const{
		if(!domains_compactified[dnum])
			return 0.0;
		if(dm == 0){
			double xshift = x[0] - xst;
			if(xshift*cc == 1.0) return std::numeric_limits<double>::infinity();
			return 2.0*cc/std::pow(1 - xshift*cc, 3.0);
		}
		else
			return 0.0;
	}

	double xc_inv(int dm, double *x) const{
		if(dm == 0){
			double xshift = x[0] - xst;
			if(xshift <= 0.0) return x[dm];
			return xst + xshift/(1 + xshift*cc);
		}
		else
			return x[dm];
	}
};


class Rotating_Boson_Star_Equations{
public:
	std::vector<double> poly_coef;

	double kappa;
	double omega;
	double m_ang;

	Rotating_Boson_Star_Equations(): poly_coef({1.0}), kappa(1.0), omega(1.0), m_ang(1.0) {};

	double vpot(const double phi) const{
		int deg = poly_coef.size()-1;
		double tmp = poly_coef[deg];
		for(int i = deg; i > 0; i--){
			tmp *= phi;
			tmp += poly_coef[i-1];
		}
		return phi*tmp;
	}

	double dvpot(const double phi) const{
		int deg = poly_coef.size()-1;
		double tmp = (deg+1)*poly_coef[deg];
		for(int i = deg; i > 0; i--){
			tmp *= phi;
			tmp += i*poly_coef[i-1];
		}
		return tmp;
	}

	//Here we specify equations, the functions go in the order
	// \gamma	phi[0]
	// \rho		phi[1]
	// \sigma	phi[2]
	// \omega	phi[3]
	// \phi		phi[4]

	double equations(int fn, double *x, double *phi, double **dphi, double **ddphi) const{
		double sin_th = std::sin(x[1]);
		double cos_th = std::cos(x[1]);
		double ctg_th = cos_th/sin_th;
		double source;

		switch(fn){
			case 0:
				source = exp(-phi[0]) * (exp(-phi[1]) * std::pow(omega + m_ang*phi[3], 2.0) - std::exp(phi[1])*m_ang*m_ang/(x[0]*x[0]*sin_th*sin_th))*phi[4]*phi[4] - vpot(phi[4]*phi[4]);
				return ddphi[0][0] + (3.0/x[0] + dphi[0][0])*dphi[0][0] + 
					(ddphi[1][0] + (2.0*ctg_th + dphi[1][0])*dphi[1][0])/(x[0]*x[0]) - 2.0*kappa*exp(2.0*phi[2])*source;
				break;
			case 1:
				source = exp(-phi[0]) * (exp(-phi[1]) * pow(omega + m_ang*phi[3], 2.0) + exp(phi[1])*m_ang*m_ang/(x[0]*x[0]*sin_th*sin_th))*phi[4]*phi[4];
				return ddphi[0][1] + (2.0/x[0] + dphi[0][0])*dphi[0][1] + 
					(ddphi[1][1] + (ctg_th + dphi[1][0])*dphi[1][1])/(x[0]*x[0]) - (dphi[0][0] + ctg_th*dphi[1][0]/x[0])/x[0] -
					exp(-2.0*phi[1])*(pow(x[0]*dphi[0][3], 2.0) + pow(dphi[1][3], 2.0))*sin_th*sin_th - 2.0*kappa*exp(2.0*phi[2])*source;
				break;
			case 2:
				source = exp(-phi[0]) * (exp(-phi[1]) * pow(omega + m_ang*phi[3], 2.0) - exp(phi[1])*m_ang*m_ang/(x[0]*x[0]*sin_th*sin_th))*phi[4]*phi[4] + exp(-2.0*phi[2])*(pow(dphi[0][4], 2.0) + pow(dphi[1][4]/x[0], 2.0));
				return ddphi[0][2] + dphi[0][2]/x[0] + ddphi[1][2]/(x[0]*x[0]) - 0.25*(2.0/x[0] + dphi[0][0])*dphi[0][0] - 0.25*(2.0*ctg_th + dphi[1][0])*dphi[1][0]/(x[0]*x[0]) - 
					0.25*(2.0/x[0] - dphi[0][1])*dphi[0][1] - 0.25*(2.0*ctg_th - dphi[1][1])*dphi[1][1]/(x[0]*x[0]) -
					0.25*exp(-2.0*phi[1])*(pow(dphi[0][3]*x[0], 2.0) + pow(dphi[1][3], 2.0))*sin_th*sin_th + kappa*exp(2.0*phi[2])*source;
				break;
			case 3:
				source = 2.0/(x[0]*x[0]*sin_th*sin_th) * exp(phi[1] - phi[0])*m_ang*(omega + m_ang*phi[3])*phi[4]*phi[4];
				return ddphi[0][3] + (4.0/x[0] + dphi[0][0] - 2.0 * dphi[0][1])*dphi[0][3] +
					(ddphi[1][3] + (3.0*ctg_th + dphi[1][0] - 2.0*dphi[1][1])*dphi[1][3])/(x[0]*x[0]) - 2.0*kappa*exp(2.0*phi[2])*source;
				break;
			case 4:
				return ddphi[0][4] + (2.0/x[0] + dphi[0][0])*dphi[0][4] + (ddphi[1][4] + (ctg_th + dphi[1][0])*dphi[1][4])/(x[0]*x[0]) + 
					exp(2.0*phi[2])*(exp(-phi[0]) * (exp(-phi[1]) * pow(omega + m_ang*phi[3], 2.0) - exp(phi[1])*m_ang*m_ang/(x[0]*x[0]*sin_th*sin_th)) - dvpot(phi[4]*phi[4]))*phi[4];
				break;
			default:
				std::cout << "Wrong equation number" << std::endl;
				throw("Wrong equation number");
		}
	}

	double boundary_condition_left(int dm, int fn, double *phi, double *dphi) const{
		switch(fn){
			case 0:
				if(dm == 0) return dphi[0];
				else if(dm == 1) return dphi[0];
				break;
			case 1:
				if(dm == 0) return dphi[1];
				else if(dm == 1) return dphi[1];
				break;
			case 2:
				if(dm == 0) return 2.0 * phi[2] + phi[1] - phi[0];
				else if(dm == 1) return 2.0 * phi[2] + phi[1] - phi[0];
				break;
			case 3:
				if(dm == 0) return dphi[3];
				else if(dm == 1) return dphi[3];
				break;
			case 4:
				if(dm == 0) return phi[4];
				else if(dm == 1) return phi[4];
				break;
		}
		std::cout << "Wrong equation number" << std::endl;
		throw("Wrong equation number");
	}

	double boundary_condition_right(int dm, int fn, double *phi, double *dphi) const{
		switch(fn){
			case 0:
				if(dm == 0) return phi[0];
				else if(dm == 1) return dphi[0];
				break;
			case 1:
				if(dm == 0) return phi[1];
				else if(dm == 1) return dphi[1];
				break;
			case 2:
				if(dm == 0) return phi[2];
				else if(dm == 1) return dphi[2];
				break;
			case 3:
				if(dm == 0) return phi[3];
				else if(dm == 1) return dphi[3];
				break;
			case 4:
				if(dm == 0) return phi[4];
				else if(dm == 1) return dphi[4];
				break;
		}
		std::cout << "Wrong equation number" << std::endl;
		throw("Wrong equation number");
	}

	double residuals(int fn, double *x, double *phi, double **dphi, double **ddphi) const{
		double sin_th = std::sin(x[1]);
		double cos_th = std::cos(x[1]);
		double ctg_th = cos_th/sin_th;
		double source;

		switch(fn){
			case 0:
				source = 2.0*dphi[0][4]*dphi[1][4];
				return ddphi[2][0] + 0.5*dphi[0][0]*dphi[1][0] + (0.5*ctg_th - dphi[1][2])*dphi[0][0] - (0.5/x[0] + dphi[0][2])*dphi[1][0] -
					0.5*ctg_th*dphi[0][1] + 0.5*dphi[0][1]*dphi[1][1] - 0.5*dphi[1][1]/x[0] - ctg_th*dphi[0][2] - dphi[1][2]/x[0] -
					0.5*exp(-2.0*phi[1])*dphi[0][3]*dphi[1][3]*x[0]*x[0]*sin_th*sin_th + kappa*source;
				break;
			case 1:
				source = 2.0*exp(-2.0*phi[2])*pow(dphi[0][4], 2.0) + vpot(phi[4]*phi[4]);
				return ddphi[0][0] + 0.5*(2.0/x[0] + dphi[0][0])*dphi[0][0] + ddphi[0][2] - dphi[0][0]*dphi[0][2] + 
					(ddphi[1][2] + (ctg_th + dphi[1][0])*dphi[1][2])/(x[0]*x[0]) - 0.5*(2.0/x[0] - dphi[0][1])*dphi[0][1] - 
					0.5*exp(-2.0*phi[1])*pow(dphi[0][3]*x[0]*sin_th, 2.0) + kappa*exp(2.0*phi[2])*source;
				break;
			default:
				std::cout << "Wrong residual number" << std::endl;
				throw("Wrong residual number");
				break;
		}
	}

	double densities(int fn, double *x, double *phi, double **dphi, double **ddphi) const{
		double sin_th = sin(x[1]);
		switch(fn){
			case 0://This is Komar density
				return exp(2.0*phi[2])*(2.0*exp(-phi[1])*omega*(omega + m_ang*phi[3])*phi[4]*phi[4] - exp(phi[0])*vpot(phi[4]*phi[4]));
				break;
			case 1://This is charge density j_\mu n^\mu
				return exp(2.0*phi[2] - phi[1])*(omega + m_ang*phi[3])*phi[4]*phi[4];
				break;
			case 2://This is ADM \rho = T_{\mu\nu} n^\mu n^\nu energy density
				return exp(-phi[0])*(exp(-phi[1])*pow(omega + m_ang*phi[3], 2.0) + exp(phi[1])*m_ang*m_ang/pow(x[0]*sin_th, 2.0))*phi[4]*phi[4] + exp(-2.0*phi[2])*(pow(dphi[0][4], 2.0) + pow(dphi[1][4]/x[0], 2.0)) + vpot(phi[4]*phi[4]);
				break;
			default:
				std::cout << "Wrong density number" << std::endl;
				throw("Wrong density number");
				break;
		}
	}

	double area_element(double *x) const{
		return x[0]*x[0]*std::sin(x[1]);
	}

	double jacoby_p(unsigned int n, double a, double x) const{
		if(n == 0) return 1.0;
		if(n == 1) return 1.0 + 0.5*a + 0.5*(a + 2.0)*(x - 1.0);
		double p0 = 1.0, p1 = 1.0 + 0.5*a + 0.5*(a + 2.0)*(x - 1.0), tmp;
		for(int i = 2; i <= n; ++i){
			tmp = (2*i+a)*((2*i+a-1)*x*p1 - 0.5*(2*i+a-2)*p0)/(2.0*i*(i+a));
			p0 = p1;
			p1 = tmp;
		}
		return p1;
	}

	double rho_mode(unsigned n, double *x, double *phi){
		double sin_th = sin(x[1]), cos_th = cos(x[1]);
		return 2.0*(2.0*n+1.0)/2.0*phi[1]*jacoby_p(n, 0, cos_th)*sin_th;
	}

	double b_mode(unsigned n, double *x, double *phi){
		double sin_th = sin(x[1]), cos_th = cos(x[1]);
		double result = 2.0*(n+1.0)*(exp(phi[0])-1.0)*jacoby_p(n, 1, cos_th)*sin_th*sin_th/M_PI;
		for(int i = 0; i < n; ++i)//Multiplying by the normalization factor
			result *= (i + 1.0)*(i + 2.0)/((i + 1.5)*(i + 1.5));
		return 2.0*result;
	}

	double omega_mode(unsigned n, double *x, double *phi){
		double sin_th = sin(x[1]), cos_th = cos(x[1]);
		return 2.0*(2.0*n+3.0)*(n+2.0)/(8.0*(n+1.0))*phi[3]*jacoby_p(n, 2, cos_th)*sin_th*sin_th*sin_th;
	}
};

#endif