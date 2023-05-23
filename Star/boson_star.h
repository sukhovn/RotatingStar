#include <container.h>
#include <parameter_file.h>
#include <collocation_routines.h>
#include <collocation_solver_full.h>
#include <equations.h>
#include <unordered_map>

#ifndef BOSON_STAR_H
#define BOSON_STAR_H

#define RAD_FRACTION 0.99

class Boson_Star: public Container, public Rotating_Boson_Star_Equations{
private:
	Parameter_File pfile;

public:
	std::unordered_map<std::string, double> params;

	std::vector<double> residual_error;
	double coordinate_komar_radius;
	double ps_length;

	//Standard constructor
	Boson_Star(const Domain_Partitioning &part);
	
	void read_parameters(const Parameter_File &reader);
	void record_parameters(Parameter_File &reader);

	void load(const std::string &folder_name){
		pfile.read_file(folder_name + "/parameters.dat");
		read_parameters(pfile);
		Container::load(folder_name);
		return;
	}

	//Load from file constructor
	Boson_Star(const std::string &folder_name);

	//Copy constructor
	Boson_Star(const Boson_Star &star);

	Boson_Star &operator=(const Boson_Star &star);

	//Saves star solution
	void save(const std::string &folder_name){
		std::cout << "Saving the solution" << std::endl;
		Container::save(folder_name);
		record_parameters(pfile);
		pfile.print_file(folder_name + "/parameters.dat");
	}

	void regrid(const Domain_Partitioning &part_new){
		Radial_Compactification comp_old(partition);
		Radial_Compactification comp_new(part_new);
		Container::regrid(part_new, comp_old, comp_new);
	}

	friend class Boson_Star_Solver;
	friend class Solution_Saver;
	friend class Extrapolation;
};

class Boson_Star_Solver: private Radial_Compactification, public Collocation_Solver_Full<Radial_Compactification, Rotating_Boson_Star_Equations>{
private:
	int density_extrapolation_resolution;
	Boson_Star *star;

public:
	Boson_Star_Solver(Boson_Star *_star):
		Radial_Compactification(_star->partition),
		Collocation_Solver_Full<Radial_Compactification, Rotating_Boson_Star_Equations>(*this, _star, _star),
		star(_star),
		density_extrapolation_resolution(50) {};

	void change_star(Boson_Star *stari){
		star = stari;
		eq_system = stari;
		change_container(*stari);
	}

	void find_star_parameters(void){
		star->params["omega"] = star->omega;
		std::vector<double> approx_min;

		auto factor_function = [](int fn, double *x, double *phi, double **dphi, double **ddphi)
		{
			return std::exp(0.5*(phi[0] - phi[1]));
		};

		Evaluator metric_prefactor(evaluate_function(factor_function, 1, 1.0));

		//Calculating density functions
		auto density_function = [this](int fn, double *x, double *phi, double **dphi, double **ddphi)
		{
			if(x[0] == 0.0 || x[1] == 0) return 0.0;
			double xc = xc_inv(0, x);
			double factor = 1.0;
			if(fn == 1) factor = star->m_ang;
			if(fn == 2) factor = std::exp(2.0 * phi[2] + 0.5 * (phi[0] - phi[1]));
			return star->kappa*factor*star->densities(fn, x, phi, dphi, ddphi)*star->area_element(x)*dxc(0, &xc);
		};

		Container densities = evaluate_function(density_function, 3, 0.0);
		densities.set_doubling_begin(0, {false, false, false});
		densities.set_doubling_begin(1, {true, true, true});
		Evaluator densities_eval(densities);

		//Calculating mass, momentum and energy
		std::vector<double> charges = densities_eval.integrate();
		star->params["komar_mass"] = charges[0];
		star->params["komar_momentum"] = charges[1];
		star->params["adm_energy"] = charges[2];

		//Printing out masses and radii
		std::cout << "Komar mass is: " << star->params["komar_mass"] << std::endl;
		std::cout << "Komar angular momentum is: " << star->params["komar_momentum"] << std::endl;
		std::cout << "ADM energy is: " << star->params["adm_energy"] << std::endl;

		//Calculating maximal scalar field value
		struct{
			double xmax;
			double ymax;
			double phi_max = 0.0;

			double operator()(int fn, double *x, double *phi, double **dphi, double **ddphi){
				if(phi[4] > phi_max){
					xmax = x[0];
					ymax = x[1];
					phi_max = phi[4];
				}
				return -phi[4];
			}
		} neg_phi;

		Container sc_field = evaluate_function(neg_phi,1,0.0);
		sc_field.set_doubling_begin(0, {true});
		sc_field.set_doubling_begin(1, {true});
		Evaluator sc_field_eval(sc_field);
		Collocation_Nelder_Mead sc_field_minimizer(sc_field_eval);
		sc_field_minimizer.set_simplex({neg_phi.xmax, neg_phi.ymax});
		star->params["phi_max"] = -sc_field_minimizer.minimize();

		std::cout << "Maximal value of the scalar field is: " << star->params["phi_max"] << std::endl;

		//Calculating minimal g_tt value
		struct{
			double xmin;
			double ymin;
			double gtt_min = 1.0;

			double operator()(int fn, double *x, double *phi, double **dphi, double **ddphi){
				double tmp = x[0]*std::sin(x[1])*phi[3];
				tmp = std::exp(phi[0]+phi[1]) - std::exp(phi[0]-phi[1])*tmp*tmp;
				if(tmp < gtt_min){
					xmin = x[0];
					ymin = x[1];
					gtt_min = tmp;
				}
				return tmp;
			}
		} g_tt_function;

		Container g_tt = evaluate_function(g_tt_function,1,1.0);
		g_tt.set_doubling_begin(0, {false});
		g_tt.set_doubling_begin(1, {false});
		Evaluator g_tt_eval(g_tt);
		Collocation_Nelder_Mead g_tt_minimizer(g_tt_eval);
		g_tt_minimizer.set_simplex({g_tt_function.xmin, g_tt_function.ymin});
		star->params["min_gtt"] = g_tt_minimizer.minimize();
		double gtt_min_x = g_tt_minimizer.min_x[0];

		std::cout << "Minimal value of g_tt component " << star->params["min_gtt"] << " reached at r = " << xc(0, &gtt_min_x) << std::endl;
		if(star->params["min_gtt"] < 0)
			std::cout << "The star has an ergoregion" << std::endl;
		else
			std::cout << "The star does not have an ergoregion" << std::endl;
	}

	int solve(void){
		std::cout << "Solving star equations" << std::endl;
		int err = Collocation_Solver_Full<Radial_Compactification, Rotating_Boson_Star_Equations>::solve();
		if(err) return err;

		//Calculating maximal residual values
		star->residual_error = find_function_maximum(&Rotating_Boson_Star_Equations::residuals, 2);
		std::cout << "First residual maximum is: " << star->residual_error[0] << std::endl;
		std::cout << "Second residual maximum is: " << star->residual_error[1] << std::endl;

		find_star_parameters();
		return 0;
	}

	// void plot_polar_cut(const char *folder_name, int nr, int nz);

	friend class Extrapolation;
};

Boson_Star gaussian_seed_initial_configuration(const std::string &parameter_file, const Domain_Partitioning &part, bool verbose = false);

#endif