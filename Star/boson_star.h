#include <container.h>
#include <parameter_file.h>
#include <collocation_routines.h>
#include <collocation_solver_full.h>
#include <equations.h>
#include <min_root.h>
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

class Boson_Star_Solver: public Radial_Compactification, public Collocation_Solver_Full<Radial_Compactification, Rotating_Boson_Star_Equations>{
private:
	int density_extrapolation_resolution;
	Boson_Star *star;

public:
	Boson_Star_Solver(Boson_Star *_star):
		Radial_Compactification(_star->partition),
		Collocation_Solver_Full<Radial_Compactification, Rotating_Boson_Star_Equations>(*this, _star, _star),
		star(_star),
		density_extrapolation_resolution(100) {};

	Boson_Star *return_star_ptr(){return star;};

	void change_star(Boson_Star *stari){
		star = stari;
		eq_system = stari;
		change_container(*stari);
	}

	void find_star_parameters(void);

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

	void star_save(const std::string &folder_name){
		star->save(folder_name);

		auto residuals_function = [this](int fn, double *x, double *phi, double **dphi, double **ddphi)
		{
			if(x[0] == 0.0 || x[1] == 0) return 0.0;
			return star->residuals(fn, x, phi, dphi, ddphi);
		};

		Container residuals = evaluate_function(residuals_function, 2, {0.0, 0.0});
		residuals.save(folder_name + "/residuals");

		plot_polar_cut(folder_name + "/plot", 200, 100);
	}

	void plot_polar_cut(const std::string &folder_name, int nr, int nz);

	template<class Parameter_Calculator>
	int parameter_change(Parameter_Calculator &calc, double par_new){
		double par = calc.calculate_parameter(star);
		double par_past = par;
		double mult_coef, direction;
		int ord = 1;

		if(par_new > par){
			direction = 1.0;
			mult_coef = 1.05;
		}
		else{
			direction = -1.0;
			mult_coef = 0.95;
		}

		std::cout << "Starting at multiplication coefficient " << mult_coef << std::endl;

		while(1){
			std::cout << std::endl << "Solving for " + calc.name + ": " << par << std::endl;
			if(solve() != 0){
				if(par == par_past){
					std::cout << "Star incorrectly initialized" << std::endl;
					return -1;
				}
				mult_coef = std::sqrt(mult_coef);
				par = par_past * mult_coef;
				calc.update_parameter(star, par);
				std::cout << std::endl << "Changing multiplication coefficient to " << mult_coef << std::endl;
				ord++;
			}
			else{
				par_past = par;
				par *= mult_coef;
				calc.update_parameter(star, par);
			}
			if(ord > 10){
				std::cout << std::endl << "Minimal increment in " + calc.name + " reached, terminating" << std::endl;
				return 1;
			}
			if((par_new - par) * direction < 0.0){
				par = par_new;
				calc.update_parameter(star, par);
				std::cout << std::endl << "Solving for " + calc.name + " = " << par << std::endl;
				if(solve() == 0)
					break;
				else{
					while((par_new - par_past * mult_coef) * direction < 0.0)
						mult_coef = std::sqrt(mult_coef);
					par = par_past * mult_coef;
					calc.update_parameter(star, par);
				}
			}
		}

		std::cout << std::endl << "We reached " + calc.name + " value of " << par << std::endl;
		return 0;
	}

	friend class Extrapolation;
};

Boson_Star gaussian_seed_initial_configuration(const std::string &parameter_file, const Domain_Partitioning &part, bool verbose = false);

#endif