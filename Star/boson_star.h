#include <container.h>
#include <parameter_file.h>
// #include <collocation_routines.h>
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

	int solve(void){
		std::cout << "Solving star equations" << std::endl;
		int err = Collocation_Solver_Full<Radial_Compactification, Rotating_Boson_Star_Equations>::solve();
		if(err) return err;

		//Calculating maximal residual values
		star->residual_error = find_function_maximum(&Rotating_Boson_Star_Equations::residuals, 2);
		std::cout << "First residual maximum is: " << star->residual_error[0] << std::endl;
		std::cout << "Second residual maximum is: " << star->residual_error[1] << std::endl;

		// find_star_parameters();
		return 0;
	}

	// void plot_polar_cut(const char *folder_name, int nr, int nz);

	friend class Extrapolation;
};

Boson_Star gaussian_seed_initial_configuration(const std::string &parameter_file, const Domain_Partitioning &part, bool verbose = false);

#endif