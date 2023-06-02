#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include <timer.h>
#include <boson_star.h>
#include <supplement.h>

int main(int argc, char const *argv[]){
	std::string folder_in;
	std::string folder_out;
	std::string domains_file;

	bool change_sigma = false;
	double sigma;

	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-in") == 0){
			folder_in = argv[i+1];
			i++;
		}
		else if(strcmp(argv[i], "-out") == 0){
			folder_out = argv[i+1];
			i++;
		}
		else if(strcmp(argv[i], "-domains") == 0){
			domains_file = argv[i+1];
			i++;
		}
		else if(strcmp(argv[i], "-sigma") == 0){
			change_sigma = true;
			sigma = atof(argv[i+1]);
			i++;
		}
		else{
			std::cout << "Unrecognized flag: ";
			printf("%s\n", argv[i]);
			return -1;
		}
	}

	if(folder_in.empty()){
		std::cout << "No input file provided" << std::endl;
		return 0;
	}

	Boson_Star star(folder_in);
	if(!domains_file.empty()){
		Domain_Partitioning part = read_partitioning(domains_file);
		// Domain_Partitioning part(domains_file);

		std::cout << "New domains configuration:\n" << part << std::endl;

		star.regrid(part);
		std::cout << "New solution size: " << star.size() << '\n' << std::endl;
	}

	Boson_Star_Solver solver(&star);

	struct{
		std::string name = "sigma";

		double calculate_parameter(const Boson_Star *star){
			return std::sqrt(-4.0*star->poly_coef[0]/star->poly_coef[1]);
		}

		void update_parameter(Boson_Star *star, double sigma){
			double mass = std::sqrt(star->poly_coef[0]);
			star->poly_coef[1] = -4.0*mass*mass/(sigma*sigma);
			star->poly_coef[2] = 4.0*mass*mass/std::pow(sigma, 4.0);
		}
	} sigma_calculator;

	Timer timer;

	if(change_sigma)
		solver.parameter_change(sigma_calculator, sigma);
	else
		solver.solve();

	timer.record();
	std::cout << "The calculation took " << timer.last_period << std::endl;

	if(!folder_out.empty()){
		solver.star_save(folder_out);
		solver.print_domain_2D(folder_out + "/domains.pbm");
	}

	return 0;
}