#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include <boson_star.h>

int main(int argc, char const *argv[]){
	std::string setup_parameters_file = "setup_data.dat";
	std::string setup_domains_file = "setup_domains.dat";
	std::string folder_name = "../test";
	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-setup_params") == 0){
			setup_parameters_file = argv[i+1];
			i++;
		}
		if(strcmp(argv[i], "-setup_domains") == 0){
			setup_domains_file = argv[i+1];
			i++;
		}
		if(strcmp(argv[i], "-out") == 0){
			folder_name = argv[i+1];
			i++;
		}
	}

	std::cout << "Running initial guess solver\n" << std::endl;
	Boson_Star star = gaussian_seed_initial_configuration(setup_parameters_file, Domain_Partitioning(setup_domains_file), true);
	Boson_Star_Solver solver(&star);

	solver.solve();

	if(!folder_name.empty()){
		solver.star_save(folder_name);
	}

	return 0;
}