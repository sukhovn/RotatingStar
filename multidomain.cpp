#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include <boson_star.h>
// #include <extrapolation.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void star_test(int argc, char const *argv[]){
	std::string folder_in;
	std::string folder_out;
	std::string domains_file;
	// bool change_sigma = false;
	// bool change_m = false;
	// bool change_omega = false;	
	// int plsize = 51;

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
		// else if(strcmp(argv[i], "-plot_size") == 0){
		// 	plsize = std::stoi(argv[i+1]);
		// 	i++;
		// }
		else{
			std::cout << "Unrecognized flag: ";
			printf("%s\n", argv[i]);
			return;
		}
	}

	if(folder_in.empty()) return;
	Boson_Star star(folder_in);
	if(!domains_file.empty()){
		Domain_Partitioning part(domains_file);
		star.regrid(part);
	}

	Boson_Star_Solver solver(&star);
	if(!domains_file.empty())
		solver.solve();

	// sigma_change(solver, 0.5);
	// omega_change(solver, 0.8);
	// m_change(solver, 3);

	if(!folder_out.empty()){
		star.save(folder_out);
		// solver.plot_polar_cut(folder_out + "_plot", plsize, plsize);
	}
}


int main(int argc, char const *argv[]){
	// std::string test_folder = "../results";

	// struct stat sttus = {0};
	// if (stat(test_folder.c_str(), &sttus) == -1) {
	// 	mkdir(test_folder.c_str(), 0700);
	// }

	star_test(argc, argv);

	return 0;
}