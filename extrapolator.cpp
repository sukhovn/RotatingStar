#include <cstdio>
#include <cmath>
#include <string.h>
#include <boson_star.h>
#include <extrapolation.h>
#include <schedule.h>
#include <timer.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int main(int argc, char const *argv[]){
	double sigma = 1.0;
	double increment = 0.01;
	double ps_end = 0.5;
	std::string work_folder_name = "../Test_tracing";
	std::string load_stop, save_stop;

	for(int i = 1; i < argc; i++){
		if(strcmp(argv[i], "-wf") == 0){
			work_folder_name = argv[i+1];
			i++;
		}
		else if(strcmp(argv[i], "-ps_end") == 0){
			ps_end = atof(argv[i+1]);
			i++;
		}
		else if(strcmp(argv[i], "-increment") == 0){
			increment = atof(argv[i+1]);
			i++;
		}
		else if(strcmp(argv[i], "-load_stop") == 0){
			load_stop = argv[i+1];
			i++;
		}
		else if(strcmp(argv[i], "-save_stop") == 0){
			save_stop = argv[i+1];
			i++;
		}
		else{
			std::cout << "Wrong key option received" << std::endl;
			return -1;
		}
	}

	Boson_Star star(work_folder_name + "/Setup");

	Tracing_Schedule schedule;
	try{
		schedule.load_schedule(work_folder_name + "/schedule.dat");
	}
	catch(const char *error){
		if(strcmp(error, "file not found") == 0){
			std::cout << "No schedule.dat file provided" << std::endl;
		}
		else{
			std::cout << "Extrapolation aborted" << std::endl;
			return 1;
		}
	}

	Extrapolation extra(star, 5, work_folder_name);
	extra.stepper_increment = increment;

	//Checking if the tracing is resumed from a stop
	if(load_stop.empty()){
		//Checking scheduler entry if initial regridding is needed
		if(schedule.schedule.find(0.0) != schedule.schedule.end() && !(schedule.schedule.at(0.0) == extra.sol_current->partition)){
			Domain_Partitioning &part = schedule.schedule.at(0.0);

			std::cout << "Regridding starting solution to new domain configuration:" << std::endl;
			for(const Bulk_Domain *dom = part.first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
				std::cout << "Domain: " << dom->bulk_domain_number << std::endl;
				std::cout << *dom << std::endl;
			}

			extra.sol_current->regrid(part);
			extra.change_star(extra.sol_current);
			extra.adjust_compactification(part, true);
		}

		extra.solve();

		if(extra.initial_fill(-0.01)){
			std::cout << "Error in initial filling\n";
			return 1;
		}

		schedule.setup_scheduling(0.0);
	}
	else{
		extra.load_tracing_state(load_stop);
		schedule.setup_scheduling(extra.sol_current->ps_length);
	}

	#include <conditions.h>
	check_maximum_mass.check_previous_max_save(work_folder_name + "/Maximal_mass");
	extra.conditions.push_back(&check_maximum_mass);
	extra.conditions.push_back(&schedule);
	extra.conditions.push_back(&maximal_residual);

	Timer timer;
	extra.pslenght_tracer(ps_end);
	timer.record();
	std::cout << "The run took " << timer.last_period << std::endl;

	extra.save_current_solution(work_folder_name + "/Final");

	if(!load_stop.empty()){
		std::string old_tracing_file = work_folder_name + "/tracing.csv";
		std::string new_tracing_file = work_folder_name + "/tracing_cont.csv";
		if(std::remove(old_tracing_file.c_str())){
			std::cout << "Could not remove old tracing.csv file" << std::endl;
			throw("file removing error");
		}
		if(std::rename(new_tracing_file.c_str(), old_tracing_file.c_str())){
			std::cout << "Could not rename tracing.csv file" << std::endl;
			throw("file removing error");
		}
	}

	if(!save_stop.empty())
		extra.save_tracing_state(save_stop);

	return 0;
}