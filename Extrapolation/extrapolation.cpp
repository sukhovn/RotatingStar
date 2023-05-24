#include <iostream>
#include <iomanip>
#include <fstream>
#include <new>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <interpolation.h>
#include <extrapolation.h>

#define ZERO_CHECK(a) (a < 1.0e-5 ? 1.0e-5 : a)
#define MAX_FAILS 5

//Here we define parameter recorder

//This routine initializes a structure for parameter recorder, it creates or opens a folder with files with phi, w, msmass and rad
//If append != 0 it appends solutions to already existing files
Solution_Saver::Solution_Saver(const std::string &file_name, bool append) : is_on(false), records(0),
	save_parameters({"omega", "phi_max", "komar_mass", "komar_momentum", "radius_charge", "radius_komar", "min_gtt"})
{	
	//If append == 0 we are overwriting old files, if not we are appending new data to the old files

	std::string head_line = "ps_length";
	for(auto par : save_parameters)
		head_line = head_line + ',' + par;

	//Creating or opening output files
	if(append){
		std::ifstream f_test(file_name);
		if(f_test.good() == false) append = false;
		else{
			std::string line;
			std::getline(f_test, line);
			if(line.compare(head_line) != 0){
				std::cout << "Table save parameters mismatch" << std::endl;
				throw("append error");
			}
			f_test.close();
		}
	}

	if(append){
		out.open(file_name, std::ofstream::app);
		if(!out.is_open()){
			std::cout << "Solution saver file could not be opened" << std::endl;
			throw("solution saver file open error");
		}
	}
	else{
		out.open(file_name);
		if(!out.is_open()){
			std::cout << "Solution saver file could not be opened" << std::endl;
			throw("solution saver file open error");
		}
		out << head_line << std::endl;
	}

	out << std::setprecision(16) << std::fixed;
}

//This routine stores solution for parameter recorder
Solution_Saver& Solution_Saver::operator<< (const Boson_Star& star){
	if(is_on){
		std::cout << "Recording star solution\n";
		out << star.ps_length;
		for(std::string par : save_parameters)
			out << ',' << star.params.at(par);
		out << std::endl;
		records++;}
	return *this;
}


Extrapolation::Extrapolation(const Boson_Star &sol, int num, const std::string &_save_folder):
	Boson_Star_Solver(new Boson_Star(sol)), sol_num(num), is_initialized(false), saver(nullptr),
	ps_length_components({"omega", "komar_mass", "phi_max"}),
	save_folder(_save_folder)
{
	sol_current = star;
	sol_array = new Boson_Star*[num];

	stepper_increment = 0.01;
}

int Extrapolation::initial_fill(const double dpar){
	if(is_initialized){
		std::cout << "The extrapolation has already been initialized, aborting" << std::endl;
		return 1;
	}

	double mult_coef_init;

	if(dpar > 0.0){
		mult_coef_init = 1.0 + stepper_increment;
	}
	else{
		mult_coef_init = 1.0 - stepper_increment;
	}

	double omega_past = sol_current->omega, par_last = -1.0e5;
	int ord = 1, sol_record = 0;

	for(int i = 0; i < 70; ++i)
		std::cout << '*';
	std::cout << std::endl;
	std::cout << "Initializing extrapolation" << std::endl;;

	while(1){
		for(int i = 0; i < 70; ++i)
			std::cout << '-';
		std::cout << std::endl;
		std::cout << "Solving for omega = " << sol_current->omega << "\n";
		if(solve()){
			for(int i = 0; i < 70; ++i)
				std::cout << '-';
			std::cout << std::endl;
			mult_coef_init = std::sqrt(mult_coef_init);
			std::cout << "Changing multiplication coefficient to " << mult_coef_init << "\n";
			ord++;
		}
		else{
			omega_past = sol_current->omega;
			if(fabs(sol_current->ps_length - par_last) > fabs(dpar)){
				if(sol_record == 0)
					sol_current->ps_length = 0.0;
				else{
					double tmp, dlen = 0.0;
					for(auto par : ps_length_components){
						tmp = (sol_current->params.at(par) - sol_array[sol_record-1]->params.at(par))/ZERO_CHECK(sol_array[sol_record-1]->params.at(par));
						dlen += tmp*tmp;
					}
					sol_current->ps_length = sol_array[sol_record-1]->ps_length + std::sqrt(dlen);
				}

				std::cout << "Recording the solution " << sol_record+1 << " out of " << sol_num << " with pseudolengh " << sol_current->ps_length << " for extrapolation\n";
				sol_array[sol_record] = new Boson_Star(*sol_current);

				for(int i = 0; i < 70; ++i)
					std::cout << '-';
				std::cout << std::endl;

				sol_record++;

				if(sol_record > sol_num-1)
					break;
			}
		}
		if(ord > 20){
			printf("Minimal increment in omega reached, terminating\n");
			return 1;
			for(int i = 0; i < sol_record; i++) delete sol_array[i];
		}
		sol_current->omega = omega_past * mult_coef_init;
	}

	std::cout << "Extrapolation initialized\n";
	for(int i = 0; i < 70; ++i)
		std::cout << '*';
	std::cout << std::endl;

	//Setting up tracing coefficient
	mult_coef = mult_coef_max = 1.0 + stepper_increment;
	fails = 1;
	easy_steps = 0;

	//Creating output folder
	struct stat sttus = {0};
	if(stat(save_folder.c_str(), &sttus) == -1){
		mkdir(save_folder.c_str(), 0700);
	}

	//Setting up parameter saver
	saver = new Solution_Saver(save_folder + "/tracing.csv", false);

	is_initialized = true;

	return 0;
}

void Extrapolation::load_tracing_state(const std::string &folder_name){
	if(is_initialized){
		std::cout << "The extrapolation has already been initialized, aborting" << std::endl;
		return;
	}

	for(int i = 0; i < 70; ++i)
		std::cout << '*';
	std::cout << std::endl;
	std::cout << "Loading extrapolator state from folder: " + folder_name << std::endl;

	//Counting number of saved solutions
	std::string path;
	int nfunc_load = 0;
	struct stat sttus = {0};
	while(1){
		path = save_folder + '/' + folder_name + "/Solution_" + std::to_string(nfunc_load);
		if(stat(path.c_str(), &sttus) == 0) nfunc_load++;
		else break;
	}
	if(nfunc_load == 0){
		std::cout << "No saved solutions found in folder: " << folder_name << std::endl;
		throw("reading error");
	}
	std::cout << "Found " << nfunc_load << " solutions" << std::endl;
	delete [] sol_array;
	sol_num = nfunc_load;
	sol_array = new Boson_Star*[sol_num];

	//Reading solutions
	std::ifstream inp_file;
	std::string::size_type sztp;
	std::string line, token;
	for(int i = 0; i < sol_num; ++i){
		path = save_folder + '/' + folder_name + "/Solution_" + std::to_string(i);
		sol_array[i] = new Boson_Star(path);
		sol_array[i]->ps_length = -1.0;
		//Reading pseudolengths
		inp_file.open(path + "/parameters.dat");
		while(std::getline(inp_file, line)){
			std::istringstream iss(line);
			std::getline(iss, token, ':');
			if(token.compare("The pseudolength is") == 0){
				std::getline(iss, token, ':');
				sol_array[i]->ps_length = std::stod(token, &sztp);
				break;
			}
		}
		inp_file.close();
		if(sol_array[i]->ps_length == -1.0){
			std::cout << "Pseudolength for solution " << i << " was not read" << std::endl;
			throw("reading error");
		}
	}

	delete sol_current;
	sol_current = new Boson_Star(*sol_array[sol_num-1]);
	change_star(sol_current);
	adjust_compactification(sol_current->partition, true);
	std::cout << "Continuing tracing from pseudolength: " << sol_current->ps_length << std::endl;

	//Loading multiplication coefficient
	mult_coef = 0.0;
	fails = -1;
	easy_steps = -1;
	inp_file.open(save_folder + '/' + folder_name + "/save_parameters.dat");
	if(!inp_file.is_open()){
		std::cout << "Save parameters file not found" << std::endl;
		throw("reading error");
	}
	while(std::getline(inp_file, line)){
		std::istringstream iss(line);
		std::getline(iss, token, ':');
		if(token.compare("The multiplication coefficient is") == 0){
			std::getline(iss, token, ':');
			mult_coef = std::stod(token, &sztp);
		}
		if(token.compare("The fails parameter is") == 0){
			std::getline(iss, token, ':');
			fails = std::stoi(token, &sztp);
		}
		if(token.compare("The easy steps parameter is") == 0){
			std::getline(iss, token, ':');
			easy_steps = std::stoi(token, &sztp);
		}
	}
	if(mult_coef == 0.0 || fails == -1 || easy_steps == -1){
		std::cout << "The extrapolation parameters were not loaded" << std::endl;
		throw("reading error");
	}
	else
		std::cout << "Loading multiplication coefficient: " << mult_coef << std::endl;
	inp_file.close();

	//Preparing parameter recorder file
	inp_file.open(save_folder + "/tracing.csv");
	std::ofstream out_file(save_folder + "/tracing_cont.csv");
	double tmp;
	int line_read = 0, values_found = 0, prev_found;
	while(std::getline(inp_file, line)){
		++line_read;
		out_file << line << std::endl;
		//Checking if line corresponds to one of the loaded solutions
		std::istringstream iss(line);
		std::getline(iss, token, ',');
		try{
			tmp = std::stod(token, &sztp);
		}
		catch(std::invalid_argument const& ex){
			continue;
		}
		if(std::fabs(sol_array[values_found]->ps_length - tmp) < 1.0e-12){
			if(values_found && prev_found != line_read-1){
				std::cout << "Solutions " << values_found-1 << " and " << values_found << " are not continuous in the recorder file" << std::endl;
				throw("reading error");
			}
			else{
				prev_found = line_read;
				++values_found;
			}
			if(values_found == sol_num)
				break;
		}
	}

	if(values_found != sol_num){
		std::cout << "Not enough solutions found in the recorder file" << std::endl;
		throw("reading error");
	}

	out_file.close();
	inp_file.close();
	saver = new Solution_Saver(save_folder + "/tracing_cont.csv", true);

	for(int i = 0; i < 70; ++i)
		std::cout << '*';
	std::cout << std::endl;

	is_initialized = true;
}

int Extrapolation::extrapolate_solution(double *param_array, const double param, bool swap_solutions){
	if(!is_initialized){
		std::cout << "The extrapolation has not been initialized, aborting\n";
		return 1;
	}
	std::cout << "Attempting extrapolation to pseudolength " << param << "\n";

	double tmp, dyi;
	Boson_Star *sol_pt;
	std::vector<double> tmparr(sol_num);

	std::cout << "Extrapolating solution\n";

	for(int i = 0; i < sol_num; i++)
		tmparr[i] = sol_array[i]->omega;
	if(polint(param_array, tmparr.data(), sol_num, param, &tmp, &dyi) != 0){
		std::cout << "Could not extrapolate omega, extrapolation failed\n";
		return 1;
	}

	sol_current->omega = tmp;
	std::cout << "Extrapolated omega value is " << tmp << "\n";

	for(int i = 0; i < sol_current->size(); i++){
		for(int ii = 0; ii < sol_num; ii++)
			tmparr[ii] = (*sol_array[ii])[i];
		if(polint(param_array, tmparr.data(), sol_num, param, &tmp, &dyi) != 0){
			std::cout << "Could not extrapolate solution, extrapolation failed\n";
			return 1;
		}
		(*sol_current)[i] = tmp;
	}

	std::cout << "Extrapolation is done, getting a solution from extrapolation\n";
	if(solve() != 0){
		std::cout << "Could not get a solution from extrapolation\n";
		sol_current->Boson_Star::operator=(*sol_array[sol_num-1]);
		return 1;
	}

	double dlen = 0.0;
	for(auto par : ps_length_components){
		tmp = (sol_current->params.at(par) - sol_array[sol_num-1]->params.at(par))/ZERO_CHECK(sol_array[sol_num-1]->params.at(par));
		dlen += tmp*tmp;
	}
	sol_current->ps_length = sol_array[sol_num-1]->ps_length + std::sqrt(dlen);
	std::cout << "The pseudolength is " << sol_current->ps_length << "\n";

	//Interchanging solutions
	if(swap_solutions){
		delete sol_array[0];

		for(int i = 0; i < sol_num-1; i++){
			sol_array[i] = sol_array[i+1];
		}

		sol_array[sol_num-1] = new Boson_Star(*sol_current);
	}

	return 0;
}

int Extrapolation::regrid_solutions(const Domain_Partitioning &part_new){
	if(!is_initialized){
		std::cout << "The extrapolation has not been initialized, aborting\n";
		return 1;
	}

	for(int i = 0; i < 70; ++i)
		std::cout << '*';
	std::cout << std::endl;

	std::cout << "New domains configuration:" << std::endl;
	for(const Bulk_Domain *dom = part_new.first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
		std::cout << "Domain: " << dom->bulk_domain_number << std::endl;
		std::cout << *dom << std::endl;
	}
	adjust_compactification(part_new, true);

	Boson_Star **sol_backup = new Boson_Star*[sol_num];
	for(int i = 0; i < sol_num; ++i){
		std::cout << "\nRegridding solution: " << i << std::endl;
		sol_backup[i] = new Boson_Star(*sol_array[i]);
		sol_array[i]->regrid(part_new);
		change_star(sol_array[i]);
		if(solve() != 0){
			std::cout << "Regridding failed" << std::endl;
			for(int j = 0; j <= i; ++j){
				delete sol_array[j];
				sol_array[j] = sol_backup[j];
			}
			change_star(sol_current);
			adjust_compactification(sol_current->partition, false);
			delete [] sol_backup;
			return 1;
		}
	}
	std::cout << "\nRegridding was successful" << std::endl;
	for(int i = 0; i < 70; ++i)
		std::cout << '*';
	std::cout << std::endl;

	delete sol_current;
	sol_current = new Boson_Star(*sol_array[sol_num-1]);
	change_star(sol_current);
	return 0;
}

int Extrapolation::param_change(double param_new){
	if(!is_initialized){
		std::cout << "The extrapolation has not been initialized, aborting\n";
		return 1;
	}

	if(param_new <= sol_current->ps_length){
		std::cout << "Required pseudolength is reached, no tracing is needed" << std::endl;
		return 0;
	}

	double param_last = sol_current->ps_length, param_try;
	std::vector<double> param_array(sol_num);

	while(1){
		param_try = param_last * mult_coef;
		for(int i = 0; i < sol_num; i++)
			param_array[i] = sol_array[i]->ps_length;

		for(int i = 0; i < 70; ++i)
			std::cout << '-';
		std::cout << std::endl;
		if(extrapolate_solution(param_array.data(), param_try, true) != 0){
			for(int i = 0; i < 70; ++i)
				std::cout << '-';
			std::cout << std::endl;

			for(int i = 0; i < 70; ++i)
				std::cout << '+';
			std::cout << std::endl;
			mult_coef = std::sqrt(mult_coef);
			std::cout << "Changing multiplication coefficient to " << mult_coef << "\n";
			fails++;
			for(int i = 0; i < 70; ++i)
				std::cout << '+';
			std::cout << std::endl;
		}
		else{
			*saver << *sol_current;
			for(int i = 0; i < 70; ++i)
				std::cout << '-';
			std::cout << std::endl;

			if(iterations > 2){
				mult_coef = std::sqrt(mult_coef);
				std::cout << "Too many iterations, changing multiplication coefficient to " << mult_coef << "\n";
			}
			param_last = sol_current->ps_length;

			if(iterations == 1) easy_steps++;
			else easy_steps = 0;

			for(auto cond = conditions.begin(); cond != conditions.end(); ++cond)
				(*cond)->check_condition(*this);
		}

		//No further extrapolation termination
		if(fails > MAX_FAILS){
			for(int i = 0; i < 70; ++i)
				std::cout << '+';
			std::cout << std::endl;
			std::cout << "Minimal increment is reached, terminating\n";
			for(int i = 0; i < 70; ++i)
				std::cout << '+';
			std::cout << std::endl;
			return 1;
		}

		//If there were too many easy steps we increase the stepsize
		if(easy_steps > 10){
			if(mult_coef < mult_coef_max - 1.0e-6){
				mult_coef = std::pow(mult_coef, 2.0);
				std::cout << "Increasing interpolation step, changing multiplication coefficient to " << mult_coef << std::endl;
			}
			easy_steps = 0;
			fails = 0;
		}

		//Termination condition
		if(param_new - param_last * mult_coef < 0.0){
			if(param_new - param_last < 0.0){
				break;
			}
			param_try = param_new;
			for(int i = 0; i < sol_num; i++)
				param_array[i] = sol_array[i]->ps_length;

			for(int i = 0; i < 70; ++i)
				std::cout << '-';
			std::cout << std::endl;
			if(extrapolate_solution(param_array.data(), param_try, false) == 0){
				*saver << *sol_current;
				for(int i = 0; i < 70; ++i)
					std::cout << '-';
				std::cout << std::endl;
				break;
			}
			else{
				for(int i = 0; i < 70; ++i)
					std::cout << '-';
				std::cout << std::endl;

				fails++;
				while(param_new - (param_last * mult_coef) < 0.0)
					mult_coef = std::sqrt(mult_coef);
			}
		}
	}

	std::cout << "The extrapolation is done, omega is: " << sol_current->omega << "\n";
	return 0;
}

int Extrapolation::pslenght_tracer(double pslength_end){
	if(!is_initialized){
		std::cout << "The extrapolation has not been initialized, aborting\n";
		return 1;
	}

	saver->on();
	std::cout << "Tracing pslength\n";
	std::cout << "Maximum multiplication coefficient set to: " << mult_coef_max << std::endl;
	if(param_change(pslength_end)){
		std::cout << "\nTracing interrupted" << std::endl;
		for(auto cond = conditions.begin(); cond != conditions.end(); ++cond)
			(*cond)->final_action(*this, false);
		return 1;
	}

	std::cout << "\nTracing is finished" << std::endl;

	for(auto cond = conditions.begin(); cond != conditions.end(); ++cond)
		(*cond)->final_action(*this, true);

	return 0;
}

void Extrapolation::save_tracing_state(const std::string &folder_name){
	if(!is_initialized){
		std::cout << "The extrapolation has not been initialized, aborting\n";
		return;
	}

	std::cout << "Saving extrapolator state to folder: " + folder_name << std::endl;
	std::string folder_path = save_folder + '/' + folder_name;
	//Creating output folder
	struct stat sttus = {0};
	if(stat(folder_path.c_str(), &sttus) == -1){
		mkdir(folder_path.c_str(), 0700);
	}

	//Saving solutions
	for(int i = 0; i < sol_num; ++i){
		save_solution(folder_path + "/Solution_" + std::to_string(sol_num-i-1),i);
	}

	//Saving pseudolength to parameter files
	std::ofstream param_file;
	for(int i = 0; i < sol_num; ++i){
		param_file.open(folder_path + "/Solution_" + std::to_string(i) + "/parameters.dat", std::ofstream::app);
		param_file << std::setprecision(16);
		param_file << "The pseudolength is: " << sol_array[i]->ps_length << std::endl;
		param_file.close();
	}

	//Saving multiplication constant
	param_file.open(folder_path + "/save_parameters.dat");
	param_file << std::setprecision(16);
	param_file << "The multiplication coefficient is: " << mult_coef << std::endl;
	param_file << "The fails parameter is: " << fails << std::endl;
	param_file << "The easy steps parameter is: " << easy_steps << std::endl;
	param_file.close();
}