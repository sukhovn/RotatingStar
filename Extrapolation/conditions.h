struct: public Extrapolation_Condition{
	double previous_local_maximum_observed = 0.0;
	double global_maximum = 0.0;

	void check_previous_max_save(const std::string &max_path){
		struct stat sttus = {0};
		std::string::size_type sztp;
		std::string line, token;
		if(stat(max_path.c_str(), &sttus) == -1) return;
		std::ifstream inp_file(max_path + "/parameters.dat");
		if(!inp_file.is_open()){
			std::cout << "No parameters file found in previous maximal mass solution" << std::endl;
			throw("reading error");
		}
		while(std::getline(inp_file, line)){
			std::istringstream iss(line);
			std::getline(iss, token, ':');
			if(token.compare("The Komar mass is") == 0){
				std::getline(iss, token, ':');
				global_maximum = std::stod(token, &sztp);
				break;
			}
		}
	}

	void check_condition(Extrapolation &extra){
		if(previous_local_maximum_observed > extra.sol_current->params.at("komar_mass")){
			if(previous_local_maximum_observed > global_maximum){
				global_maximum = previous_local_maximum_observed;
				std::cout << "Local mass maximum reached at " << previous_local_maximum_observed << std::endl;
				extra.save_solution(extra.save_folder + "/Maximal_mass",1);
			}
		}
		else
			previous_local_maximum_observed = extra.sol_current->params.at("komar_mass");
	}

	void final_action(Extrapolation &extra, bool normal_exit){
		if(previous_local_maximum_observed > global_maximum){
			global_maximum = previous_local_maximum_observed;
			std::cout << "A maximum in mass was reached at the end: " << global_maximum << std::endl;
			extra.save_current_solution(extra.save_folder + "/Maximal_mass");
		}
		else{
			std::cout << "A maximum in mass was reached at: " << global_maximum << '\n' << std::endl;
		}
	}

} check_maximum_mass;

struct: public Extrapolation_Condition{
	std::vector<std::vector<double>> residual_max = std::vector<std::vector<double>>(2, std::vector<double>(2, 0.0));
	std::vector<double> residual_max_ps_length = std::vector<double>(2, 0.0);

	void check_condition(Extrapolation &extra){
		for(int i = 0; i < 2; ++i){
			if(extra.sol_current->residual_error[i] > residual_max[i][i]){
				residual_max[i] = extra.sol_current->residual_error;
				residual_max_ps_length[i] = extra.sol_current->ps_length;
			}
		}
	}

	void final_action(Extrapolation &extra, bool normal_exit){
		std::cout << "The first residual reached maximum at pseudolength: " << residual_max_ps_length[0] << std::endl;
		std::cout << "At that point the residuals maxima were:" << std::endl;
		std::cout << "First residual maximum is: " << residual_max[0][0] << std::endl;
		std::cout << "Second residual maximum is: " << residual_max[0][1] << '\n' << std::endl;

		std::cout << "The second residual reached maximum at pseudolength: " << residual_max_ps_length[1] << std::endl;
		std::cout << "At that point the residuals maxima were:" << std::endl;
		std::cout << "First residual maximum is: " << residual_max[1][0] << std::endl;
		std::cout << "Second residual maximum is: " << residual_max[1][1] << '\n' << std::endl;
	}

} maximal_residual;