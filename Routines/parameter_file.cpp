#include <parameter_file.h>

Parameter_File::Parameter_File(const std::string &input): num(0), sz(0), pfile_name("default"){
	std::istringstream file(input);
	std::string line; 
	while(std::getline(file, line)){
		std::istringstream iss(line);
		if(iss.peek() != EOF){
			std::string token;
			std::getline(iss, token, ':');
			if(iss.peek() != EOF){ //Treating fields with ranged parameters
				if(token == "Name"){
					std::getline(iss, token, '#');
					while(token[0] == ' ') token.erase(0, 1);
					pfile_name = token;
					continue;
				}
				size_t index_st = token.find("%");
				size_t index_end = token.find(" ", index_st);
				if(index_st != std::string::npos && token[index_end-1] == 'd'){
					std::string token_tmp;
					std::string::size_type sztp;
					int ind_st = 0;//Getting index start and end
					token_tmp = token.substr(index_st+1, index_end-index_st-2);
					if(token_tmp.size() != 0) ind_st = std::stoi(token_tmp, &sztp);
					int st, end;//Getting range start and end
					std::getline(iss, token_tmp, '#');
					std::istringstream iss_range(token_tmp);
					std::getline(iss_range, token_tmp, '-');
					st = std::stoi(token_tmp, &sztp); //Finding range start
					if(iss_range.peek() == EOF){ //Finding range end
						end = st;
					}
					else{
						std::getline(iss_range, token_tmp, '#');
						end = std::stoi(token_tmp, &sztp);
						if(end < st) end = st;
					}
					int replace_size = index_end-index_st; //Generating entries
					for(int i = 0; i <= end-st; i++){
						token_tmp = std::to_string(ind_st+i);
						token.replace(index_st, replace_size, token_tmp);
						replace_size = token_tmp.length();
						fields.push_back(token);
						places.push_back(st+i);
						num++;
					}
				}
				else{
					fields.push_back(token);
					std::getline(iss, token, '#');
					std::string::size_type sztp;
					places.push_back(std::stoi(token, &sztp));
					num++;
				}
			}
		}
	}

	for(int i = 0; i < places.size(); i++) if(sz < places[i]) sz = places[i];
	values.assign(++sz, 0.0);
}

Parameter_File::Parameter_File(const char *input, int input_size): num(0), sz(0), pfile_name("default"){
	std::string str(input, input_size);
	std::istringstream file(str);
	std::string line; 
	while(std::getline(file, line)){
		std::istringstream iss(line);
		if(iss.peek() != EOF){
			std::string token;
			std::getline(iss, token, ':');
			if(iss.peek() != EOF){ //Treating fields with ranged parameters
				if(token == "Name"){
					std::getline(iss, token, '#');
					while(token[0] == ' ') token.erase(0, 1);
					pfile_name = token;
					continue;
				}
				size_t index_st = token.find("%");
				size_t index_end = token.find(" ", index_st);
				if(index_st != std::string::npos && token[index_end-1] == 'd'){
					std::string token_tmp;
					std::string::size_type sztp;
					int ind_st = 0;//Getting index start and end
					token_tmp = token.substr(index_st+1, index_end-index_st-2);
					if(token_tmp.size() != 0) ind_st = std::stoi(token_tmp, &sztp);
					int st, end;//Getting range start and end
					std::getline(iss, token_tmp, '#');
					std::istringstream iss_range(token_tmp);
					std::getline(iss_range, token_tmp, '-');
					st = std::stoi(token_tmp, &sztp); //Finding range start
					if(iss_range.peek() == EOF){ //Finding range end
						end = st;
					}
					else{
						std::getline(iss_range, token_tmp, '#');
						end = std::stoi(token_tmp, &sztp);
						if(end < st) end = st;
					}
					int replace_size = index_end-index_st; //Generating entries
					for(int i = 0; i <= end-st; i++){
						token_tmp = std::to_string(ind_st+i);
						token.replace(index_st, replace_size, token_tmp);
						replace_size = token_tmp.length();
						fields.push_back(token);
						places.push_back(st+i);
						num++;
					}
				}
				else{
					fields.push_back(token);
					std::getline(iss, token, '#');
					std::string::size_type sztp;
					places.push_back(std::stoi(token, &sztp));
					num++;
				}
			}
		}
	}

	for(int i = 0; i < places.size(); i++) if(sz < places[i]) sz = places[i];
	values.assign(++sz, 0.0);
}

Parameter_File::Parameter_File(const char *file_name): num(0), sz(0), pfile_name("default"){
	std::ifstream file(file_name);
	if(!file){
		std::cout << "No input file in constructor\n";
		throw("No input file in constructor");
	}
	std::string line; 
	while(std::getline(file, line)){
		std::istringstream iss(line);
		if(iss.peek() != EOF){
			std::string token;
			std::getline(iss, token, ':');
			if(iss.peek() != EOF){ //Treating fields with ranged parameters
				if(token == "Name"){
					std::getline(iss, token, '#');
					while(token[0] == ' ') token.erase(0, 1);
					pfile_name = token;
					continue;
				}
				size_t index_st = token.find("%");
				size_t index_end = token.find(" ", index_st);
				if(index_st != std::string::npos && token[index_end-1] == 'd'){
					std::string token_tmp;
					std::string::size_type sztp;
					int ind_st = 0;//Getting index start and end
					token_tmp = token.substr(index_st+1, index_end-index_st-2);
					if(token_tmp.size() != 0) ind_st = std::stoi(token_tmp, &sztp);
					int st, end;//Getting range start and end
					std::getline(iss, token_tmp, '#');
					std::istringstream iss_range(token_tmp);
					std::getline(iss_range, token_tmp, '-');
					st = std::stoi(token_tmp, &sztp); //Finding range start
					if(iss_range.peek() == EOF){ //Finding range end
						end = st;
					}
					else{
						std::getline(iss_range, token_tmp, '#');
						end = std::stoi(token_tmp, &sztp);
						if(end < st) end = st;
					}
					int replace_size = index_end-index_st; //Generating entries
					for(int i = 0; i <= end-st; i++){
						token_tmp = std::to_string(ind_st+i);
						token.replace(index_st, replace_size, token_tmp);
						replace_size = token_tmp.length();
						fields.push_back(token);
						places.push_back(st+i);
						num++;
					}
				}
				else{
					fields.push_back(token);
					std::getline(iss, token, '#');
					std::string::size_type sztp;
					places.push_back(std::stoi(token, &sztp));
					num++;
				}
			}
		}
	}

	for(int i = 0; i < places.size(); i++) if(sz < places[i]) sz = places[i];
	values.assign(++sz, 0.0);
}

void Parameter_File::read_file(const std::string &file_name){
	std::ifstream file(file_name);
	if(!file){
		std::cout << "No input file in read_file\n";
		throw("No input file in read_file");
	}
	std::string line; 
	while(std::getline(file, line)){
		std::istringstream iss(line);
		if(iss.peek() != EOF){
			std::string token;
			std::getline(iss, token, ':');
			if(iss.peek() != EOF){
				for(int i = 0; i < num; i++)
					if(token == fields[i]){
						std::getline(iss, token, ':');
						std::string::size_type sztp;
						values[places[i]] = std::stod(token, &sztp);
						break;
					}
			}
		}
	}
	return;
}

void Parameter_File::print_file(const std::string &file_name) const{
	std::ofstream out;
	out.open(file_name);

	out << std::setprecision(15);
	for(int i = 0; i < num; i++)
		out << fields[i] << ": " << values[places[i]] << std::endl;

	out.close();
	return;
}