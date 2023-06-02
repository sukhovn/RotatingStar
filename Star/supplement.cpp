#include <vector>
#include <fstream>
#include <iomanip>
#include <ios>

#include <supplement.h>

Domain_Partitioning read_partitioning(const std::string &file_name){
	std::ifstream file(file_name);
    if(!file.is_open()){
		std::cout << "No domain file found" << std::endl;
		throw("read fail");
	}

	int domains_read = 0;
	std::string line, token;
	std::string::size_type sztp;
	std::vector<std::vector<double>> domains_start;
	std::vector<std::vector<double>> domains_end;
	std::vector<std::vector<int>> dim_nums;
	
	while(std::getline(file, line)){
		for(char c : line){//Skipping commented lines
			if(c == ' ') continue;
			if(c == '#') line = "";
			else break;
		}

		if(line.size() == 0)
			continue;

		std::istringstream iss(line);
		//Checking the first half of the entry
		std::getline(iss, token, ':');
		if(token.substr(0, 6).compare("Domain") != 0){
			std::cout << "Corrupt schedule file, the line should start with keyword Domain" << std::endl;
			std::cout << line << std::endl;
			throw("corrupt file");
		}
		int ist = 6, dom_read;
		while(token[ist] == ' ') ++ist;
		try{
			dom_read = std::stoi(token.substr(ist), &sztp);
		}
		catch(std::invalid_argument const& ex){
			std::cout << "Corrupt schedule file, domain number corrupted" << std::endl;
			std::cout << line << std::endl;
			throw("corrupt file");
		}

		//Checking the second half of the entry
		std::getline(iss, token, ':');
		token.erase(remove(token.begin(), token.end(), ' '), token.end());
		std::istringstream iss_data(token);
		int swcase = 0, nx, ny;
		double xst, xend, yst, yend;
		while(std::getline(iss_data, token, ',')){
			try{
				switch(swcase){
					case 0:
						if(token.substr(0, 2).compare("[[")) throw(0);
						xst = std::stod(token.substr(2));
						break;
					case 1:
						if(token.substr(token.size()-1).compare("]")) throw(1);
						xend = std::stod(token.substr(0, token.size()-1));
						break;
					case 2:
						if(token.substr(0, 1).compare("[")) throw(2);
						yst = std::stod(token.substr(1));
						break;
					case 3:
						if(token.substr(token.size()-2).compare("]]")) throw(3);
						yend = std::stod(token.substr(0, token.size()-2));
						break;
					case 4:
						if(token.substr(0, 1).compare("[")) throw(4);
						nx = std::stoi(token.substr(1));
						break;
					case 5:
						if(token.substr(token.size()-1).compare("]")) throw(1);
						ny = std::stoi(token.substr(0, token.size()-1));
						break;
					default:
						throw(swcase);
						break;
				}
			}
			catch(std::invalid_argument const& ex){
				std::cout << "Corrupt schedule file, corrupted entry in line" << std::endl;
				std::cout << line << std::endl;
				throw("corrupt file");
			}
			catch(int exi){
			    std::cout << "Corrupt schedule file, corrupted line format in line" << std::endl;
				std::cout << line << std::endl;
				throw("corrupt file");
			}
			++swcase;
		}
		domains_start.push_back({xst, yst});
		domains_end.push_back({xend, yend});
		dim_nums.push_back({nx, ny});

		//Checking domain ordering
		if(dom_read == domains_read) ++domains_read;
		else{
			std::cout << "Corrupt schedule file, the domain " << dom_read << " is out of order" << std::endl;
			std::cout << line << std::endl;
			throw("corrupt file");		
		}
	}

	return Domain_Partitioning(domains_start, domains_end, dim_nums);
}

std::ostream &operator<<(std::ostream &os, const Domain_Partitioning &part){
	std::ios oldState(nullptr);
	oldState.copyfmt(os);

	os << std::setprecision(16);
	for(const Bulk_Domain *dom = part.first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
		os << "Domain " << dom->bulk_domain_number << ": [";
		for(int dim = 0, is_first = 1; dim < part.n_dim; ++dim){
			if(is_first) is_first = 0;
			else os << ", ";
			os << '[' << dom->domain_start(dim) << ", " << dom->domain_end(dim) << ']'; 
		}
		os << "], [";
		for(int dim = 0, is_first = 1; dim < part.n_dim; ++dim){
			if(is_first) is_first = 0;
			else os << ", ";
			os << dom->dim_num(dim);
		}
		os << "]" << std::endl;
	}

	os.copyfmt(oldState);
	return os;
}