#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <map>
#include <boundaries.h>
#include <extrapolation.h>

class Tracing_Schedule: public Extrapolation_Condition{
public:
	std::map<double, Domain_Partitioning>::iterator curr;
	std::map<double, Domain_Partitioning> schedule;

	Tracing_Schedule(): curr(schedule.end()) {};

	void load_schedule(const std::string &file_name){
		double curr_pslength, tmp;
		std::vector<std::vector<double>> domains_start;
		std::vector<std::vector<double>> domains_end;
		std::vector<std::vector<int>> dim_nums;
		int domains_read = 0;

		std::ifstream file(file_name);
		if(!file){
			throw("file not found");
		}
		std::string line;
		std::string token;
		std::string::size_type sztp;

		while(std::getline(file, line)){
			if(line.size() == 0)
				continue;
			std::istringstream iss(line);
			std::getline(iss, token, ':');

			//Checking the first half of the entry
			if(token.compare("Pseudolength") == 0){
				if(domains_start.size() != 0){
					try{
						schedule.insert({curr_pslength, Domain_Partitioning(domains_start, domains_end, dim_nums)});
					}
					catch(const char *ex){
						throw("corrupt file");
					}
					domains_read = 0; domains_start.clear(); domains_end.clear(); dim_nums.clear();
				}
				std::getline(iss, token, ':');
				try{
					curr_pslength = std::stod(token, &sztp);
				}
				catch(std::invalid_argument const& ex){
					std::cout << "Corrupt pseudolength value" << std::endl;
					std::cout << line << std::endl;
					throw("corrupt file");
				}
				continue;
			}
			else if(token.substr(0, 6).compare("Domain") != 0){
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
		if(domains_start.size() != 0)
			schedule.insert({curr_pslength, Domain_Partitioning(domains_start, domains_end, dim_nums)});
		else{
			std::cout << "Missing domain data for the last pseudolength value" << std::endl;
			throw("corrupt file");
		}
	}

	void setup_scheduling(double ps_length){
		curr = schedule.begin();
		while(curr != schedule.end()){
			if(curr->first > ps_length)
				break;
			++curr;
		}
	}

	void check_condition(Extrapolation &extra){
		if(curr == schedule.end()) return;
		if(extra.sol_current->ps_length > curr->first){
			extra.regrid_solutions(curr->second);
			++curr;
		}
	};

	void final_action(Extrapolation &extra, bool normal_exit){return;}
	friend std::ostream &operator<<(std::ostream &os, const Tracing_Schedule &sch);
};

std::ostream &operator<<(std::ostream &os, const Tracing_Schedule &sch){
	bool first;
	for(auto el = sch.schedule.begin(); el != sch.schedule.end(); ++el){
		if(el != sch.schedule.begin()) os << '\n';
		os << "Pseudolength: " << el->first << std::endl;
		os << "Domains:\n";

		for(const Bulk_Domain *dom = (el->second).first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			if(dom->bulk_domain_number) os << std::endl;
			os << "Domain: " << dom->bulk_domain_number << std::endl;
			os << *dom;
		}
	}
	return os;
}