#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <iomanip>
#include <filesystem>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <domains.h>
#include <boundaries.h>

std::ostream& operator<<(std::ostream& os, const Generic_Domain &dom){
	dom.print_domain(os);
	return os;
}

bool check_adjacency(const Generic_Domain *dom1, const Generic_Domain *dom2){
	int bound_dim = dom1->ndim;

	if(dom1->ndim != dom2->ndim){
		std::cout << "Domains dimensions don't match in adjacency checker" << std::endl;
		throw("Domains dimensions don't match in adjacency checker");
	}

	for(int dm = 0; dm < dom1->ndim; dm++){
		if(dom1->dom_end[dm] == dom2->dom_start[dm] || dom1->dom_start[dm] == dom2->dom_end[dm]){
			bool candidate = true;
			for(int dm2 = 0; dm2 < dom1->ndim; dm2++){
				if(dm2 == dm) continue;
				if(MAX(dom1->dom_start[dm2], dom2->dom_start[dm2]) > MIN(dom1->dom_end[dm2], dom2->dom_end[dm2])){
					candidate = false;
					break;
				}
			}
			if(candidate){
				if(bound_dim == dom1->ndim) bound_dim = dm;
				else{
					std::cout << "Warning! At least two possible boundary dimensions in adjacency checker" << std::endl;
					throw("At least two possible boundary dimensions in adjacency checker");
					break;
				}
			}
		}
	}

	if(bound_dim == dom1->ndim){
		return false;
	}
	else return true;
}

void drop_zero_volume_domains(std::vector<Generic_Domain*> &list){
	std::vector<Generic_Domain*> new_list;
	for(int i = 0; i < list.size(); i++){
		if(!list[i]->is_zero_volume()) new_list.push_back(list[i]);
	}
	list = new_list;
	return;
}

void Generic_Domain_Partitioning::verify_ordering(bool verbose) const{
	std::vector<Generic_Domain*> dom_list = {first_dom};
	Generic_Domain *tmp = first_dom;

	if(verbose) std::cout << "Verifying domain ordering:\n\n";

	int count = 0;
	if(verbose) std::cout << "Counting domain " << count++ << std::endl;

	while(tmp->dom_next != NULL){
		dom_list.push_back(tmp->dom_next);
		if(verbose) std::cout << "Counting domain " << count++ << std::endl;
		tmp = tmp->dom_next;
	}

	if(verbose) std::cout << "\nCounting backwards:\n" << std::endl;

	for(int i = dom_list.size()-1; i >= 0; i--){
		if(dom_list[i] != tmp){
			std::cout << "Ordering is inconsistent" << std::endl;
			throw("Ordering is inconsistent");
		}

		if(verbose) std::cout << "Counting domain " << --count << std::endl;
		tmp = tmp->dom_prev;
	}

	if(verbose) std::cout << std::endl;
}

void Generic_Domain_Partitioning::save_domains(const std::string &file_name) const{
	if(first_dom == NULL){
		std::cout << "The partitioning has not been initialized\n";
		throw("The partitioning has not been initialized");
	}

	std::ofstream out;
	out.open(file_name);

	int domain_num = 0;
	Generic_Domain *domain_tmp = first_dom;

	out << "Number of dimensions: " << ndim << "\n" << std::endl;
	out << "Partition sizes: " << std::endl;
	for(int dm = 0; dm < ndim; dm++){
		out << "Dimension: " << dm << ", start: " << part_start[dm] << ", end: " << part_end[dm] << std::endl;
	}
	out << std::endl;

	while(domain_tmp != NULL){
		out << "Domain: " << domain_num++ << std::endl;
		out << *domain_tmp << std::endl;
		domain_tmp = domain_tmp->dom_next;
	}

	out.close();
	return;
}

Generic_Domain::Generic_Domain(const std::vector<std::string> &instruct): n_dim(ndim){
	std::vector<double> dom_start_load;
	std::vector<double> dom_end_load;

	int nlines = instruct.size();
	for(int ln = 0; ln < nlines; ln++){
		if(instruct[ln].compare(0, 9, "Dimension") == 0){
			std::istringstream iss(instruct[ln]);
			std::string piece;
			bool run = true;
			while(run){
				run = (std::getline(iss, piece, ',').peek() != EOF);
				std::string::size_type sztp;
				if(piece.compare(1, 5, "start") == 0){
					std::string token = piece.substr(7, piece.length()-7);
					dom_start_load.push_back(std::stod(token, &sztp));
				}
				if(piece.compare(1, 3, "end") == 0){
					std::string token = piece.substr(5, piece.length()-5);
					dom_end_load.push_back(std::stod(token, &sztp));
				}
			}
		}
	}

	if(dom_start_load.size() != dom_end_load.size()){
		std::cout << "Domain reading error\n";
		throw("Domain reading error");
	}

	ndim = dom_start_load.size();
	dom_prev = NULL;
	dom_next = NULL;
	
	dom_start = dom_start_load;
	dom_end = dom_end_load;
}

Generic_Domain* Generic_Domain_Partitioning::load_single_domain(const std::vector<std::string> &domain_info, Generic_Domain *prev_dom){
	Generic_Domain *domain_tmp = new Generic_Domain(domain_info);
	if(prev_dom != NULL){
		domain_tmp->dom_prev = prev_dom;
		prev_dom->dom_next = domain_tmp;
	}
	return domain_tmp;
}

void Generic_Domain_Partitioning::load_domains(const std::string &file_name){
	if(first_dom != NULL){
		std::cout << "The partitioning has already been initialized\n";
		throw("The partitioning has already been initialized");
	}

	std::ifstream file(file_name);
	if(file.is_open() == false){
		std::cout << "No domains.dat file found\n";
		throw("No domains.dat file found");
	}

	Generic_Domain *domain_tmp;
	bool record_domain_info = false;
	std::vector<std::string> domain_info;

	std::string line;
	bool last_domain = false;
	while(1){
		if(!std::getline(file, line)) last_domain = true;
		if(line.compare(0, 7, "Domain:") == 0 || last_domain){
			if(record_domain_info){
				if(first_dom == NULL){
					first_dom = domain_tmp = load_single_domain(domain_info, NULL);
					ndim = first_dom->ndim;
					
					part_start = first_dom->dom_start;
					part_end = first_dom->dom_end;
				}
				else{
					domain_tmp = load_single_domain(domain_info, domain_tmp);

					for(int i = 0; i < ndim; i++){
						part_start[i] = (domain_tmp->dom_start[i] < part_start[i] ? domain_tmp->dom_start[i] : part_start[i]);
						part_end[i] = (domain_tmp->dom_end[i] > part_end[i] ? domain_tmp->dom_end[i] : part_end[i]);
					}
					
					if(ndim != domain_tmp->ndim){
						std::cout << "Inconsistent number of dimensions in domain file\n";
						throw("Inconsistent number of dimensions in domain file");
					}
				}
			}
			domain_info.clear();
			record_domain_info = true;
			if(last_domain) break;
		}
		else if(record_domain_info)
			domain_info.push_back(line);
	}

	file.close();
	verify_consistency();
	return;
}

void Generic_Domain_Partitioning::list_domains(int dm, std::vector<double> &posit, std::vector<Generic_Domain*> &layer) const{
	layer.resize(0);
	Generic_Domain *domain_tmp = first_dom;

	while(domain_tmp != NULL){
		bool in_layer = true;
		for(int i = 0; i < ndim; i++){
			if(i == dm) continue;
			if(posit[i] < domain_tmp->dom_start[i] || posit[i] >= domain_tmp->dom_end[i]){
				in_layer = false;
				break;
			}
		}
		if(in_layer) layer.push_back(domain_tmp);
		domain_tmp = domain_tmp->dom_next;
	}

	if(layer.size() == 0){
		std::cout << "Domain partitioning error, no domains found in the layer\n\n";
		throw("Domain partitioning error, no domains found in the layer");
	}

	struct comp_dom_str{
		int dim;
		comp_dom_str(int dimi): dim(dimi) {};
	  	bool operator()(const Generic_Domain *domain_i, const Generic_Domain *domain_j) {
	  		if(domain_i->dom_start[dim] < domain_j->dom_start[dim]) return true;
	  		if(domain_i->dom_start[dim] == domain_j->dom_start[dim] && domain_i->dom_start[dim] == domain_i->dom_end[dim]) return true;
	  		return false;
	  	}
	} compare_function(dm);

	std::sort(layer.begin(), layer.end(), compare_function);

	//Checking consistency
	int ndomains = layer.size();
	if(part_start[dm] != layer[0]->dom_start[dm] || part_end[dm] != layer[ndomains-1]->dom_end[dm]){
		std::cout << "Domain partitioning error, domains don't cover the whole partition\n";
		for(int j = 0; j < ndomains; j++){
			std::cout << "Domain:\n";
			std::cout << *layer[j];
		}
		throw("Domain partitioning error, domains don't cover the whole partition");
	}

	for(int i = 0; i < ndomains-1; i++){
		if(layer[i]->dom_end[dm] != layer[i+1]->dom_start[dm]){
			std::cout << "Domain partitioning error, domains clash\n";
			for(int j = 0; j < ndomains; j++){
				std::cout << "Domain:\n";
				std::cout << *layer[j];
			}
			throw("Domain partitioning error, domains clash");
		}
	}
}

void Generic_Domain_Partitioning::verify_consistency(void) const{
	std::vector<double> current_coord(part_start);
	std::vector<double> next_coord(part_end);
	std::vector<Generic_Domain*> layer;

	int dim_change = 0;

	while(1){
		list_domains(ndim-1, current_coord, layer);
		int nlayer = layer.size();
		for(int dm = 0; dm < ndim-1; dm++)
			for(int i = 0; i < nlayer; i++) 
				if(next_coord[dm] > layer[i]->dom_end[dm] && layer[i]->dom_end[dm] != layer[i]->dom_start[dm])
					next_coord[dm] = layer[i]->dom_end[dm];

		for(dim_change = ndim-2; dim_change >= 0; dim_change--){
			if(next_coord[dim_change] != part_end[dim_change]) break;
		}
		if(dim_change < 0) break;
		for(int dm = ndim-2; dm > dim_change; dm--)
			current_coord[dm] = part_start[dm];
		current_coord[dim_change] = next_coord[dim_change];
		for(int dm = ndim-2; dm >= dim_change; dm--)
			next_coord[dm] = part_end[dm];
	}

	// std::vector<std::vector<double>> begins, ends;
	// for(Generic_Domain *dom = first_dom; dom != NULL; dom = dom->dom_next){
	// 	begins.push_back(dom->dom_start);
	// 	ends.push_back(dom->dom_end);
	// }
	// check_consistency<double>(begins, ends, part_start, part_end);
}

void Generic_Domain_Partitioning::print_domain_2D(const std::string &file_name, int size_x, int size_y, int dim_i, int dim_j, double *posit_st) const{
	std::ofstream out;
	out.open(file_name);

	int icurr = 0, inext;
	std::vector<Generic_Domain*> layer;
	std::vector<double> posit(ndim);
	for(int i = 0; i < ndim; i++){
		if(i == dim_i || i == dim_j) continue;
		posit[i] = posit_st[icurr++];
	}
	posit[dim_i] = part_start[dim_i];
	posit[dim_j] = part_start[dim_j];
	icurr = 0;

	double coord_next;

	out << "P1\n" << size_x << " " << size_y << std::endl;

	while(posit[dim_i] < part_end[dim_i]){
		list_domains(dim_j, posit, layer);
		drop_zero_volume_domains(layer);
		int nlayer = layer.size();

		std::vector<int> jbound(nlayer+1);
		int ist, iend;

		//Getting y boundaries
		for(int dmi = 0; dmi < nlayer; dmi++){
			jbound[dmi] = (int) size_y*(layer[dmi]->dom_start[dim_j] - part_start[dim_j])/(part_end[dim_j] - part_start[dim_j]);
		}
		jbound[nlayer] = size_y-1;

		//Removing duplicates
		jbound.erase(std::unique(jbound.begin(),jbound.end()), jbound.end());
		nlayer = jbound.size()-1;
		
		//Drawing y boundaries
		for(int dmi = 0; dmi < nlayer; dmi++){
			out << "1";
			if(layer[dmi]->dom_start[dim_i] == posit[dim_i]){
				for(int j = jbound[dmi]+1; j < jbound[dmi+1]; j++) out << "1";
			}
			else{
				for(int j = jbound[dmi]+1; j < jbound[dmi+1]; j++) out << "0";	
			}
		}
		out << "1\n";

		//Getting the x coordinate of next y boundaries
		coord_next = layer[0]->dom_end[dim_i];
		for(int i = 1; i < layer.size(); i++)
			if(layer[i]->dom_end[dim_i] < coord_next && layer[i]->dom_start[dim_i] != layer[i]->dom_end[dim_i])
				coord_next = layer[i]->dom_end[dim_i]; 
		inext = (int) size_x*(coord_next - part_start[dim_i])/(part_end[dim_i] - part_start[dim_i]);

		//Drawing y boundaries
		for(int i = icurr; i < inext; i++){
			for(int dmi = 0; dmi < nlayer; dmi++){
				out << "1";
				for(int j = jbound[dmi]+1; j < jbound[dmi+1]; j++) out << "0";	
			}
			out << "1\n";
		}

		posit[dim_i] = coord_next;
		icurr = inext;
	}

	for(int i = 0; i < size_y-1; i++) out << "1";
	out << "1" << std::endl;

	out.close();
	return;
}

void Generic_Domain_Partitioning::print_domain_3D(const std::string &folder_name, int size_x, int size_y, int dim_i, int dim_j, int dim_k) const{
	if(ndim != 3){
		std::cout << "Wrong dimensionality\n";
		throw("Wrong dimensionality");
	}

	//Creating output folder
	struct stat sttus = {0};
	if (stat(folder_name.c_str(), &sttus) == -1) {
		mkdir(folder_name.c_str(), 0700);
	}
	std::ofstream out;

	int icurr = 0, inext, nplot = 0;
	std::vector<Generic_Domain*> layer;
	std::vector<double> posit(part_start);

	std::string file_name;
	double coord_i_next;

	while(posit[dim_i] < part_end[dim_i]){
		double coord_j_next;

		file_name = folder_name +  "/layer_" + std::to_string(nplot++) +  ".pbm";
		out.open(file_name);
		out << "P1\n" << size_x << " " << size_y << std::endl;

		coord_i_next = part_end[dim_i];

		while(posit[dim_j] < part_end[dim_j]){
			list_domains(dim_k, posit, layer);
			drop_zero_volume_domains(layer);
			int nlayer = layer.size();

			std::vector<int> kbound(nlayer+1);
			int ist, iend;

			//Getting k boundaries
			for(int dmi = 0; dmi < nlayer; dmi++){
				kbound[dmi] = (int) size_y*(layer[dmi]->dom_start[dim_k] - part_start[dim_k])/(part_end[dim_k] - part_start[dim_k]);
			}
			kbound[nlayer] = size_y-1;

			//Removing duplicates
			kbound.erase(std::unique(kbound.begin(),kbound.end()), kbound.end());
			nlayer = kbound.size()-1;
			
			//Drawing k boundaries
			for(int dmi = 0; dmi < nlayer; dmi++){
				out << "1";
				if(layer[dmi]->dom_start[dim_j] == posit[dim_j]){
					for(int j = kbound[dmi]+1; j < kbound[dmi+1]; j++) out << "1";
				}
				else{
					for(int j = kbound[dmi]+1; j < kbound[dmi+1]; j++) out << "0";	
				}
			}
			out << "1\n";

			for(int i = 0; i < layer.size(); i++)
				if(layer[i]->dom_end[dim_i] < coord_i_next && layer[i]->dom_start[dim_i] != layer[i]->dom_end[dim_i])
					coord_i_next = layer[i]->dom_end[dim_i]; 

			//Getting the j coordinate of next k boundaries
			coord_j_next = layer[0]->dom_end[dim_j];
			for(int i = 1; i < layer.size(); i++)
				if(layer[i]->dom_end[dim_j] < coord_j_next && layer[i]->dom_start[dim_j] != layer[i]->dom_end[dim_j])
					coord_j_next = layer[i]->dom_end[dim_j]; 
			inext = (int) size_x*(coord_j_next - part_start[dim_j])/(part_end[dim_j] - part_start[dim_j]);

			//Drawing k boundaries
			for(int i = icurr; i < inext; i++){
				for(int dmi = 0; dmi < nlayer; dmi++){
					out << "1";
					for(int j = kbound[dmi]+1; j < kbound[dmi+1]; j++) out << "0";	
				}
				out << "1\n";
			}

			posit[dim_j] = coord_j_next;
			icurr = inext;
		}

		for(int i = 0; i < size_y-1; i++) out << "1";
		out << "1" << std::endl;

		out.close();

		posit[dim_i] = coord_i_next;
		posit[dim_j] = part_start[dim_j];
		icurr = 0;
	}

	return;
}

void Generic_Domain_Partitioning::move_domain_wall(const std::vector<double> &wall_posit, double new_posit){
	int wall_dim_left = -1, wall_dim_right = -1;
	std::vector<Generic_Domain *> domains_left, domains_right;
	for(Generic_Domain *dom = first_dom; dom != NULL; dom = dom->dom_next){
		bool is_inside = true;
		for(int i = 0; i < ndim; i++){
			if(dom->dom_start[i] - wall_posit[i] > COMP_TOL || wall_posit[i] - dom->dom_end[i] > COMP_TOL){
				is_inside = false;
				break;
			}
		}
		if(is_inside){
			for(int i = 0; i < ndim; i++){
				if(std::fabs(dom->dom_start[i] - wall_posit[i]) < COMP_TOL){
					if(wall_dim_right >= 0){
						std::cout << "Domain wall mover hit a corner" << std::endl;
						throw("Domain wall move fail");
					}
					wall_dim_right = i;
					if(domains_right.size() == 0){
						domains_right.push_back(dom);
					}
					else{
						std::cout << "Partition is broken, the wall cannot be moved" << std::endl;
						throw("Domain wall move fail");
					}
				}
				if(std::fabs(wall_posit[i] - dom->dom_end[i]) < COMP_TOL){
					if(wall_dim_left >= 0){
						std::cout << "Domain wall mover hit a corner" << std::endl;
						throw("Domain wall move fail");
					}
					wall_dim_left = i;
					if(domains_left.size() == 0){
						domains_left.push_back(dom);
					}
					else{
						std::cout << "Partition is broken, the wall cannot be moved" << std::endl;
						throw("Domain wall move fail");
					}
				}
			}
		}
	}
	if(domains_left.size() == 0 || domains_right.size() == 0){
		std::cout << "The outer boundary cannot be moved" << std::endl;
		throw("Domain wall move fail");
	}
	if(wall_dim_left < 0 || wall_dim_right < 0 || wall_dim_left != wall_dim_right){
		std::cout << "Partition is broken, the wall cannot be moved" << std::endl;
		throw("Domain wall move fail");
	}

	while(1){
		bool is_added = false;
		for(Generic_Domain *dom = first_dom; dom != NULL; dom = dom->dom_next){
			if(std::fabs(dom->dom_end[wall_dim_left] - wall_posit[wall_dim_left]) < COMP_TOL){
				bool is_adjacent = false, is_old = false;
				for(int i = 0; i < domains_left.size(); i++){
					if(dom == domains_left[i]){
						is_old = true;
						break;
					}
				}
				if(is_old) continue;
				for(int i = 0; i < domains_right.size(); i++){
					if(check_adjacency(dom, domains_right[i])){
						is_adjacent = true;
						break;
					}
				}
				if(is_adjacent)
					domains_left.push_back(dom);
			}
		}

		for(Generic_Domain *dom = first_dom; dom != NULL; dom = dom->dom_next){
			if(std::fabs(dom->dom_start[wall_dim_left] - wall_posit[wall_dim_left]) < COMP_TOL){
				bool is_adjacent = false, is_old = false;
				for(int i = 0; i < domains_right.size(); i++){
					if(dom == domains_right[i]){
						is_old = true;
						break;
					}
				}
				if(is_old) continue;
				for(int i = 0; i < domains_left.size(); i++){
					if(check_adjacency(dom, domains_left[i])){
						is_adjacent = true;
						break;
					}
				}
				if(is_adjacent){
					domains_right.push_back(dom);
					is_added = true;
				}
			}
		}

		if(!is_added) break;
	}

	for(int i = 0; i < domains_left.size(); i++){
		if(new_posit - (domains_left[i])->dom_start[wall_dim_left] <= 0.0){
			std::cout << "Attempted domain move destroys one of the domains" << std::endl;
			throw("Domain wall move fail");
		}
	}
	for(int i = 0; i < domains_right.size(); i++){
		if(new_posit - (domains_right[i])->dom_end[wall_dim_left] >= 0.0){
			std::cout << "Attempted domain move destroys one of the domains" << std::endl;
			throw("Domain wall move fail");
		}
	}

	for(int i = 0; i < domains_left.size(); i++)
		(domains_left[i])->dom_end[wall_dim_left] = new_posit;
	for(int i = 0; i < domains_right.size(); i++)
		(domains_right[i])->dom_start[wall_dim_left] = new_posit;

	Generic_Domain_Partitioning::verify_ordering();
	verify_consistency();
}