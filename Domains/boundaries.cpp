#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <string>
#include <iomanip>

#include <domains.h>
#include <boundaries.h>

Bulk_Domain::Bulk_Domain(const std::vector<std::string> &instruct): Bulk_Domain(){
	std::vector<double> dom_start_load;
	std::vector<double> dom_end_load;
	std::vector<int> dim_n_load;

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
				if(piece.compare(1, 6, "points") == 0){
					std::string token = piece.substr(8, piece.length()-8);
					dim_n_load.push_back(std::stoi(token, &sztp));
				}
			}
		}
	}

	if(dom_start_load.size() != dom_end_load.size()){
		std::cout << "Domain reading error\n";
		throw("Domain reading error");
	}
	
	ndim = dom_start_load.size();
	dom_start = dom_start_load;
	dom_end = dom_end_load;

	dim_n.resize(ndim);
	for(int dm = 0; dm < ndim; dm++) dim_n[dm] = 0;
	if(dim_n_load.size() == 0) return;

	if(ndim != dim_n_load.size()){
		std::cout << "Incomplete number of points informaion\n";
		throw("Incomplete number of points informaion");
	}

	sz = 1;
	for(int dm = 0; dm < ndim; dm++){
		dim_n[dm] = dim_n_load[dm];
		sz *= dim_n[dm];
	}
}

Bulk_Domain::Bulk_Domain(const std::string &instruct): Bulk_Domain(){
	std::vector<double> dom_start_load;
	std::vector<double> dom_end_load;
	std::vector<int> dim_n_load;
	
	std::istringstream iss(instruct);
	int num_read = 0;
	std::string piece;
	bool run = true;
	
	while(run){
		run = (std::getline(iss, piece, ',').peek() != EOF);
		std::string::size_type sztp;
		switch((num_read++)%3){
			case 0:
				dim_n_load.push_back(std::stoi(piece, &sztp));
				break;
			case 1:
				dom_start_load.push_back(std::stod(piece, &sztp));
				break;
			case 2:
				dom_end_load.push_back(std::stod(piece, &sztp));
				break;
		}
	}

	if(dim_n_load.size() != dom_start_load.size() || dom_start_load.size() != dom_end_load.size()){
		std::cout << "Domain reading error\n";
		throw("Domain reading error");
	}

	sz = 1;
	ndim = dom_start_load.size();
	dim_n.resize(ndim);
	dom_start = dom_start_load;
	dom_end = dom_end_load;
	for(int dm = 0; dm < ndim; dm++){
		dim_n[dm] = dim_n_load[dm];
		sz *= dim_n[dm];
	}
}

void Domain_Partitioning::save_domains(const std::string &file_name){
	if(first_dom == NULL){
		std::cout << "The partitioning has not been initialized\n";
		throw("The partitioning has not been initialized");
	}

	std::ofstream out;
	out.open(file_name);

	int domain_num = 0;
	Bulk_Domain *domain_tmp = first_bulk_dom;

	out << "Number of dimensions: " << ndim << "\n" << std::endl;
	out << "Partition sizes: " << std::endl;
	for(int dm = 0; dm < ndim; dm++){
		out << "Dimension: " << dm << ", start: " << part_start[dm] << ", end: " << part_end[dm] << std::endl;
	}
	out << std::endl;

	while(domain_tmp != NULL){
		out << "Domain: " << domain_num++ << std::endl;
		out << *domain_tmp << std::endl;
		domain_tmp = domain_tmp->dom_bulk_next;
	}

	out.close();
	return;
}

Generic_Domain* Domain_Partitioning::load_single_domain(const std::vector<std::string> &domain_info, Generic_Domain *prev_dom){
	Bulk_Domain *prev_bulk_dom;
	if((prev_bulk_dom = dynamic_cast<Bulk_Domain*>(prev_dom)) == NULL && prev_dom != NULL){
		std::cout << "Incorrect previous domain type in single domain loader" << std::endl;
		throw("Incorrect previous domain type in single domain loader");
	}
	Bulk_Domain *domain_tmp = new Bulk_Domain(domain_info);
	if(prev_dom != NULL){
		domain_tmp->dom_prev = prev_bulk_dom;
		prev_bulk_dom->dom_next = domain_tmp;

		domain_tmp->dom_bulk_prev = prev_bulk_dom;
		prev_bulk_dom->dom_bulk_next = domain_tmp;
	}
	return domain_tmp;
}

void Domain_Partitioning::load_domains(const std::string &file_name){
	Generic_Domain_Partitioning::load_domains(file_name);
	first_bulk_dom = dynamic_cast<Bulk_Domain*>(first_dom);
	Bulk_Domain *bulk_dom_tmp = first_bulk_dom;
	int num = 0;
	while(bulk_dom_tmp != NULL){
		bulk_dom_tmp->b_dom_number = num++;
		bulk_dom_tmp = bulk_dom_tmp->dom_bulk_next;
	}
}

void Domain_Partitioning::load_domains(const std::vector<std::string> &instruct){
	if(first_dom != NULL){
		std::cout << "The partitioning has already been initialized" << std::endl;
		throw("The partitioning has already been initialized");
	}
	if(instruct.size() == 0){
		std::cout << "Instructions given are empty" << std::endl;
		throw("Instructions given are empty");
	}

	Bulk_Domain *domain_tmp;
	for(int dm = 0; dm < instruct.size(); dm++){
		if(first_dom == NULL){
			first_dom = first_bulk_dom = domain_tmp = new Bulk_Domain(instruct[dm]);
			ndim = domain_tmp->ndim;
			part_start = domain_tmp->dom_start;
			part_end = domain_tmp->dom_end;
		}
		else{
			domain_tmp->dom_next = domain_tmp->dom_bulk_next = new Bulk_Domain(instruct[dm]);
			(domain_tmp->dom_bulk_next)->dom_prev = (domain_tmp->dom_bulk_next)->dom_bulk_prev = domain_tmp;
			domain_tmp = domain_tmp->dom_bulk_next;
		}
		domain_tmp->b_dom_number = dm;
		for(int i = 0; i < ndim; i++){
			part_start[i] = (domain_tmp->dom_start[i] < part_start[i] ? domain_tmp->dom_start[i] : part_start[i]);
			part_end[i] = (domain_tmp->dom_end[i] > part_end[i] ? domain_tmp->dom_end[i] : part_end[i]);
		}
	}

	verify_ordering();
	verify_consistency();
	return;
}

void Domain_Partitioning::verify_ordering(bool verbose){
	Generic_Domain_Partitioning::verify_ordering(verbose);

	if(verbose) std::cout << "Verifying bulk domain ordering:\n\n";

	std::vector<Bulk_Domain*> dom_list = {first_bulk_dom};
	Bulk_Domain *tmp = first_bulk_dom;

	int count = 0;
	if(verbose) std::cout << "Counting domain " << count++ << std::endl;

	while(tmp->dom_bulk_next != NULL){
		dom_list.push_back(tmp->dom_bulk_next);
		if(verbose) std::cout << "Counting domain " << count++ << std::endl;
		tmp = tmp->dom_bulk_next;
	}

	if(verbose) std::cout << "\nCounting backwards:\n" << std::endl;

	for(int i = dom_list.size()-1; i >= 0; i--){
		if(dom_list[i] != tmp){
			std::cout << "Ordering is inconsistent" << std::endl;
			throw("Ordering is inconsistent");
		}

		if(verbose) std::cout << "Counting domain " << --count << std::endl;
		tmp = tmp->dom_bulk_prev;
	}
	if(verbose) std::cout << std::endl;
}

void Domain_Partitioning::verify_boundary_point_splitting(const std::unordered_set<int> &exclude_left, const std::unordered_set<int> &exclude_right){
	if(!boundary_domains_inserted){
		std::cout << "Warning! Trying to verify boundary point splitting for a partition with no boundaries" << std::endl;
		return;
	}
	const Boundary_Domain *boundary;
	std::vector<int> area_begin(ndim-1), area_end(ndim-1);

	for(const Bulk_Domain *bulk = first_bulk_dom; bulk != NULL; bulk = bulk->dom_bulk_next){
		for(int dm = 0; dm < ndim; dm++){
			for(int i = 0, indx = 0; i < ndim; i++){
				if(i == dm) continue;

				if(i < dm && !(exclude_left.count(i) && std::fabs(part_start[i] - bulk->dom_start[i]) < COMP_TOL))
					area_begin[indx] = 1;
				else
					area_begin[indx] = 0;
				if(i < dm && !(exclude_right.count(i) && std::fabs(part_end[i] - bulk->dom_end[i]) < COMP_TOL))
					area_end[indx++] = bulk->dim_n[i]-1;
				else
					area_end[indx++] = bulk->dim_n[i];
			}

			std::vector<std::vector<int>> begins_left;
			std::vector<std::vector<int>> ends_left;
			std::vector<std::vector<int>> begins_right;
			std::vector<std::vector<int>> ends_right;

			for(const Generic_Domain *tmp_dom = first_dom; tmp_dom != NULL; tmp_dom = tmp_dom->next_domain()){
				boundary = dynamic_cast<const Boundary_Domain *>(tmp_dom);
				if(boundary == NULL || boundary->bound_dim != dm) continue;
				if(boundary->dom_left_side == bulk){
					std::vector<int> end(boundary->left_n_end);
					for(int dm = 0; dm < ndim-1; dm++){
						if(end[dm] != area_end[dm]) end[dm]++;
					}
					begins_right.push_back(boundary->left_n_st);
					ends_right.push_back(end);
				}
				if(boundary->dom_right_side == bulk){
					std::vector<int> end(boundary->right_n_end);
					for(int dm = 0; dm < ndim-1; dm++){
						if(end[dm] != area_end[dm]) end[dm]++;
					}
					begins_left.push_back(boundary->right_n_st);
					ends_left.push_back(end);
				}
			}
			try{
				check_consistency<int>(begins_left, ends_left, area_begin, area_end);
			}
			catch(const char *err){
				std::cout << "Begins left:" << std::endl;
				for(auto el : begins_left){
					for(int i = 0; i < el.size(); ++i){
						if(i) std::cout << ',';
						std::cout << el[i];
					}
					std::cout << std::endl;
				}
				std::cout << "Ends left:" << std::endl;
				for(auto el : ends_left){
					for(int i = 0; i < el.size(); ++i){
						if(i) std::cout << ',';
						std::cout << el[i];
					}
					std::cout << std::endl;
				}
				std::cout << "Problem detected while verifying left boundary along dimension " << dm << " of the domain " << bulk->bulk_domain_number << std::endl;
				throw(err);
			}
			try{
				check_consistency<int>(begins_right, ends_right, area_begin, area_end);
			}
			catch(const char *err){
				std::cout << "Begins right:" << std::endl;
				for(auto el : begins_right){
					for(int i = 0; i < el.size(); ++i){
						if(i) std::cout << ',';
						std::cout << el[i];
					}
					std::cout << std::endl;
				}
				std::cout << "Ends right:" << std::endl;
				for(auto el : ends_right){
					for(int i = 0; i < el.size(); ++i){
						if(i) std::cout << ',';
						std::cout << el[i];
					}
					std::cout << std::endl;
				}
				std::cout << "Problem detected while verifying right boundary along dimension " << dm << " of the domain " << bulk->bulk_domain_number << std::endl;
				throw(err);
			}
		}
	}
}

void print_partitioning(const Generic_Domain_Partitioning &part){
	int i = 0;
	const Bulk_Domain *bulk_dom_tmp;
	const Boundary_Domain *bound_dom_tmp;
	for(const Generic_Domain *dom = part.first_domain; dom != NULL; dom = dom->next_domain()){
		std::cout << "Domain: " << i << std::endl;
		if((bulk_dom_tmp = dynamic_cast<const Bulk_Domain*>(dom)) != NULL){
			std::cout << "Domain type: Bulk Domain " << bulk_dom_tmp->bulk_domain_number << std::endl;
		}
		if((bound_dom_tmp = dynamic_cast<const Boundary_Domain*>(dom)) != NULL){
			std::cout << "Domain type: Boundary Domain" << std::endl;
		}
		std::cout << *dom << std::endl;
		i++;
	}
}