#include <cstdio>
#include <vector>
#include <iostream>
#include <cmath>

#ifndef DOMAINS_H
#define DOMAINS_H

#ifndef MIN
#define MIN(a, b) (a < b ? a : b)
#endif
#ifndef MAX
#define MAX(a, b) (a > b ? a : b)
#endif

#define COMP_TOL 1.0e-15

class Generic_Domain{
private:
	int ndim;
	// Domains are ordered through pointers to previous and next domains 
	// Initial domain has previous domain pointer NULL 
	// Final domain has next domain pointer NULL
	Generic_Domain *dom_prev;
	Generic_Domain *dom_next;

	std::vector<double> dom_start;
	std::vector<double> dom_end;

	// Prints domain starting and ending coordinates along every dimension
	virtual void print_domain(std::ostream& os) const{
		for(int dm = 0; dm < ndim; dm++){
			os << "Dimension: " << dm << ", start: " << dom_start[dm] << ", end: " << dom_end[dm] << std::endl;
		}
	}

	// Splits domains into two at some postition after the start
	virtual void domain_split(int split_dim, double split_coord){
		if(split_dim >= ndim){
			std::cout << "Invalid domain splitting dimension\n";
			throw("Invalid domain splitting dimension");
		}
		if(split_coord < dom_start[split_dim] || split_coord > dom_end[split_dim]){
			std::cout << "Invalid domain splitting position\n";
			throw("Invalid domain splitting");
		}

		Generic_Domain *new_dom = new Generic_Domain(*this);

		new_dom->dom_next = dom_next;
		if(dom_next != NULL) dom_next->dom_prev = new_dom;

		new_dom->dom_prev = this;
		dom_next = new_dom;

		new_dom->dom_start[split_dim] = split_coord;
		dom_end[split_dim] = new_dom->dom_start[split_dim];
	}

	// Default constructor
	Generic_Domain(): n_dim(ndim), ndim(0), dom_prev(NULL), dom_next(NULL) {};

	// Basic domain constructor
	Generic_Domain(int ndimi): n_dim(ndim), ndim(ndimi), dom_start(ndimi), dom_end(ndimi), dom_prev(NULL), dom_next(NULL) {};

	// Copy constructor
	Generic_Domain(const Generic_Domain &g_dom): Generic_Domain(g_dom.ndim){
		dom_start = g_dom.dom_start;
		dom_end = g_dom.dom_end;
	};

	// Standard domain constructor using starting position and lengths
	Generic_Domain(const std::vector<double> &posit, const std::vector<double> &lengths): Generic_Domain(posit.size()){
		if(posit.size() != lengths.size()){
			std::cout << "Dimensions of arrays are inconsistent in Generic Domain constructor" << std::endl;
			throw("Dimensions of arrays inconsistent");
		}

		for(int dm = 0; dm < ndim; dm++){
			dom_start[dm] = posit[dm];
			dom_end[dm] = posit[dm]+lengths[dm];
		}
	}

	virtual ~Generic_Domain() {};

	// Constructs domain from domain file description
	Generic_Domain(const std::vector<std::string> &instruct);

public:
	const int &n_dim;
	// Returns the number of dimensions
	double domain_start(int dm) const {return dom_start[dm];};
	double domain_end(int dm) const {return dom_end[dm];};
	const std::vector<double> &domain_start() const {return dom_start;};
	const std::vector<double> &domain_end() const {return dom_end;};

	// Checks if the domain has zero volume
	bool is_zero_volume() const{
		bool is_zero_volume = false;
		for(int dm = 0; dm < ndim; dm++){
			if(dom_start[dm] == dom_end[dm]) is_zero_volume = true;
		}
		return is_zero_volume;
	}

	// Comparison operator
	virtual bool operator==(const Generic_Domain &dom) const{
		if(ndim != dom.ndim) return false;
		for(int dm = 0; dm < ndim; dm++){
			if(dom_start[dm] != dom.dom_start[dm] || dom_end[dm] != dom.dom_end[dm]) return false;
		}
		return true;
	}

	// Checks if domains are adjacent, returns false if the domains are not adjacent
	// Throws exception if the domains adjacency region has dim-2 or less dimensions 
	friend bool check_adjacency(const Generic_Domain *dom1, const Generic_Domain *dom2);

	// Print operator invokes a virtual function specified separately for base class and dependants
	friend std::ostream& operator<<(std::ostream& os, const Generic_Domain &dom);

	// Check if the domain is initial or final
	bool is_initial() const{
		if(dom_prev == NULL) return true;
		else return false;
	}
	bool is_final() const{
		if(dom_next == NULL) return true;
		else return false;
	}
	
	// Outputs a pointer to previous and next domains
	Generic_Domain *previous_domain(){return dom_prev;}
	const Generic_Domain *previous_domain() const {return dom_prev;}
	Generic_Domain *next_domain(){return dom_next;}
	const Generic_Domain *next_domain() const {return dom_next;}

	friend class Bulk_Domain;
	friend class Boundary_Domain;
	friend class Generic_Domain_Partitioning;
	friend class Domain_Partitioning;
};


class Generic_Domain_Partitioning{
private:
	int ndim; //Number of dimensions
	Generic_Domain *first_dom;
	std::vector<double> part_start; //Leftmost coordinates of the leftmost domain of the partitioning
	std::vector<double> part_end; //Rightmost coordinates of the rightmost domain of the partitioning

	Generic_Domain_Partitioning(): n_dim(ndim), first_domain(first_dom), ndim(0), first_dom(NULL) {};

	//Lists all domains along the dimension dm, when other coordinates are fixed to be posit
	//Also checks if the domains cover all the range from part_start[dm] to part_end[dm] and throws exception if not
	void list_domains(int dm, std::vector<double> &posit, std::vector<Generic_Domain*> &layer) const;

	virtual Generic_Domain* load_single_domain(const std::vector<std::string> &domain_info, Generic_Domain *prev_dom);
	virtual void load_domains(const std::string &file_name);

	friend class Domain_Partitioning;

public:
	// Returns the number of dimensions
	const int &n_dim;

	double partition_start(int dm) const {return part_start[dm];};
	double partition_end(int dm) const {return part_end[dm];};
	const std::vector<double> &partition_start() const {return part_start;}
	const std::vector<double> &partition_end() const {return part_end;}

	//Returns first domain pointer
	Generic_Domain* const &first_domain;

	Generic_Domain_Partitioning &operator=(const Generic_Domain_Partitioning &part){
		ndim = part.ndim;
		part_start = part.part_start;
		part_end = part.part_end;

		Generic_Domain *dom = first_dom, *dom_tmp;
		first_dom = NULL;
		while(dom != NULL){
			dom_tmp = dom;
			dom = dom->dom_next;
			delete dom_tmp;
		}
		for(const Generic_Domain *dom_cpy = part.first_dom; dom_cpy != NULL; dom_cpy = dom_cpy->dom_next){
			dom_tmp = new Generic_Domain(*dom_cpy);
			if(first_dom == NULL) first_dom = dom_tmp;
			dom_tmp->dom_prev = dom;
			if(dom != NULL) dom->dom_next = dom_tmp;
			dom = dom_tmp;
		}

		verify_ordering();
		return *this;
	}

	Generic_Domain_Partitioning &operator=(Generic_Domain_Partitioning &&part){
		ndim = part.ndim;
		part_start = part.part_start;
		part_end = part.part_end;

		first_dom = part.first_dom;
		part.first_dom = NULL;

		return *this;
	}

	// Copy constructor
	Generic_Domain_Partitioning(const Generic_Domain_Partitioning &part): n_dim(ndim), first_domain(first_dom), 
		ndim(part.ndim), part_start(part.part_start), part_end(part.part_end){
		const Generic_Domain* part_old_tmp = part.first_dom;
		Generic_Domain* part_new_tmp;

		first_dom = part_new_tmp = new Generic_Domain(*part_old_tmp);
		first_dom->dom_prev = NULL;
		part_old_tmp = part_old_tmp->dom_next;

		while(part_old_tmp != NULL){
			part_new_tmp->dom_next = new Generic_Domain(*part_old_tmp);
			(part_new_tmp->dom_next)->dom_prev = part_new_tmp;
			part_new_tmp = part_new_tmp->dom_next;

			part_old_tmp = part_old_tmp->dom_next;
		}
		part_new_tmp->dom_next = NULL;
		verify_ordering();
	}

	// Move constructor
	Generic_Domain_Partitioning(Generic_Domain_Partitioning &&part): n_dim(ndim), first_domain(first_dom), 
		ndim(part.ndim), part_start(part.part_start), part_end(part.part_end), first_dom(part.first_dom)
	{
		part.first_dom = NULL;
	}

	// Initialize with a single domain constructor 
	Generic_Domain_Partitioning(const std::vector<double> &posit, const std::vector<double> &lengths): Generic_Domain_Partitioning(){
		if(posit.size() != lengths.size()){
			std::cout << "Dimensions of arrays are inconsistent in Generic Domain Partitioning constructor" << std::endl;
			throw("Dimensions of arrays inconsistent");
		}

		ndim = posit.size();
		first_dom = new Generic_Domain(posit, lengths);

		part_start.resize(ndim);
		part_end.resize(ndim);
		for(int i = 0; i < ndim; i++){
			part_start[i] = posit[i];
			part_end[i] = posit[i] + lengths[i];
		}
	}

	// Load from file constructor
	Generic_Domain_Partitioning(const std::string &file_name): Generic_Domain_Partitioning(){
		load_domains(file_name);
	}

	// Comparison operator
	bool operator==(const Generic_Domain_Partitioning &part) const{
		if(ndim != part.ndim) return false;
		if(part_start != part.part_start || part_end != part.part_end) return false;

		const Generic_Domain *tmp = first_dom;
		const Generic_Domain *comp_tmp = part.first_dom;
		while(tmp != NULL){
			if(comp_tmp == NULL) return false;
			if(!(*tmp == *comp_tmp)) return false;
			tmp = tmp->dom_next;
			comp_tmp = comp_tmp->dom_next;
		}
		if(tmp != comp_tmp) return false;
		else return true;
	}

	//Finds domain by point
	const Generic_Domain& get_domain(const std::vector<double> &posit) const{
		const Generic_Domain *domain_tmp = first_dom;
		while(domain_tmp != NULL){
			bool is_inside = true;
			for(int i = 0; i < ndim; i++){
				if(domain_tmp->dom_start[i] == domain_tmp->dom_end[i]){//Point cannot be in a domain of volume zero
					is_inside = false;
					break;
				}
				if(domain_tmp->dom_start[i] >= posit[i] || domain_tmp->dom_end[i] <= posit[i]){
					is_inside = false;
					break;
				}
			}
			if(is_inside) return *domain_tmp;
			else domain_tmp = domain_tmp->dom_next;
		}
		std::cout << "The point is not inside any of the domains\n";
		throw("The point is not inside any of the domains");
	}

	Generic_Domain& get_domain(const std::vector<double> &posit){
		return const_cast<Generic_Domain&>(const_cast<const Generic_Domain_Partitioning*>(this)->get_domain(posit));
	}

	//Splits domains into two at some postition after the start
	void domain_split(const std::vector<double> &posit, int split_dim){
		if(posit.size() < ndim){
			std::cout << "Domain split was coordinates with insufficient number of dimensions" << std::endl;
			throw("Domain split was coordinates with insufficient number of dimensions");
		}
		get_domain(posit).domain_split(split_dim, posit[split_dim]);
	}

	void move_domain_wall(const std::vector<double> &wall_posit, double new_posit);

	//Scales the whole partitioning along some dimension, keeps partitions_start fixed
	void scale_partitioning(double scale, int scale_dim){
		if(scale_dim >= ndim){
			std::cout << "Wrong scale dimension number" << std::endl;
			throw("Wrong scale dimension number");
		}
		double init = part_start[scale_dim];
		part_end[scale_dim] = init + (part_end[scale_dim] - init)*scale;
		for(Generic_Domain *dom = first_dom; dom != NULL; dom = dom->dom_next){
			dom->dom_start[scale_dim] = init + (dom->dom_start[scale_dim] - init)*scale;
			dom->dom_end[scale_dim] = init + (dom->dom_end[scale_dim] - init)*scale;
		}
		verify_consistency();
	}

	~Generic_Domain_Partitioning(){
		if(first_dom != NULL){
			Generic_Domain *domain_tmp = first_dom;
			while(domain_tmp->dom_next != NULL) domain_tmp = domain_tmp->dom_next;
			while(domain_tmp->dom_prev != NULL){
				domain_tmp = domain_tmp->dom_prev;
				delete domain_tmp->dom_next;
			}
			delete domain_tmp;
		}
	}

	//Prints out 2D domain configuration in a PBM file
	void print_domain_2D(const std::string &file_name, int size_x = 300, int size_y = 300, int dim_i = 0, int dim_j = 1, double *posit_st = NULL) const;
	void print_domain_3D(const std::string &folder_name, int size_x = 300, int size_y = 300, int dim_i = 0, int dim_j = 1, int dim_k = 2) const;

	//Functions that save domains
	void save_domains(const std::string &file_name) const;

	//Functions that verify domain consistency
	//Verifies that all domains are interconnected
	void verify_ordering(bool verbose = false) const;
	//Varifies domain consistency (no overlaps and no missing pieces)
	void verify_consistency(void) const;
};

//Generic function that verifies domain consistency by checking for intersections and comparing volume to total volume
template<class T>
void check_consistency(const std::vector<std::vector<T>> &begins, const std::vector<std::vector<T>> &ends, 
						const std::vector<T> &area_begin, const std::vector<T> &area_end){
	if(begins.size() != ends.size()){
		std::cout << "Wrong input data" << std::endl;
		throw("Wrong input data");
	}
	if(begins.size() == 0)
		return;
	int ndomains = begins.size();
	int ndim = begins[0].size();
	for(int i = 0; i < ndomains; i++){
		if(begins[i].size() != ndim || ends[i].size() != ndim){
			std::cout << "Wrong input data" << std::endl;
			throw("Wrong input data");
		}
	}
	if(area_begin.size() != ndim || area_end.size() != ndim){
		std::cout << "Wrong input data" << std::endl;
		throw("Wrong input data");
	}

	//Checking intersections
	for(int i = 0; i < ndomains-1; i++){
		for(int j = i+1; j < ndomains; j++){
			bool intersection = true;
			for(int dm = 0; dm < ndim; dm++)
				intersection *= ((begins[i][dm] < ends[j][dm]) && (begins[j][dm] < ends[i][dm]));
			if(intersection){
				std::cout << "Domain intersection detected between domains " << i << " and " << j << std::endl;
				throw("Domain intersection detected");
			}
		}
	}

	//Computing volume
	T vol_sum = 0;
	std::vector<T> min(begins[0]), max(ends[0]); 
	for(int i = 0; i < ndomains; i++){
		T tmp = 1;
		bool min_cand = true, max_cand = true;
		for(int dm = 0; dm < ndim; dm++){
			if(ends[i][dm] < begins[i][dm]){
				std::cout << "Along dimension " << dm << " domain end " << ends[i][dm];
				std::cout << " happens earlier than domain start " << begins[i][dm] << std::endl;
				throw("Domain end happens earlier than domain start");
			}
			tmp *= (ends[i][dm] - begins[i][dm]);
			if(min[dm] < begins[i][dm]) min_cand = false;
			if(max[dm] > ends[i][dm]) max_cand = false;
		}
		if(min_cand) min = begins[i];
		if(max_cand) max = ends[i];
		vol_sum += tmp;
	}
	T vol_total = 1.0;
	for(int dm = 0; dm < ndim; dm++){
		vol_total *= (max[dm] - min[dm]);
	}
	if(std::fabs(vol_sum - vol_total) > COMP_TOL){
		std::cout << "Area deficit detected" << std::endl;
		throw("Area deficit detected");
	}
	if(vol_total == 0){
		std::cout << "Total volume is zero" << std::endl;
		throw("Total volume is zero");
	}
	if(min != area_begin || max != area_end){
		std::cout << "Domain-partition mismatch detected" << std::endl;
		for(int dm = 0; dm < ndim; dm++){
			std::cout << "Along dimension " << dm << " domain range is [" << min[dm] << ", " << max[dm] << "]";
			std::cout << ", while partition range is [" << area_begin[dm] << ", " << area_end[dm] << "]" << std::endl;
		}
		throw("Domain-partition mismatch detected");
	}
}

#endif