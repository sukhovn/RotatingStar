#include <domains.h>
#include <numeric>
#include <algorithm>
#include <unordered_set>

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

class Bulk_Domain: public Generic_Domain{
private:
	Bulk_Domain *dom_bulk_prev;
	Bulk_Domain *dom_bulk_next;
	std::vector<int> dim_n;
	int sz;

	//For testing purposes;
	int b_dom_number;

	//Default constructor
	Bulk_Domain(): Generic_Domain(), dom_bulk_prev(NULL), dom_bulk_next(NULL),
		sz(0), bulk_domain_number(b_dom_number), size(sz), b_dom_number(-1) {};

	// Copy constructor from Generic_Domain
	Bulk_Domain(const Generic_Domain &g_dom): Generic_Domain(g_dom),
		dom_bulk_prev(NULL), dom_bulk_next(NULL), dim_n(g_dom.ndim, 0), sz(0),
		bulk_domain_number(b_dom_number), size(sz), b_dom_number(-1)
	{
		const Bulk_Domain *b_dom_tmp = dynamic_cast<const Bulk_Domain*>(&g_dom);
		if(b_dom_tmp != NULL){
			dim_n = b_dom_tmp->dim_n;
			sz = b_dom_tmp->sz;
		}
	}

	// Copy constructor
	Bulk_Domain(const Bulk_Domain &b_dom): Generic_Domain(b_dom), 
		dom_bulk_prev(NULL), dom_bulk_next(NULL), dim_n(b_dom.dim_n), sz(b_dom.sz),
		bulk_domain_number(b_dom_number), size(sz), b_dom_number(b_dom.b_dom_number) {};

	// Constructs domain from domain file description
	Bulk_Domain(const std::vector<std::string> &instruct);

	// Constructs domain from a single line description
	Bulk_Domain(const std::string &instruct);

	// Consturcts domain from direct specification
	Bulk_Domain(const std::vector<double> &posit, const std::vector<double> &lengths, const std::vector<int> &dim_nums): Generic_Domain(posit, lengths),
		dom_bulk_prev(NULL), dom_bulk_next(NULL), dim_n(dim_nums), sz(1),
		bulk_domain_number(b_dom_number), size(sz), b_dom_number(-1)
	{
		if(lengths.size() != dim_nums.size()){
			std::cout << "Dimensions of arrays are inconsistent in Bulk Domain constructor" << std::endl;
			throw("Dimensions of arrays inconsistent");
		}

		for(int i = 0; i < ndim; i++){
			sz *= dim_n[i];
		}
	}

	// Splits domains into two at some postition after the start
	void domain_split(int split_dim, double split_coord){
		if(split_dim >= ndim){
			std::cout << "Invalid domain splitting dimension\n";
			throw("Invalid domain splitting dimension");
		}
		if(split_coord < dom_start[split_dim] || split_coord > dom_end[split_dim]){
			std::cout << "Invalid domain splitting position\n";
			throw("Invalid domain splitting");
		}

		Bulk_Domain *new_dom = new Bulk_Domain(*this);
		new_dom->b_dom_number = b_dom_number+1;

		new_dom->dom_next = dom_next;
		new_dom->dom_bulk_next = dom_bulk_next;
		if(dom_next != NULL){
			dom_next->dom_prev = new_dom;
		}
		if(dom_bulk_next != NULL){
			dom_bulk_next->dom_bulk_prev = new_dom;
			for(Bulk_Domain *tmp = dom_bulk_next; tmp != NULL; tmp = tmp->dom_bulk_next)
				tmp->b_dom_number++;
		}

		new_dom->dom_prev = new_dom->dom_bulk_prev = this;
		dom_next = dom_bulk_next = new_dom;

		new_dom->dom_start[split_dim] = split_coord;
		dom_end[split_dim] = new_dom->dom_start[split_dim];		
	}

	friend class Domain_Partitioning;

public:
	const int &bulk_domain_number;
	const int &size;

	//Outputs a pointer to previous and next bulk domains
	Bulk_Domain *previous_bulk_domain(){return dom_bulk_prev;}
	const Bulk_Domain *previous_bulk_domain() const {return dom_bulk_prev;}
	Bulk_Domain *next_bulk_domain(){return dom_bulk_next;}
	const Bulk_Domain *next_bulk_domain() const {return dom_bulk_next;}
	int dim_num(int dm) const {return dim_n[dm];}
	const std::vector<int> &dim_num() const {return dim_n;}

	// Prints domain starting and ending coordinates along every dimension
	void print_domain(std::ostream& os) const{
		for(int dm = 0; dm < ndim; dm++){
			os << "Dimension: " << dm << ", start: " << dom_start[dm] << ", end: " << dom_end[dm];
			if(dim_n[dm] > 0) os << ", points: " << dim_n[dm];
			os << std::endl;
		}
	}

	// Comparison operator
	bool operator==(const Generic_Domain &dom) const{
		const Bulk_Domain *tmp = dynamic_cast<const Bulk_Domain *>(&dom);
		if(tmp == NULL) return false;
		if(!(this->Generic_Domain::operator==(*tmp))) return false;
		if(dim_n != tmp->dim_n) return false;
		if(sz != tmp->sz) return false;
		return true;
	}

};

class Boundary_Domain: public Generic_Domain{
private:
	int bound_dim;
	Bulk_Domain *dom_left_side;
	Bulk_Domain *dom_right_side;

	std::vector<int> left_n_st;
	std::vector<int> left_n_end;
	std::vector<int> right_n_st;
	std::vector<int> right_n_end;

	//Creates an edge boundary domain
	Boundary_Domain(Bulk_Domain &b_dom, int edge_dim, bool left = true): boundary_dimension(bound_dim),
	Generic_Domain(b_dom.ndim), bound_dim(edge_dim), dom_left_side(NULL), dom_right_side(NULL)
	{
		dom_start = b_dom.dom_start;
		dom_end = b_dom.dom_end;

		if(left){
			dom_end[bound_dim] = dom_start[bound_dim];
			dom_right_side = &b_dom;

			right_n_st.resize(ndim-1);
			right_n_end.resize(ndim-1);
			for(int dm = 0; dm < ndim-1; dm++){right_n_st[dm] = right_n_end[dm] = 0;}
		}
		else{
			dom_start[bound_dim] = dom_end[bound_dim];
			dom_left_side = &b_dom;

			left_n_st.resize(ndim-1);
			left_n_end.resize(ndim-1);
			for(int dm = 0; dm < ndim-1; dm++){left_n_st[dm] = left_n_end[dm] = 0;}
		}
	}

	//Creates a boundary domain between two bulk domains
	Boundary_Domain(Bulk_Domain &dom1, Bulk_Domain &dom2): boundary_dimension(bound_dim),
	Generic_Domain(dom1.ndim), bound_dim(dom1.ndim), dom_left_side(NULL), dom_right_side(NULL),
	left_n_st(dom1.ndim-1, 0), left_n_end(dom1.ndim-1, 0), right_n_st(dom1.ndim-1, 0), right_n_end(dom1.ndim-1, 0)
	{
		if(dom1.ndim != dom2.ndim){
			std::cout << "Domains dimensions don't match in Boundary_Domain constructor" << std::endl;
			throw("Domains dimensions don't match in Boundary_Domain constructor");
		}

		for(int dm = 0; dm < ndim; dm++){
			if(dom1.dom_end[dm] == dom2.dom_start[dm] || dom1.dom_start[dm] == dom2.dom_end[dm]){
				if(bound_dim == ndim) bound_dim = dm;
				else{
					std::cout << "Warning! At least two possible boundary dimensions in Boundary_Domain constructor" << std::endl;
					break;
				}
			}
		}
		if(bound_dim == ndim){
			std::cout << "No possible boundary dimension found in Boundary_Domain constructor" << std::endl;
			throw("No possible boundary dimension found in Boundary_Domain constructor");
		}
		if(dom1.dom_end[bound_dim] == dom2.dom_start[bound_dim]){
			dom_left_side = &dom1; dom_right_side = &dom2;
		}
		else{
			dom_left_side = &dom2; dom_right_side = &dom1;
		}

		for(int dm = 0; dm < ndim; dm++) dom_start[dm] = MAX(dom_left_side->dom_start[dm], dom_right_side->dom_start[dm]);
		for(int dm = 0; dm < ndim; dm++) dom_end[dm] = MIN(dom_left_side->dom_end[dm], dom_right_side->dom_end[dm]);
	}

	// Copy constructor
	Boundary_Domain(const Boundary_Domain &b_dom): Generic_Domain(b_dom), boundary_dimension(bound_dim),
		bound_dim(b_dom.bound_dim), dom_left_side(NULL), dom_right_side(NULL),
		left_n_st(b_dom.left_n_st), left_n_end(b_dom.left_n_end), right_n_st(b_dom.right_n_st), right_n_end(b_dom.right_n_end) {};

	friend class Domain_Partitioning;

public:
	const int &boundary_dimension;

	// Outputs a pointer to previous and next domains
	Bulk_Domain *left_side_domain(){return dom_left_side;}
	const Bulk_Domain *left_side_domain() const {return dom_left_side;}
	Bulk_Domain *right_side_domain(){return dom_right_side;}
	const Bulk_Domain *right_side_domain() const {return dom_right_side;}

	// Prints domain starting and ending coordinates along every dimension
	void print_domain(std::ostream& os) const{
		os << "Boundary dimension: " << bound_dim << std::endl;
		if(dom_left_side == NULL) os << "This boundary domain is adjacent to bulk domain " << dom_right_side->bulk_domain_number << " from the left\n";
		else if(dom_right_side == NULL) os << "This boundary domain is adjacent to bulk domain " << dom_left_side->bulk_domain_number << " from the right\n";
		else{
			os << "This boundary domain is in between bulk domain " << dom_left_side->bulk_domain_number;
			os << " and " << dom_right_side->bulk_domain_number << std::endl;
		}
		for(int dm = 0, ii = 0; dm < ndim; dm++){
			os << "Dimension: " << dm << ", start: " << dom_start[dm] << ", end: " << dom_end[dm] << std::endl;
			if(dm == bound_dim) continue;
			if(left_n_end.size() > 0)
				if(left_n_end[ii] != 0)
					os << "Left point range along dimension " << dm << " is from " << left_n_st[ii] << " to " << left_n_end[ii] << std::endl;
			if(right_n_end.size() > 0)
				if(right_n_end[ii] != 0)
					os << "Right point range along dimension " << dm << " is from " << right_n_st[ii] << " to " << right_n_end[ii] << std::endl;
			ii++;
		}
	}

	// Comparison operator
	bool operator==(const Generic_Domain &dom) const{
		const Boundary_Domain *tmp = dynamic_cast<const Boundary_Domain *>(&dom);
		if(tmp == NULL) return false;
		if(!(this->Generic_Domain::operator==(*tmp))) return false;
		if(bound_dim != tmp->bound_dim) return false;
		if(left_n_st != tmp->left_n_st) return false;
		if(left_n_end != tmp->left_n_end) return false;
		if(right_n_st != tmp->right_n_st) return false;
		if(right_n_end != tmp->right_n_end) return false;
		if(dom_left_side != NULL){
			if(tmp->dom_left_side == NULL) return false;
			else if(dom_left_side->bulk_domain_number != (tmp->dom_left_side)->bulk_domain_number) return false;
		}
		if(dom_right_side != NULL){
			if(tmp->dom_right_side == NULL) return false;
			else if(dom_right_side->bulk_domain_number != (tmp->dom_right_side)->bulk_domain_number) return false;
		}
		return true;
	}

	int left_num_st(int dm) const {return left_n_st[dm];}
	int left_num_end(int dm) const {return left_n_end[dm];}
	int right_num_st(int dm) const {return right_n_st[dm];}
	int right_num_end(int dm) const {return right_n_end[dm];}
	const std::vector<int> &left_num_st() const {return left_n_st;}
	const std::vector<int> &left_num_end() const {return left_n_end;}
	const std::vector<int> &right_num_st() const {return right_n_st;}
	const std::vector<int> &right_num_end() const {return right_n_end;}
};

class Domain_Partitioning: public Generic_Domain_Partitioning{
private:
	Bulk_Domain *first_bulk_dom;

	bool boundary_domains_inserted;

	Domain_Partitioning():  Generic_Domain_Partitioning(), first_bulk_domain(first_bulk_dom), first_bulk_dom(NULL), boundary_domains_inserted(false) {};

	virtual Generic_Domain* load_single_domain(const std::vector<std::string> &domain_info, Generic_Domain *prev_dom);
	void load_domains(const std::string &file_name);
	void load_domains(const std::vector<std::string> &instruct);

public:
	//Returns first bulk domain
	Bulk_Domain* const &first_bulk_domain;
	bool boundary_domains_present(void){return boundary_domains_inserted;}

	//Adds boundary domains around and between bulk domains
	void insert_boundary_domains(const std::unordered_set<int> &exclude_left = {}, const std::unordered_set<int> &exclude_right = {}){
		if(boundary_domains_inserted){
			std::cout << "Boundary domains are already inserted" << std::endl;
			return;
		}

		//Preparing an included dimensions list
		std::vector<int> dim_left(ndim);
		std::vector<int> dim_right(ndim);
		std::iota(dim_left.begin(), dim_left.end(), 0);
		std::iota(dim_right.begin(), dim_right.end(), 0);
		dim_left.erase(std::remove_if(dim_left.begin(), dim_left.end(), [&exclude_left](int i){return exclude_left.count(i);}), dim_left.end());
		dim_right.erase(std::remove_if(dim_right.begin(), dim_right.end(), [&exclude_right](int i){return exclude_right.count(i);}), dim_right.end());

		Generic_Domain *dom_tmp;
		Boundary_Domain *bound_dom_tmp;
		Bulk_Domain *bulk_dom_tmp = first_bulk_dom;
		while(bulk_dom_tmp != NULL){
			//Making left boundary domains
			dom_tmp = bulk_dom_tmp->dom_prev;
			bound_dom_tmp = NULL;

			for(int dm : dim_left){
				if(bulk_dom_tmp->dom_start[dm] == part_start[dm]){
					bound_dom_tmp = new Boundary_Domain(*bulk_dom_tmp, dm, true);
					bound_dom_tmp->dom_right_side = bulk_dom_tmp;
					bound_dom_tmp->dom_prev = dom_tmp;
					if(dom_tmp != NULL)
						dom_tmp->dom_next = bound_dom_tmp;
					dom_tmp = bound_dom_tmp;
				}
			}
			if(bound_dom_tmp != NULL){
				bound_dom_tmp->dom_next = bulk_dom_tmp;
				bulk_dom_tmp->dom_prev = bound_dom_tmp;
			}

			//Making cross-domain boundaries
			Bulk_Domain *bulk_dom_iter = bulk_dom_tmp->dom_bulk_next;
			dom_tmp = bulk_dom_tmp;
			bound_dom_tmp = NULL;
			while(bulk_dom_iter != NULL){
				if(check_adjacency(bulk_dom_tmp, bulk_dom_iter)){
					bound_dom_tmp = new Boundary_Domain(*bulk_dom_tmp, *bulk_dom_iter);
					bound_dom_tmp->dom_prev = dom_tmp;
					dom_tmp->dom_next = bound_dom_tmp;
					dom_tmp = bound_dom_tmp;
				}
				bulk_dom_iter = bulk_dom_iter->dom_bulk_next;
			}
			if(bound_dom_tmp != NULL){
				dom_tmp->dom_next = bulk_dom_tmp->dom_bulk_next;
				(bulk_dom_tmp->dom_bulk_next)->dom_prev = dom_tmp;
			}

			//Making right boundary domains
			bound_dom_tmp = NULL;
			for(int dm : dim_right){
				if(bulk_dom_tmp->dom_end[dm] == part_end[dm]){
					bound_dom_tmp = new Boundary_Domain(*bulk_dom_tmp, dm, false);
					bound_dom_tmp->dom_left_side = bulk_dom_tmp;
					bound_dom_tmp->dom_prev = dom_tmp;
					dom_tmp->dom_next = bound_dom_tmp;
					dom_tmp = bound_dom_tmp;
				}
			}
			if(bound_dom_tmp != NULL){
				dom_tmp->dom_next = bulk_dom_tmp->dom_bulk_next;
				if(bulk_dom_tmp->dom_bulk_next != NULL)
					(bulk_dom_tmp->dom_bulk_next)->dom_prev = dom_tmp;
			}

			//Going to the next bulk domain
			bulk_dom_tmp = bulk_dom_tmp->dom_bulk_next;
		}

		dom_tmp = first_bulk_dom;
		while(dom_tmp->dom_prev != NULL){
			dom_tmp = dom_tmp->dom_prev;
		}
		first_dom = dom_tmp;

		verify_ordering();
		verify_consistency();

		boundary_domains_inserted = true;
	}

	//Reverts previous operation and removes all boundary domains
	void clear_boundary_domains(void){
		if(!boundary_domains_inserted){
			std::cout << "Boundary domains indicated as abscent" << std::endl;
		}

		Bulk_Domain *bulk_tmp;
		for(Generic_Domain* dom_tmp = first_dom; dom_tmp != NULL; dom_tmp = dom_tmp->dom_next){
			if(dynamic_cast<Boundary_Domain *>(dom_tmp) != NULL)
				delete dom_tmp;
			else if((bulk_tmp = dynamic_cast<Bulk_Domain *>(dom_tmp)) != NULL){
				bulk_tmp->dom_prev = bulk_tmp->dom_bulk_prev;
				if(bulk_tmp->dom_bulk_prev != NULL) (bulk_tmp->dom_bulk_prev)->dom_next = bulk_tmp;
			}
		}
		first_dom = first_bulk_dom;
		bulk_tmp->dom_next = NULL;

		verify_ordering();

		boundary_domains_inserted = false;
	}

	Domain_Partitioning(const Generic_Domain_Partitioning &part): 
		first_bulk_domain(first_bulk_dom), boundary_domains_inserted(false){
		const Domain_Partitioning *dom_part = dynamic_cast<const Domain_Partitioning*>(&part);
		if(dom_part != NULL){
			*this = part;
			return;
		}

		ndim = part.ndim;
		part_start = part.part_start;
		part_end = part.part_end;

		const Generic_Domain* part_old_tmp = part.first_dom;
		Bulk_Domain *part_new_prev = NULL, *part_new_curr;

		int domain_number = 0;
		while(part_old_tmp != NULL){
			if(!(part_old_tmp->is_zero_volume())){
				part_new_curr = new Bulk_Domain(*part_old_tmp);				
				part_new_curr->b_dom_number = domain_number++;

				if(part_new_prev != NULL){//Linking domains
					part_new_curr->dom_prev = part_new_curr->dom_bulk_prev = part_new_prev;
					part_new_prev->dom_next = part_new_prev->dom_bulk_next = part_new_curr;
				}
				else first_dom = first_bulk_dom = part_new_curr;

				part_new_prev = part_new_curr;
			}
			part_old_tmp = part_old_tmp->dom_next;
		}
	}

	Domain_Partitioning &operator=(const Domain_Partitioning &part){
		boundary_domains_inserted = part.boundary_domains_inserted;
		ndim = part.ndim;
		part_start = part.part_start;
		part_end = part.part_end;

		Generic_Domain *dom_curr = first_dom, *dom_tmp;
		first_dom = NULL;
		while(dom_curr != NULL){
			dom_tmp = dom_curr;
			dom_curr = dom_curr->dom_next;
			delete dom_tmp;
		}

		Bulk_Domain *bulk_curr = NULL, *bulk_tmp;
		Boundary_Domain *bound_tmp;
		for(const Generic_Domain* part_old_tmp = part.first_dom; part_old_tmp != NULL; part_old_tmp = part_old_tmp->dom_next){
			if(dynamic_cast<const Bulk_Domain *>(part_old_tmp) != NULL){
				dom_tmp = bulk_tmp = new Bulk_Domain(*dynamic_cast<const Bulk_Domain *>(part_old_tmp));
				if(bulk_curr == NULL){
					first_bulk_dom = bulk_curr = bulk_tmp;
				}
				else{
					bulk_tmp->dom_bulk_prev = bulk_curr;
					bulk_curr->dom_bulk_next = bulk_tmp;
					bulk_curr = bulk_tmp;
				} 
			}
			else if(dynamic_cast<const Boundary_Domain *>(part_old_tmp) != NULL)
				dom_tmp = new Boundary_Domain(*dynamic_cast<const Boundary_Domain *>(part_old_tmp));
			else{
				std::cout << "Unrecognized domain type in copy constructor" << std::endl;
				throw("Unrecognized domain type in copy constructor");
			}
			if(dom_curr == NULL){
				first_dom = dom_curr = dom_tmp;
			}
			else{
				dom_tmp->dom_prev = dom_curr;
				dom_curr->dom_next = dom_tmp;
				dom_curr = dom_tmp;
			}
		}

		dom_tmp = first_dom;
		const Boundary_Domain *bound_const_tmp;
		for(const Generic_Domain* part_old_tmp = part.first_dom; part_old_tmp != NULL; part_old_tmp = part_old_tmp->dom_next, dom_tmp = dom_tmp->dom_next){
			if((bound_const_tmp = dynamic_cast<const Boundary_Domain *>(part_old_tmp)) == NULL) continue;
			bound_tmp = dynamic_cast<Boundary_Domain *>(dom_tmp);
			if(bound_const_tmp->dom_left_side != NULL){
				for(bulk_tmp = first_bulk_dom; bulk_tmp != NULL; bulk_tmp = bulk_tmp->dom_bulk_next){
					if(bulk_tmp->b_dom_number == (bound_const_tmp->dom_left_side)->b_dom_number){
						bound_tmp->dom_left_side = bulk_tmp;
					}
				}
			}
			if(bound_const_tmp->dom_right_side != NULL){
				for(bulk_tmp = first_bulk_dom; bulk_tmp != NULL; bulk_tmp = bulk_tmp->dom_bulk_next){
					if(bulk_tmp->b_dom_number == (bound_const_tmp->dom_right_side)->b_dom_number){
						bound_tmp->dom_right_side = bulk_tmp;
					}
				}
			}
		}

		verify_ordering();
		verify_consistency();
		return *this;
	}

	Domain_Partitioning &operator=(Domain_Partitioning &&part){
		this->Generic_Domain_Partitioning::operator=(std::move(part));
		boundary_domains_inserted = part.boundary_domains_inserted;
		first_bulk_dom = part.first_bulk_dom;
		part.first_bulk_dom = NULL;

		return *this;
	}

	//Copy constructor
	Domain_Partitioning(const Domain_Partitioning &part): Generic_Domain_Partitioning(), first_bulk_domain(first_bulk_dom), 
		boundary_domains_inserted(part.boundary_domains_inserted)
	{
		*this = part;
	}

	// Move constructor
	Domain_Partitioning(Domain_Partitioning &&part): Generic_Domain_Partitioning(std::move(part)),
		first_bulk_domain(first_bulk_dom), first_bulk_dom(part.first_bulk_dom), boundary_domains_inserted(part.boundary_domains_inserted)
	{
		part.first_bulk_dom = NULL;
	}

	//Creates a single domain partitioning
	Domain_Partitioning(const std::vector<double> &posit, const std::vector<double> &lengths, const std::vector<int> &dim_nums):
		Domain_Partitioning()
	{
		if(posit.size() != lengths.size() || lengths.size() != dim_nums.size()){
			std::cout << "Dimensions of arrays are inconsistent in Domain Partitioning constructor" << std::endl;
			throw("Dimensions of arrays inconsistent");
		}
		first_dom = first_bulk_dom = new Bulk_Domain(posit, lengths, dim_nums);
		first_bulk_dom->b_dom_number = 0;

		ndim = posit.size();
		part_start.resize(ndim);
		part_end.resize(ndim);
		for(int i = 0; i < ndim; i++){
			part_start[i] = posit[i];
			part_end[i] = posit[i] + lengths[i];
		}
	}

	//Loads boundary domain partitioning
	Domain_Partitioning(const std::string &file_name): Domain_Partitioning(){
		load_domains(file_name);
	}

	Domain_Partitioning(const std::vector<std::string> &instruct): Domain_Partitioning(){
		load_domains(instruct);
	}

	bool operator==(const Domain_Partitioning &part) const{
		if(ndim != part.ndim) return false;
		if(part_start != part.part_start || part_end != part.part_end) return false;

		const Generic_Domain *tmp = first_dom;
		const Generic_Domain *comp_tmp = part.first_dom;
		for(; tmp != NULL; tmp = tmp->dom_next, comp_tmp = comp_tmp->dom_next){
			if(comp_tmp == NULL) return false;
			if(!(*tmp == *comp_tmp)) return false;
		}
		if(tmp != comp_tmp) return false;
		else return true;
	}

	bool compare_bulk_domains(const Domain_Partitioning &part) const{
		if(ndim != part.ndim) return false;
		if(part_start != part.part_start || part_end != part.part_end) return false;

		const Bulk_Domain *bulk_tmp = first_bulk_dom, *bulk_comp_tmp = part.first_bulk_dom;
		for(; bulk_tmp != NULL; bulk_tmp = bulk_tmp->dom_bulk_next, bulk_comp_tmp = bulk_comp_tmp->dom_bulk_next){
			if(bulk_comp_tmp == NULL) return false;
			if(!(*bulk_tmp == *bulk_comp_tmp)) return false;
		}
		if(bulk_tmp != bulk_comp_tmp) return false;
		else return true;
	}

	//Finds bulk domain by point
	const Bulk_Domain* operator()(const std::vector<double> &posit) const{
		for(const Bulk_Domain *domain_tmp = first_bulk_dom; domain_tmp != NULL; domain_tmp = domain_tmp->dom_bulk_next){
			bool is_inside = true;
			for(int i = 0; i < ndim; i++){
				if(domain_tmp->dom_start[i] > posit[i] || domain_tmp->dom_end[i] < posit[i]){
					is_inside = false;
					break;
				}
			}
			if(is_inside) return domain_tmp;
		}
		std::cout << "The point is not inside any of the domains\n";
		throw("The point is not inside any of the domains");
	}

	Bulk_Domain* operator()(const std::vector<double> &posit){
		return const_cast<Bulk_Domain*>((*const_cast<const Domain_Partitioning*>(this))(posit));
	}

	//Splits domains into two at some postition after the start
	void domain_split(const std::vector<double> &posit, int split_dim){
		if(boundary_domains_inserted){
			std::cout << "Domain split only works on partitions with no boundary domains" << std::endl;
			throw("Domain split only works on partitions with no boundary domains");
		}
		Generic_Domain_Partitioning::domain_split(posit, split_dim);
	}

	void move_domain_wall(const std::vector<double> &wall_posit, double new_posit){
		if(boundary_domains_inserted){
			std::cout << "Domain wall move only works on partitions with no boundary domains" << std::endl;
			throw("Domain wall move only works on partitions with no boundary domains");
		}
		Generic_Domain_Partitioning::move_domain_wall(wall_posit, new_posit);
	}

	void change_domain_resolution(const std::vector<double> &wall_posit, const std::vector<int> &new_n){
		if(boundary_domains_inserted){
			std::cout << "Domain resolution can only be changed on partitions with no boundary domains" << std::endl;
			throw("Domain resolution can only be changed on partitions with no boundary domains");
		}
		Bulk_Domain *dom = (*this)(wall_posit);
		dom->dim_n = new_n;
		int sz = 1;
		for(int i = 0; i < ndim; i++) sz *= new_n[i];
		dom->sz = sz;
	}

	//Distributes domain points along the boundaries using a function func(start, end, n, position)
	//that calculates point split position based on boundary start, end, number of points and position of the split
	template<class Cut_Func>
	void fill_boundary_point_range(const Cut_Func &func, const std::unordered_set<int> &exclude_left = {}, const std::unordered_set<int> &exclude_right = {}){
		if(!boundary_domains_inserted){
			std::cout << "Warning! Trying to fill boundary point splitting for a partition with no boundaries" << std::endl;
			return;
		}

		const Bulk_Domain *bulk;
		int dim;

		for(Generic_Domain *g_dom = first_domain; g_dom != NULL; g_dom = g_dom->next_domain()){
			Boundary_Domain *dom = dynamic_cast<Boundary_Domain *>(g_dom);
			if(dom == NULL) continue;

			for(int dim = 0, indx = 0; dim < ndim; dim++){
				if(dom->boundary_dimension == dim) continue;

				if((bulk = dom->left_side_domain()) != NULL){
					if(bulk->domain_start(dim) == dom->domain_start(dim))
						if(dim < dom->boundary_dimension && !(exclude_left.count(dim) && std::fabs(part_start[dim] - dom->dom_start[dim]) < COMP_TOL))
							dom->left_n_st[indx] = 1;
						else
							dom->left_n_st[indx] = 0;
					else
						dom->left_n_st[indx] = func(dim, bulk->domain_start(dim), bulk->domain_end(dim), bulk->dim_n[dim], dom->domain_start(dim)) + 1;
					if(bulk->domain_end(dim) == dom->domain_end(dim))
						if(dim < dom->boundary_dimension && !(exclude_right.count(dim) && std::fabs(part_end[dim] - dom->dom_end[dim]) < COMP_TOL))
							dom->left_n_end[indx] = bulk->dim_n[dim]-2;
						else
							dom->left_n_end[indx] = bulk->dim_n[dim]-1;
					else
						dom->left_n_end[indx] = func(dim, bulk->domain_start(dim), bulk->domain_end(dim), bulk->dim_n[dim], dom->domain_end(dim));
				}
				if((bulk = dom->right_side_domain()) != NULL){
					if(bulk->domain_start(dim) == dom->domain_start(dim))
						if(dim < dom->boundary_dimension && !(exclude_left.count(dim) && std::fabs(part_start[dim] - dom->dom_start[dim]) < COMP_TOL))
							dom->right_n_st[indx] = 1;
						else
							dom->right_n_st[indx] = 0;
					else
						dom->right_n_st[indx] = func(dim, bulk->domain_start(dim), bulk->domain_end(dim), bulk->dim_n[dim], dom->domain_start(dim)) + 1;
					if(bulk->domain_end(dim) == dom->domain_end(dim))
						if(dim < dom->boundary_dimension && !(exclude_right.count(dim) && std::fabs(part_end[dim] - dom->dom_end[dim]) < COMP_TOL))
							dom->right_n_end[indx] = bulk->dim_n[dim]-2;
						else
							dom->right_n_end[indx] = bulk->dim_n[dim]-1;
					else
						dom->right_n_end[indx] = func(dim, bulk->domain_start(dim), bulk->domain_end(dim), bulk->dim_n[dim], dom->domain_end(dim));
				}
				indx++;
			}
		}

		verify_boundary_point_splitting(exclude_left, exclude_right);
	}

	//Functions that verify domain consistency
	//Verifies that all the domains are interconnected
	void verify_ordering(bool verbose = false);
	//Verifies that the boundaries have points consistent with bulk domains
	void verify_boundary_point_splitting(const std::unordered_set<int> &exclude_left, const std::unordered_set<int> &exclude_right);
	
	//Functions that load and save domains
	void save_domains(const std::string &file_name);
};

//Debug functions
void print_partitioning(const Generic_Domain_Partitioning &part);

#endif