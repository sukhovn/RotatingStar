#include <vector>
#include <boundaries.h>
#include <spectral.h>

#define EIGEN_SUPERLU_SUPPORT
#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#include <eigen3/Eigen/Dense>

#ifndef CONTAINER_H
#define CONTAINER_H

int count_func_files(const std::string &folder_name);
Domain_Partitioning get_partitioning_from_data(const std::string &folder_name);

class Indexer;
class Collocation_Matrix;
class Container;

class Evaluator: public Domain_Partitioning{
private:
	int nfunc;

	int npoints;
	double *data;
	std::vector<int> index_st;

	std::vector<std::vector<bool>> parity_st; //false - even, true - odd, empty - no doubling
	std::vector<std::vector<bool>> parity_end;

public:
	Evaluator(const Container &cont);

	Evaluator(Evaluator &&evali): Domain_Partitioning(std::move(evali)),
		nfunc(evali.nfunc), npoints(evali.npoints), data(evali.data), index_st(evali.index_st)
	{
		evali.data = NULL;
	}

	~Evaluator(){
		if(data != NULL) delete [] data;
	}

	double operator()(const std::vector<double> &coord, int fn) const;

	double d(int dim, const std::vector<double> &coord, int fn) const;

	std::vector<double> integrate() const;

	double test(const Container &cont, int der = -1);

	friend class Container;
};

class Container: protected Domain_Partitioning{
private:
	int nfunc;

	int npoints;
	double *data;
	std::vector<int> index_st;

	std::vector<std::vector<bool>> parity_st; //false - even, true - odd, empty - no doubling
	std::vector<std::vector<bool>> parity_end;

	void reallocate_partitioning(const Domain_Partitioning &part){
		this->Domain_Partitioning::operator=(part);
		int ndomains = 0, npoints_old = npoints;
		npoints = 0;
		index_st.resize(0);
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			index_st.push_back(npoints);
			npoints += nfunc*dom->size;
			if(dom->bulk_domain_number != ndomains++)
				std::cout << "Warning! Container numbering and partition numbering deviate" << std::endl;
		}
		if(npoints_old != npoints){
			delete [] data;
			data = new double[npoints];
		}
	}

	//Calculates coordinate based on index in simple ordering, only works when started with ind = 0
	void update_coordinate(const Bulk_Domain *dom, int ind, double *coords){
		int tmp;
		for(int dim_ch = dom->n_dim-1; dim_ch >= 0; dim_ch--){
			double dom_st = dom->domain_start(dim_ch);
			double dom_end = dom->domain_end(dim_ch);
			int dom_num = dom->dim_num(dim_ch);
			tmp = (ind % dom->dim_num(dim_ch));
			if(!parity_st[dim_ch].empty() && std::fabs(dom->domain_start(dim_ch) - partition_start(dim_ch)) < COMP_TOL)
				coords[dim_ch] = ch_x(2.0*dom_st-dom_end, dom_end, 2*dom_num-1, dom_num+tmp);
			else if(!parity_end[dim_ch].empty() && std::fabs(dom->domain_end(dim_ch) - partition_end(dim_ch)) < COMP_TOL)
				coords[dim_ch] = ch_x(dom_st, 2*dom_end-dom_st, 2*dom_num-1, tmp);
			else
				coords[dim_ch] = ch_x(dom_st, dom_end, dom_num-1, tmp);
			if(tmp) break;
			ind /= dom->dim_num(dim_ch);
		}
	}

public:
	const Domain_Partitioning &partition;

	Container(int nfunci, const Domain_Partitioning &parti):
		Domain_Partitioning(parti), partition(*this), nfunc(nfunci), npoints(0),
		parity_st(parti.n_dim), parity_end(parti.n_dim)
	{
		int ndomains = 0;
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			index_st.push_back(npoints);
			npoints += nfunc*dom->size;
			if(dom->bulk_domain_number != ndomains++)
				std::cout << "Warning! Container numbering and partition numbering deviate" << std::endl;
		}
		data = new double[npoints];
		for(int i = 0; i < npoints; i++) data[i] = 0.0;
	};

	Container(const std::string &folder_name):
		Domain_Partitioning(get_partitioning_from_data(folder_name)), partition(*this), npoints(0),
		nfunc(count_func_files(folder_name))
	{
		parity_st.resize(n_dim);
		parity_end.resize(n_dim);
		int ndomains = 0;
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			index_st.push_back(npoints);
			npoints += nfunc*dom->size;
			if(dom->bulk_domain_number != ndomains++)
				std::cout << "Warning! Container numbering and partition numbering deviate" << std::endl;
		}
		data = new double[npoints];
		load(folder_name);
	}

	Container(const Container &cont):
		Domain_Partitioning(cont), partition(*this), nfunc(cont.nfunc), npoints(cont.npoints), index_st(cont.index_st),
		parity_st(cont.parity_st), parity_end(cont.parity_end)
	{
		data = new double[npoints];
		for(int i = 0; i < npoints; i++) data[i] = cont.data[i];
	};

	void set_doubling_begin(unsigned int dim, const std::vector<bool> &parity){
		if(dim >= n_dim){
			std::cout << "Wrong dimension number" << std::endl;
			throw("Wrong domain number");
		}
		if(parity.size() != nfunc){
			std::cout << "Wrong number of functions" << std::endl;
			throw("Wrong number of functions");
		}
		if(!parity_end[dim].empty()){
			for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
				if(std::fabs(dom->domain_start(dim) - partition_start(dim)) < COMP_TOL && std::fabs(dom->domain_end(dim) - partition_end(dim)) < COMP_TOL){
					std::cout << "Cannot set doubling on both sides of the same domain" << std::endl;
					throw("Cannot set doubling on both sides of the same domain");
				}
			}
		}
		parity_st[dim] = parity;
	}

	void set_doubling_end(unsigned int dim, const std::vector<bool> &parity){
		if(dim >= n_dim){
			std::cout << "Wrong dimension number" << std::endl;
			throw("Wrong domain number");
		}
		if(parity.size() != nfunc){
			std::cout << "Wrong number of functions" << std::endl;
			throw("Wrong number of functions");
		}
		if(!parity_st[dim].empty()){
			for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
				if(std::fabs(dom->domain_start(dim) - partition_start(dim)) < COMP_TOL && std::fabs(dom->domain_end(dim) - partition_end(dim)) < COMP_TOL){
					std::cout << "Cannot set doubling on both sides of the same domain" << std::endl;
					throw("Cannot set doubling on both sides of the same domain");
				}
			}
		}
		parity_end[dim] = parity;
	}

	void reset_doubling(){
		for(std::vector<bool> &el : parity_st)
			el.clear();
		for(std::vector<bool> &el : parity_end)
			el.clear();
	}

	Container &operator=(const Container &cont){
		this->Domain_Partitioning::operator=(cont);
		nfunc = cont.nfunc;
		if(npoints != cont.npoints){
			delete [] data;
			npoints = cont.npoints;
			data = new double[npoints];
		}
		index_st = cont.index_st;
		parity_st.assign(cont.parity_st.begin(), cont.parity_st.end());
		parity_end.assign(cont.parity_end.begin(), cont.parity_end.end());
		for(int i = 0; i < npoints; i++) data[i] = cont.data[i];
		return *this;
	}

	~Container(){
		delete [] data;
	}

	int number_of_functions() const {return nfunc;}
	int size() const {return npoints;}

	double norm() const{
		double tmp = 0.0;
		for(int i = 0; i < npoints; i++){
			tmp += data[i]*data[i];
		}
		return std::sqrt(tmp);
	}

	double max_abs() const{
		double max = 0.0, tmp;
		for(int i = 0; i < npoints; i++){
			tmp = std::fabs(data[i]);
			if(tmp > max) max = tmp;
		}
		return max;
	}

	const double &operator[](int i) const {return data[i];}
	double &operator[](int i){return data[i];}

	template<class Callable>
	void fill_function(const Callable &func){
		int indx = 0;
		double *coords = new double[n_dim];
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			for(int ind = 0; ind < dom->size; ind++){
				update_coordinate(dom, ind, coords);
				for(int fn = 0; fn < nfunc; fn++)
					data[indx++] = func(fn, coords);
			}
		}
		delete [] coords;
	}

	template<class Callable>
	double compare_function(const Callable &func){
		int indx = 0;
		double *coords = new double[n_dim];
		double tmp, max_diff = 0.0;
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			for(int ind = 0; ind < dom->size; ind++){
				update_coordinate(dom, ind, coords);
				for(int fn = 0; fn < nfunc; fn++){
					tmp = std::fabs(data[indx++] - func(fn, coords));
					if(tmp > max_diff) max_diff = tmp;
				}
			}
		}
		delete [] coords;
		return max_diff;
	}

	void subtract_delta(const Collocation_Matrix &colm, double coeff = 1.0);

	void copy_to_vector(const Indexer &indexer, Eigen::VectorXd &vec) const;

	//Changes partitioning while zeroing solution
	void adjust_partitioning(const Domain_Partitioning &part_init){
		if(this->compare_bulk_domains(part_init)) return;
		Domain_Partitioning part(part_init);
		if(part.boundary_domains_present()) part.clear_boundary_domains();
		reallocate_partitioning(part);
		for(int i = 0; i < npoints; i++) data[i] = 0.0;
	}

	void save(const std::string &folder_name) const;
	void load(const std::string &folder_name);

	template<class Compactification>
	void regrid(const Domain_Partitioning &part_new, const Compactification &comp_old, const Compactification &comp_new, double value_at_infinity = 0.0){
		if(part_new.n_dim != n_dim){
			std::cout << "Regrid cannot change dimensions" << std::endl;
			throw("Regrid dimension change");
		}

		Domain_Partitioning part(part_new);
		if(part.boundary_domains_present()) part.clear_boundary_domains();

		if(this->compare_bulk_domains(part)){
			std::cout << "No container regridding is needed" << std::endl;
			return;
		}

		Evaluator eval(*this);

		reallocate_partitioning(part);

		int indx = 0, dim;
		std::vector<double> coords(n_dim), coords_old(n_dim);
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			for(int ind = 0; ind < dom->size; ind++){
				update_coordinate(dom, ind, coords.data());
				for(dim = 0; dim < n_dim; ++dim){
					coords_old[dim] = comp_new.xc(dim, coords.data());
					if(std::isinf(coords_old[dim]))
						break;
					else
						coords_old[dim] = comp_old.xc_inv(dim, coords_old.data());
				}
				if(dim < n_dim){
					for(int fn = 0; fn < nfunc; fn++)
						data[indx++] = value_at_infinity;
				}
				else{
					for(int fn = 0; fn < nfunc; fn++)
						data[indx++] = eval(coords_old, fn);
				}
			}
		}
	}

	friend class Indexer;
	friend class Collocation_Matrix;
	friend class Evaluator;
	friend class Point;
	template<class>
	friend class Point_Compactified;
	template<class>
	friend class Collocation_Solution;
	template<class, class>
	friend class Collocation_Solver_Full;
};


class Indexer{
public:
	int ndim;
	int nfunc;
	int ndomains;
	std::vector<int> index_st;
	std::vector<std::vector<int>> dim_nums;
	int npoints;

	std::vector<std::vector<double>> doms_start;
	std::vector<std::vector<double>> doms_end;

	//Parity data
	std::vector<double> part_start; //Leftmost coordinates of the leftmost domain of the partitioning
	std::vector<double> part_end; //Rightmost coordinates of the rightmost domain of the partitioning

	std::vector<std::vector<bool>> parity_st;
	std::vector<std::vector<bool>> parity_end;

public:
	Indexer(const Container &part, int nfunci): ndim(part.n_dim), nfunc(nfunci), npoints(0), ndomains(0),
		part_start(part.partition_start()), part_end(part.partition_end()), parity_st(part.parity_st), parity_end(part.parity_end)
	{
		for(const Bulk_Domain *dom = part.first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			dim_nums.push_back(dom->dim_num());
			index_st.push_back(npoints);
			npoints += nfunc*dom->size;
			doms_start.push_back(dom->domain_start());
			doms_end.push_back(dom->domain_end());
			if(dom->bulk_domain_number != ndomains++)
				std::cout << "Warning! Indexer numbering and partition numbering deviate" << std::endl;
		}
	}

	int number_of_functions() const {return nfunc;}
	int size() const {return npoints;};

	int index(unsigned int dm, const std::vector<unsigned int> &ind, unsigned int fn) const{
		if(dm >= ndomains){
			std::cout << "Wrong domain number" << std::endl;
			throw("Wrong domain number");
		}
		for(int i = 0; i < ndim; i++){
			if(ind[i] >= dim_nums[dm][i]){
				std::cout << "Wrong point number" << std::endl;
				throw("Wrong point number");
			}
		}
		if(fn >= nfunc){
			std::cout << "Wrong function number" << std::endl;
			throw("Wrong function number");
		}

		int posit = ind[0];
		for(int i = 1; i < ndim; i++){
			posit = dim_nums[dm][i]*posit + ind[i];
		}

		return index_st[dm] + nfunc*posit + fn;
	}

	std::pair<int, bool> check_parity(unsigned int dm, unsigned int dim, int fn){
		if(!parity_st[dim].empty() && std::fabs(doms_start[dm][dim] - part_start[dim]) < COMP_TOL)
			return std::make_pair(1, parity_st[dim][fn]);
		else if(!parity_end[dim].empty() && std::fabs(doms_end[dm][dim] - part_end[dim]) < COMP_TOL)
			return std::make_pair(2, parity_end[dim][fn]);
		else
			return std::make_pair(0, false);
	}

	void indexer_regrid(const Container &part){
		if(part.n_dim != ndim){
			std::cout << "Regrid cannot change dimensions" << std::endl;
			throw("Regrid dimension change");
		}

		npoints = 0;
		ndomains = 0;
		index_st.resize(0);
		dim_nums.resize(0);
		doms_start.resize(0);
		doms_end.resize(0);
		for(const Bulk_Domain *dom = part.first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			dim_nums.push_back(dom->dim_num());
			index_st.push_back(npoints);
			npoints += nfunc*dom->size;
			doms_start.push_back(dom->domain_start());
			doms_end.push_back(dom->domain_end());
			if(dom->bulk_domain_number != ndomains++)
				std::cout << "Warning! Indexer numbering and partition numbering deviate" << std::endl;
		}
		part_start = part.partition_start();
		part_end = part.partition_end();
		parity_st = part.parity_st;
		parity_end = part.parity_end;
	}
};


class Collocation_Matrix: public Indexer{
public:
	Eigen::MatrixXd matr;
	Eigen::VectorXd rhs;
	Eigen::VectorXd delta;

	std::vector<std::vector<bool>> parity_st; //false - even, true - odd, empty - no doubling
	std::vector<std::vector<bool>> parity_end;

public:
	Collocation_Matrix(const Container &cont): Indexer(cont, cont.nfunc),
		parity_st(cont.parity_st), parity_end(cont.parity_end)
	{
		matr.resize(npoints, npoints);
		rhs.resize(npoints);
		delta.resize(npoints);
		matr.setZero();
		rhs.setZero();
		delta.setZero();
	}

	//This functions zeroes the matrix
	void setZero_jacobian(void){
		matr.setZero();
	}

	void setZero_rhs(void){
		rhs.setZero();
	}

	//This function zeroes the matrix row
	void setZero_row(unsigned int dm, const std::vector<unsigned int> &ind, unsigned int fn){
		matr.row(index(dm, ind, fn)).setZero();
	}

	//This function acceses right hand side array
	double &right_hand_side(unsigned int dm, const std::vector<unsigned int> &ind, unsigned int fn){
		return rhs(index(dm, ind, fn));
	}

	//This function adds a differention operator to the (row_dm, row_ind) row, multiplying it by coef
	//der1 and der2 indicate derivative dimensions, der1 = -1 and/or der2 = -1 are used to indicate no differention
	//This is the main function used to construct the Jacobian
	void row_add_Dn_factor(unsigned int row_dm, const std::vector<unsigned int> &row_ind, unsigned int row_fn,
		int der1, int der2, unsigned int dm, const std::vector<unsigned int> &ind, unsigned int fn, double coef);


	//Computes the next step given the Jacobian and the equation difference
	void compute_delta(void){
		delta = matr.partialPivLu().solve(rhs);
	}

	//Changes domains configuration
	void regrid(const Container &part){
		indexer_regrid(part);
		matr.resize(npoints, npoints);
		rhs.resize(npoints);
		delta.resize(npoints);
		matr.setZero();
		rhs.setZero();
		delta.setZero();
	}

	//Testing functions
	//Makes matrix a differentiation matrix
	void fill_Dn_matrix(int der1, int der2){
		matr.setZero();
		std::vector<unsigned int> ind(ndim, 0);
		for(int dm = 0; dm < ndomains; ++dm){
			int indx_max = std::accumulate(dim_nums[dm].begin(), dim_nums[dm].end(), 1, std::multiplies<double>());
			for(int indx = 0; indx < indx_max; indx++){
				for(int i = ndim-1, tmpi = indx; i >= 0; i--){
					if((ind[i] = tmpi % dim_nums[dm][i]) > 1) break;
					tmpi /= dim_nums[dm][i];
				}
				for(int fn = 0; fn < nfunc; ++fn)
					row_add_Dn_factor(dm, ind, fn, der1, der2, dm, ind, fn, 1.0);
			}
		}
	}
};

#endif