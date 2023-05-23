#include <fftw3.h>
#include <container.h>
#include <point_class.h>

#ifndef COLLOCATION_SOLUTION_H
#define COLLOCATION_SOLUTION_H

template<class Compactification>
class Collocation_Solution: public Domain_Partitioning{
protected:
	Point_Compactified<Compactification> pnt;
	Container *cont;

public:
	void fill_boundary_domains(Domain_Partitioning &part){
		std::unordered_set<int> exclude_left;
		std::unordered_set<int> exclude_right;
		for(int dim = 0; dim < cont->n_dim; ++dim){
			if(!cont->parity_st[dim].empty()) exclude_left.insert(dim);
			if(!cont->parity_end[dim].empty()) exclude_right.insert(dim);
		}
		part.insert_boundary_domains(exclude_left, exclude_right);

		auto boundary_function = [this](int dim, double xst, double xend, int n, double x) -> int{
			if(!cont->parity_st[dim].empty() && std::fabs(xst - cont->partition_start(dim)) < COMP_TOL){
				return ch_x_inv_low(2*xst-xend, xend, 2*n-1, x)-n;
			}
			else if(!cont->parity_end[dim].empty() && std::fabs(xend -  cont->partition_end(dim)) < COMP_TOL){
				return ch_x_inv_low(xst, 2*xend-xst, 2*n-1, x);
			}
			else{
				return ch_x_inv_low(xst, xend, n-1, x);
			}
		};
		part.fill_boundary_point_range(boundary_function, exclude_left, exclude_right);
	}

public:
	const Domain_Partitioning &partition(){return *this;}

	Collocation_Solution(const Compactification &compact, Container *conti):
		Domain_Partitioning(*conti), cont(conti), pnt(compact, conti) {};

	bool verify_container_solution_correspondence(void){
		if(compare_bulk_domains(*cont)) return false;
		Domain_Partitioning part(*cont);
		if(part.boundary_domains_present())
			part.clear_boundary_domains();
		fill_boundary_domains(part);
		this->Domain_Partitioning::operator=(part);
		return true;
	}

	void change_container(Container &conti){
		cont = &conti;
		pnt.change_container(cont);
		verify_container_solution_correspondence();
	}

	template<class Callable>
	void fill_seed(Callable func, double value_at_infinity = 0.0){
		verify_container_solution_correspondence();

		int indx = 0, dim;
		std::vector<double> coords_comp(n_dim), coords(n_dim);
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			for(int ind = 0; ind < dom->size; ind++){
				cont->update_coordinate(dom, ind, coords_comp.data());
				for(dim = 0; dim < n_dim; ++dim){
					coords[dim] = pnt.compact.xc(dim, coords_comp.data());
					if(std::isinf(coords[dim]))
						break;
				}
				if(dim < n_dim){
					for(int fn = 0; fn < cont->nfunc; fn++)
						cont->data[indx++] = value_at_infinity;
				}
				else{
					for(int fn = 0; fn < cont->nfunc; fn++)
						cont->data[indx++] = func(fn, coords.data());
				}
			}
		}
	}

	template<class Callable>
	Container evaluate_function(Callable &func, int nfunc_fill, double value_at_infinity = 0.0){
		verify_container_solution_correspondence();

		Container result(nfunc_fill, cont->partition);
		for(int dim = 0; dim < n_dim; ++dim){
			if(!cont->parity_st[dim].empty()) result.parity_st[dim] = std::vector<bool>(nfunc_fill, false);
		}
		for(int dim = 0; dim < n_dim; ++dim){
			if(!cont->parity_end[dim].empty()) result.parity_end[dim] = std::vector<bool>(nfunc_fill, false);
		}

		int indx = 0;
		for(const Bulk_Domain *dom = first_bulk_domain; dom != NULL; dom = dom->next_bulk_domain()){
			std::vector<unsigned int> ind(n_dim);

			for(int dindx = 0; dindx < dom->size; dindx++){
				for(int i = n_dim-1, tmpi = dindx; i >= 0; i--){
					if(ind[i] = (tmpi % dom->dim_num(i))) break;
					tmpi /= dom->dim_num(i);
				}
				if(pnt.get_point(dom->bulk_domain_number, ind, true)){
					for(int fn = 0; fn < nfunc_fill; ++fn)
						result.data[indx++] = value_at_infinity;
				}
				else{
					for(int fn = 0; fn < nfunc_fill; ++fn){
						result.data[indx++] = func(fn, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
					}
				}
			}
		}

		return result;
	}
};

#endif