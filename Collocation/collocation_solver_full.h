#include <collocation_solution.h>
#include <ios>

#ifndef COLLOCATION_SOLVER_FULL_H
#define COLLOCATION_SOLVER_FULL_H

//This functions find maximum and minimum values out of a and b.
#ifndef FMAX
#define FMAX(a, b) (a > b ? a : b)
#endif

//Difference with respect to which the equations are numerically differentiated
#define H_DERIV 1.0e-7

//Collocation iteration parameters
#define ALPHA 1.0e-4 //Lambda check parameter
#define TOLX 1.0e-10 //Minimum step parameter

#define STPMX 100.0 //Maximal Newton step

#define COLLOC_THRSHLD 1.0e-12
#define MAX_ITER 100 //Limit for maximal number of iterations

template<class Compactification, class Equations_Class>
class Collocation_Solver_Full: public Collocation_Solution<Compactification>{
public:
	Collocation_Matrix colm;
	double eqn_error;
	int iterations;
	const Equations_Class *eq_system;

	using Domain_Partitioning::n_dim;
	using Domain_Partitioning::first_domain;
	using Domain_Partitioning::first_bulk_domain;
	using Domain_Partitioning::partition_start;
	using Domain_Partitioning::partition_end;

private:
	using Collocation_Solution<Compactification>::cont;
	using Collocation_Solution<Compactification>::pnt;
	using Collocation_Solution<Compactification>::fill_boundary_domains;
	using Collocation_Solution<Compactification>::verify_container_solution_correspondence;

	typedef double (Equations_Class::*EqnPointer)(int, double *, double *, double **, double **) const;
	typedef double (Equations_Class::*BcPointer)(int, int, double *, double *) const;

	double deq_d_phi(const EqnPointer eqn_func, int fn_i, int fn_j){
			double tmp;
			pnt.phi[fn_j] += H_DERIV;
			tmp = (eq_system->*eqn_func)(fn_i, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
			pnt.phi[fn_j] -= 2.0*H_DERIV;
			tmp -= (eq_system->*eqn_func)(fn_i, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
			pnt.phi[fn_j] += H_DERIV;
			return tmp/(2.0*H_DERIV);
	}
	
	double deq_d_dphi(const EqnPointer eqn_func, int fn_i, int dim_j, int fn_j){
			double tmp;
			pnt.dphi[dim_j][fn_j] += H_DERIV;
			tmp = (eq_system->*eqn_func)(fn_i, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
			pnt.dphi[dim_j][fn_j] -= 2.0*H_DERIV;
			tmp -= (eq_system->*eqn_func)(fn_i, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
			pnt.dphi[dim_j][fn_j] += H_DERIV;
			return tmp/(2.0*H_DERIV);
	}

	double deq_d_ddphi(const EqnPointer eqn_func, int fn_i, int dim_j, int fn_j){
			double tmp;
			pnt.ddphi[dim_j][fn_j] += H_DERIV;
			tmp = (eq_system->*eqn_func)(fn_i, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
			pnt.ddphi[dim_j][fn_j] -= 2.0*H_DERIV;
			tmp -= (eq_system->*eqn_func)(fn_i, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
			pnt.ddphi[dim_j][fn_j] += H_DERIV;
			return tmp/(2.0*H_DERIV);
	}

	double dbc_d_phi(const BcPointer bc_func, int dim, int fn_i, int fn_j){
			double tmp;
			pnt.phi[fn_j] += H_DERIV;
			tmp = (eq_system->*bc_func)(dim, fn_i, pnt.phi, pnt.dphi[dim]);
			pnt.phi[fn_j] -= 2.0*H_DERIV;
			tmp -= (eq_system->*bc_func)(dim, fn_i, pnt.phi, pnt.dphi[dim]);
			pnt.phi[fn_j] += H_DERIV;
			return tmp/(2.0*H_DERIV);
	}

	double dbc_d_dphi(const BcPointer bc_func, int dim, int fn_i, int fn_j){
			double tmp;
			pnt.dphi[dim][fn_j] += H_DERIV;
			tmp = (eq_system->*bc_func)(dim, fn_i, pnt.phi, pnt.dphi[dim]);
			pnt.dphi[dim][fn_j] -= 2.0*H_DERIV;
			tmp -= (eq_system->*bc_func)(dim, fn_i, pnt.phi, pnt.dphi[dim]);
			pnt.dphi[dim][fn_j] += H_DERIV;
			return tmp/(2.0*H_DERIV);
	}

	template<bool FILL_RHS, bool FILL_JACOBIAN>
	void fill_left_boundary_condition(const Bulk_Domain *bulk, const Boundary_Domain *boundary){
		double tmp;
		int dm = bulk->bulk_domain_number;
		int dim = boundary->boundary_dimension;
		std::vector<unsigned int> ind(n_dim);

		int iter_size = 1;
		std::vector<int> shift(boundary->right_num_st()), margins(n_dim-1);
		for(int i = 0; i < n_dim-1; ++i){
			margins[i] = boundary->right_num_end(i)-shift[i]+1;
			iter_size *= margins[i];
		}

		ind[dim] = 0;
		for(int indx = 0; indx < iter_size; indx++){
			for(int i = n_dim-1, tmpi = indx, itr = n_dim-2; i >= 0; --i){
				if(i == dim) continue;
				if((ind[i] = tmpi % margins[itr] + shift[itr]) > shift[itr]) break;
				tmpi /= margins[itr];
				--itr;
			}

			pnt.get_point(dm, ind, false);
			for(int fn = 0; fn < cont->nfunc; fn++){
				if(FILL_RHS) colm.right_hand_side(dm, ind, fn) = eq_system->boundary_condition_left(dim, fn, pnt.phi, pnt.dphi[dim]);

				if(FILL_JACOBIAN){
					for(int fnd = 0; fnd < cont->nfunc; fnd++){
						tmp = dbc_d_phi(&Equations_Class::boundary_condition_left, dim, fn, fnd);
						if(tmp != 0) colm.row_add_Dn_factor(dm, ind, fn, -1, -1, dm, ind, fnd, tmp);
						tmp = dbc_d_dphi(&Equations_Class::boundary_condition_left, dim, fn, fnd);
						if(tmp != 0) colm.row_add_Dn_factor(dm, ind, fn, -1, dim, dm, ind, fnd, tmp);
					}
				}
			}
		}
	}

	template<bool FILL_RHS, bool FILL_JACOBIAN>
	void fill_right_boundary_condition(const Bulk_Domain *bulk, const Boundary_Domain *boundary){
		double tmp;
		int dm = bulk->bulk_domain_number;
		int dim = boundary->boundary_dimension;
		std::vector<unsigned int> ind(n_dim);

		int iter_size = 1;
		std::vector<int> shift(boundary->left_num_st()), margins(n_dim-1);
		for(int i = 0; i < n_dim-1; ++i){
			margins[i] = boundary->left_num_end(i)-shift[i]+1;
			iter_size *= margins[i];
		}

		ind[dim] = bulk->dim_num(dim)-1;
		for(int indx = 0; indx < iter_size; indx++){
			for(int i = n_dim-1, tmpi = indx, itr = n_dim-2; i >= 0; --i){
				if(i == dim) continue;
				if((ind[i] = tmpi % margins[itr] + shift[itr]) > shift[itr]) break;
				tmpi /= margins[itr];
				--itr;
			}

			pnt.get_point(dm, ind, false);
			for(int fn = 0; fn < cont->nfunc; fn++){
				if(FILL_RHS) colm.right_hand_side(dm, ind, fn) = eq_system->boundary_condition_right(dim, fn, pnt.phi, pnt.dphi[dim]);

				if(FILL_JACOBIAN){
					for(int fnd = 0; fnd < cont->nfunc; fnd++){
						tmp = dbc_d_phi(&Equations_Class::boundary_condition_right, dim, fn, fnd);
						if(tmp != 0) colm.row_add_Dn_factor(dm, ind, fn, -1, -1, dm, ind, fnd, tmp);
						tmp = dbc_d_dphi(&Equations_Class::boundary_condition_right, dim, fn, fnd);
						if(tmp != 0) colm.row_add_Dn_factor(dm, ind, fn, -1, dim, dm, ind, fnd, tmp);
					}
				}
			}
		}
	}

	template<bool FILL_RHS, bool FILL_JACOBIAN>
	void fill_interdomain_boundary(const Bulk_Domain *bulk, const Boundary_Domain *boundary){
		int dm = bulk->bulk_domain_number;
		int dim = boundary->boundary_dimension;
		int indx, n_end, n_opp_range, dm_opp;
		std::vector<unsigned int> ind(n_dim), ind2(n_dim), margins(n_dim-1);
		std::vector<double> coef(n_dim-1);
		std::vector<double> coef_alt(n_dim-1);
		double coef_tmp, coef_tmp2;
		double *val_tmp, *dom_val_copy;
		const Bulk_Domain *bulk_opp;
		std::vector<std::pair<bool,bool>> doubling_opp(n_dim);

		if(boundary->left_side_domain() == bulk){
			bulk_opp = boundary->right_side_domain();
			dm_opp = bulk_opp->bulk_domain_number;
		}
		else{
			bulk_opp = boundary->left_side_domain();
			dm_opp = bulk_opp->bulk_domain_number;
		}

		for(int i = 0; i < n_dim; ++i){
			if(!cont->parity_st[i].empty() && std::fabs(bulk_opp->domain_start(i) - cont->partition_start(i)) < COMP_TOL)
				doubling_opp[i] = {true, true};
			else if(!cont->parity_end[i].empty() && std::fabs(bulk_opp->domain_end(i) - cont->partition_end(i)) < COMP_TOL)
				doubling_opp[i] = {true, false};
			else
				doubling_opp[i] = {false, false};
		}

		if(boundary->left_side_domain() == bulk){
			int iter_size = 1, iter_opp_size = 1;
			for(int i = 0; i < n_dim-1; i++) iter_size *= (margins[i] = (boundary->left_num_end(i) - boundary->left_num_st(i) + 1));
			for(int i = 0; i < n_dim; i++){if(i != dim) iter_opp_size *= bulk_opp->dim_num(i);}

			ind2[dim] = 0;

			if(FILL_RHS){
				val_tmp = new double[cont->nfunc];
				dom_val_copy = new double[cont->nfunc*iter_opp_size];
				if(FILL_RHS)
				for(int indx = 0, itr = 0; indx < iter_opp_size; indx++){
					for(int i = n_dim-1, tmpi = indx; i >= 0; i--){
						if(i == dim) continue;
						if((ind2[i] = tmpi % bulk_opp->dim_num(i))) break;
						tmpi /= bulk_opp->dim_num(i);
					}

					pnt.get_point(dm_opp, ind2, false);
					for(int fn = 0; fn < cont->nfunc; fn++)
						dom_val_copy[itr++] = pnt.phi[fn];
				}
			}

			ind[dim] = bulk->dim_num(dim)-1;
			for(int indx = 0; indx < iter_size; indx++){
				for(int i = n_dim-1, ii = n_dim-2, tmpi = indx; i >= 0; i--){
					if(i == dim) continue;
					bool quit = (ind[i] = tmpi % margins[ii]);
					ind[i] += boundary->left_num_st(ii);
					if(quit) break;
					tmpi /= margins[ii--];
				}

				pnt.get_point(dm, ind, false);
				if(FILL_RHS) for(int fn = 0; fn < cont->nfunc; fn++) val_tmp[fn] = 0.0;

				for(int indx2 = 0, itr = 0; indx2 < iter_opp_size; indx2++){
					for(int i = n_dim-1, ii = n_dim-2, tmpi = indx2; i >= 0; i--){
						if(i == dim) continue;
						bool quit = (ind2[i] = tmpi % bulk_opp->dim_num(i));
						if(doubling_opp[i].first){
							if(doubling_opp[i].second){
								coef_tmp = cheb_Ex_el(2*bulk_opp->dim_num(i)-1, bulk_opp->dim_num(i)+ind2[i], 2*bulk_opp->domain_start(i)-bulk_opp->domain_end(i), bulk_opp->domain_end(i), pnt.x_c[i]);
								coef_tmp2 = cheb_Ex_el(2*bulk_opp->dim_num(i)-1, bulk_opp->dim_num(i)-ind2[i]-1, 2*bulk_opp->domain_start(i)-bulk_opp->domain_end(i), bulk_opp->domain_end(i), pnt.x_c[i]);
							}
							else{
								coef_tmp = cheb_Ex_el(2*bulk_opp->dim_num(i)-1, ind2[i], bulk_opp->domain_start(i), 2*bulk_opp->domain_end(i)-bulk_opp->domain_start(i), pnt.x_c[i]);
								coef_tmp2 = cheb_Ex_el(2*bulk_opp->dim_num(i)-1, 2*bulk_opp->dim_num(i)-ind2[i]-1, bulk_opp->domain_start(i), 2*bulk_opp->domain_end(i)-bulk_opp->domain_start(i), pnt.x_c[i]);
							}
							coef[ii] = coef_tmp + coef_tmp2;
							coef_alt[ii--] = coef_tmp - coef_tmp2;
						}
						else
							coef[ii--] = cheb_Ex_el(bulk_opp->dim_num(i)-1, ind2[i], bulk_opp->domain_start(i), bulk_opp->domain_end(i), pnt.x_c[i]);
						if(quit) break;
						tmpi /= bulk_opp->dim_num(i);
					}

					for(int fn = 0; fn < cont->nfunc; fn++){
						coef_tmp = 1.0;
						for(int i = 0, ii = 0; i < n_dim; i++){
							if(i == dim) continue;
							if(doubling_opp[i].first && ((!cont->parity_st[i].empty() && cont->parity_st[i][fn]) || (!cont->parity_end[i].empty() && cont->parity_end[i][fn])))
								coef_tmp *= coef_alt[ii];
							else
								coef_tmp *= coef[ii];
							++ii;
						}
						if(FILL_JACOBIAN) colm.row_add_Dn_factor(dm, ind, fn, -1, -1, dm_opp, ind2, fn, -coef_tmp);
						if(FILL_RHS) val_tmp[fn] += coef_tmp*dom_val_copy[itr++];
					}
				}

				for(int fn = 0; fn < cont->nfunc; fn++){
					if(FILL_RHS) colm.right_hand_side(dm, ind, fn) = pnt.phi[fn] - val_tmp[fn];
					if(FILL_JACOBIAN) colm.row_add_Dn_factor(dm, ind, fn, -1, -1, dm, ind, fn, 1.0);
				}
			}
		}
		else{
			int iter_size = 1, iter_opp_size = 1;
			for(int i = 0; i < n_dim-1; i++) iter_size *= (margins[i] = (boundary->right_num_end(i) - boundary->right_num_st(i) + 1));
			for(int i = 0; i < n_dim; i++){if(i != dim) iter_opp_size *= bulk_opp->dim_num(i);}

			ind2[dim] = bulk_opp->dim_num(dim)-1;

			if(FILL_RHS){
				val_tmp = new double[cont->nfunc];
				dom_val_copy = new double[cont->nfunc*iter_opp_size];
				for(int indx = 0, itr = 0; indx < iter_opp_size; indx++){
					for(int i = n_dim-1, tmpi = indx; i >= 0; i--){
						if(i == dim) continue;
						if((ind2[i] = tmpi % bulk_opp->dim_num(i))) break;
						tmpi /= bulk_opp->dim_num(i);
					}

					pnt.get_point(dm_opp, ind2, false);
					for(int fn = 0; fn < cont->nfunc; fn++)
						dom_val_copy[itr++] = pnt.dphi[dim][fn];
				}
			}

			ind[dim] = 0;
			for(int indx = 0; indx < iter_size; indx++){
				for(int i = n_dim-1, ii = n_dim-2, tmpi = indx; i >= 0; i--){
					if(i == dim) continue;
					bool quit = (ind[i] = tmpi % margins[ii]);
					ind[i] += boundary->right_num_st(ii);
					if(quit) break;
					tmpi /= margins[ii--];
				}

				pnt.get_point(dm, ind, false);

				if(FILL_RHS) for(int fn = 0; fn < cont->nfunc; fn++) val_tmp[fn] = 0.0;

				for(int indx2 = 0, itr = 0; indx2 < iter_opp_size; indx2++){
					for(int i = n_dim-1, ii = n_dim-2, tmpi = indx2; i >= 0; i--){
						if(i == dim) continue;
						bool quit = (ind2[i] = tmpi % bulk_opp->dim_num(i));
						if(doubling_opp[i].first){
							if(doubling_opp[i].second){
								coef_tmp = cheb_Ex_el(2*bulk_opp->dim_num(i)-1, bulk_opp->dim_num(i)+ind2[i], 2*bulk_opp->domain_start(i)-bulk_opp->domain_end(i), bulk_opp->domain_end(i), pnt.x_c[i]);
								coef_tmp2 = cheb_Ex_el(2*bulk_opp->dim_num(i)-1, bulk_opp->dim_num(i)-ind2[i]-1, 2*bulk_opp->domain_start(i)-bulk_opp->domain_end(i), bulk_opp->domain_end(i), pnt.x_c[i]);
							}
							else{
								coef_tmp = cheb_Ex_el(2*bulk_opp->dim_num(i)-1, ind2[i], bulk_opp->domain_start(i), 2*bulk_opp->domain_end(i)-bulk_opp->domain_start(i), pnt.x_c[i]);
								coef_tmp2 = cheb_Ex_el(2*bulk_opp->dim_num(i)-1, 2*bulk_opp->dim_num(i)-ind2[i]-1, bulk_opp->domain_start(i), 2*bulk_opp->domain_end(i)-bulk_opp->domain_start(i), pnt.x_c[i]);
							}
							coef[ii] = coef_tmp + coef_tmp2;
							coef_alt[ii--] = coef_tmp - coef_tmp2;
						}
						else
							coef[ii--] = cheb_Ex_el(bulk_opp->dim_num(i)-1, ind2[i], bulk_opp->domain_start(i), bulk_opp->domain_end(i), pnt.x_c[i]);
						if(quit) break;
						tmpi /= bulk_opp->dim_num(i);
					}

					for(int fn = 0; fn < cont->nfunc; fn++){
						coef_tmp = 1.0;
						for(int i = 0, ii = 0; i < n_dim; i++){
							if(i == dim) continue;
							if(doubling_opp[i].first && ((!cont->parity_st[i].empty() && cont->parity_st[i][fn]) || (!cont->parity_end[i].empty() && cont->parity_end[i][fn])))
								coef_tmp *= coef_alt[ii];
							else
								coef_tmp *= coef[ii];
							++ii;
						}
						if(FILL_JACOBIAN) colm.row_add_Dn_factor(dm, ind, fn, -1, dim, dm_opp, ind2, fn, -coef_tmp);
						if(FILL_RHS) val_tmp[fn] += coef_tmp*dom_val_copy[itr++];
					}
				}

				for(int fn = 0; fn < cont->nfunc; fn++){
					if(FILL_RHS) colm.right_hand_side(dm, ind, fn) = pnt.dphi[dim][fn] - val_tmp[fn];
					if(FILL_JACOBIAN) colm.row_add_Dn_factor(dm, ind, fn, -1, dim, dm, ind, fn, 1.0);
				}
			}
		}

		if(FILL_RHS){
			delete [] val_tmp;
			delete [] dom_val_copy;
		}
	}

	template<bool FILL_RHS, bool FILL_JACOBIAN>
	void fill_bulk(const Bulk_Domain *bulk){
		double tmp, deq, ddeq;
		int dm = bulk->bulk_domain_number;
		std::vector<unsigned int> ind(n_dim);

		int iter_size = 1;
		std::vector<int> shift(n_dim, 1), margins(bulk->dim_num());
		for(int i = 0; i < n_dim; i++){
			if(!cont->parity_st[i].empty() && std::fabs(bulk->domain_start(i) - cont->partition_start(i)) < COMP_TOL){
				shift[i] = 0;
				iter_size *= --margins[i];
			}
			else if(!cont->parity_end[i].empty() && std::fabs(bulk->domain_end(i) -  cont->partition_end(i)) < COMP_TOL)
				iter_size *= --margins[i];
			else
				iter_size *= (margins[i] -= 2);
		}
		for(int indx = 0; indx < iter_size; indx++){
			for(int i = n_dim-1, tmpi = indx; i >= 0; i--){
				if((ind[i] = (tmpi % margins[i]) + shift[i]) > shift[i]) break;
				tmpi /= margins[i];
			}

			pnt.get_point(dm, ind, false);
			for(int fn = 0; fn < cont->nfunc; fn++){
				if(FILL_RHS) colm.right_hand_side(dm, ind, fn) = eq_system->equations(fn, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi);
				
				if(FILL_JACOBIAN){
					//Adding field derivatives to the Jacobian
					for(int fnd = 0; fnd < cont->nfunc; fnd++){
						tmp = deq_d_phi(&Equations_Class::equations, fn, fnd);
						colm.row_add_Dn_factor(dm, ind, fn, -1, -1, dm, ind, fnd, tmp);
					}

					//Adding field derivatives' derivatives to the Jacobian
					for(int dim = 0; dim < n_dim; dim++){
						for(int fnd = 0; fnd < cont->nfunc; fnd++){
							deq = deq_d_dphi(&Equations_Class::equations, fn, dim, fnd);
							ddeq = deq_d_ddphi(&Equations_Class::equations, fn, dim, fnd);

							if(deq != 0.0 || (ddeq != 0.0 && pnt.ddx[dim] != 0)){
								tmp = deq/pnt.dx[dim] - ddeq*pnt.ddx[dim]/std::pow(pnt.dx[dim], 3.0);
								colm.row_add_Dn_factor(dm, ind, fn, -1, dim, dm, ind, fnd, tmp);
							}
							if(ddeq != 0) colm.row_add_Dn_factor(dm, ind, fn, dim, dim, dm, ind, fnd, ddeq/(pnt.dx[dim]*pnt.dx[dim]));
						}
					}
				}
			}
		}
	}

public:
	Collocation_Solver_Full(const Compactification &compact, const Equations_Class *eq_sys, Container *conti):
		Collocation_Solution<Compactification> (compact, conti), eq_system(eq_sys), colm(*conti)
	{
		fill_boundary_domains(*this);
	};

	bool verify_container_solver_correspondence(void){
		if(verify_container_solution_correspondence()){
			colm.regrid(*cont);
			return true;
		}
		else return false;
	}

	void change_container(Container &conti){
		cont = &conti;
		pnt.change_container(cont);
		verify_container_solver_correspondence();
	}

	template<bool FILL_RHS, bool FILL_JACOBIAN>
	void collocation_fill(){
		verify_container_solver_correspondence();
		const Boundary_Domain *boundary;

		if(FILL_RHS) colm.setZero_rhs();
		if(FILL_JACOBIAN) colm.setZero_jacobian();

		for(const Bulk_Domain *bulk = first_bulk_domain; bulk != NULL; bulk = bulk->next_bulk_domain()){
			for(const Generic_Domain *iterator = first_domain; iterator != NULL; iterator = iterator->next_domain()){
				boundary = dynamic_cast<const Boundary_Domain *>(iterator);
				if(boundary != NULL){
					try{
						if(boundary->left_side_domain() == NULL && boundary->right_side_domain() == bulk)
							fill_left_boundary_condition<FILL_RHS, FILL_JACOBIAN>(bulk, boundary);
						else if(boundary->right_side_domain() == NULL && boundary->left_side_domain() == bulk)
							fill_right_boundary_condition<FILL_RHS, FILL_JACOBIAN>(bulk, boundary);
						else if(boundary->right_side_domain() == bulk || boundary->left_side_domain() == bulk)
							fill_interdomain_boundary<FILL_RHS, FILL_JACOBIAN>(bulk, boundary);
					}
					catch(const char *err){
						std::cout << "\nA problem occured with a boundary domain" << std::endl;
						boundary->print_domain(std::cout);
						throw(err);
					}
				}
				else if(iterator == bulk){
					try{
						fill_bulk<FILL_RHS, FILL_JACOBIAN>(bulk);
					}
					catch(const char *err){
						std::cout << "\nA problem occured with a bulk domain " << bulk->bulk_domain_number << std::endl;
						bulk->print_domain(std::cout);
						throw(err);
					}
				}
			}
		}
	}

	std::vector<double> find_function_maximum(EqnPointer func, int nfunc = 1, bool interior_check_only = true){
		verify_container_solver_correspondence();
		double tmp;
		std::vector<double> max_diff(nfunc, 0.0);
		std::vector<unsigned int> ind(n_dim);

		for(const Bulk_Domain *bulk = first_bulk_domain; bulk != NULL; bulk = bulk->next_bulk_domain()){
			int iter_size = 1;
			std::vector<int> shift(n_dim, 0), margins(bulk->dim_num());
			if(interior_check_only){
				for(int i = 0; i < n_dim; i++){
					if(!cont->parity_st[i].empty() && std::fabs(bulk->domain_start(i) - cont->partition_start(i)) < COMP_TOL)
						iter_size *= --margins[i];
					else if(!cont->parity_end[i].empty() && std::fabs(bulk->domain_end(i) -  cont->partition_end(i)) < COMP_TOL){
						shift[i] = 1;
						iter_size *= --margins[i];
					}
					else{
						shift[i] = 1;
						iter_size *= (margins[i] -= 2);
					}
				}
			}
			else
				iter_size = std::accumulate(margins.begin(), margins.end(), 1, std::multiplies<double>());

			for(int indx = 0; indx < iter_size; indx++){
				for(int i = n_dim-1, tmpi = indx; i >= 0; i--){
					if((ind[i] = (tmpi % margins[i]) + shift[i]) > shift[i]) break;
					tmpi /= margins[i];
				}

				pnt.get_point(bulk->bulk_domain_number, ind, true);
				for(int fn = 0; fn < nfunc; fn++){
					tmp = std::fabs((eq_system->*func)(fn, pnt.x, pnt.phi, pnt.dphi, pnt.ddphi));
					if(tmp > max_diff[fn]) max_diff[fn] = tmp;
				}
			}
		}
		return max_diff;
	}

	Domain_Partitioning get_partitioning(void){
		verify_container_solver_correspondence();
		Domain_Partitioning part(*this);
		if(part.boundary_domains_present()) part.clear_boundary_domains();
		return part;
	}

private:
	int lnsearch(double stepmax){
		double temp = colm.delta.norm(); //Limiting maximal step 
		if(temp > stepmax)
			colm.delta *= stepmax/temp;

		double test = 0.0; //Adjusting minimal lambda
		int itr = 0;
		std::vector<unsigned int> ind(n_dim);
		for(const Bulk_Domain *bulk = first_bulk_domain; bulk != NULL; bulk = bulk->next_bulk_domain()){
			for(int indx = 0; indx < bulk->size; indx++){
				for(int i = n_dim-1, tmpi = indx; i >= 0; i--){
					if(ind[i] = (tmpi % bulk->dim_num(i))) break;
					tmpi /= bulk->dim_num(i);
				}
				for(int fn = 0; fn < cont->nfunc; fn++){
					temp = std::fabs(colm.delta(colm.index(bulk->bulk_domain_number, ind, fn)))/FMAX(std::fabs(cont->data[itr++]), 1.0);
					if(temp > test) test = temp;
				}
			}
		}

		double alamin = TOLX/test;
		double alam1, alam2, alamc = 1.0;
		double f1, f2, fc;
		double f = eqn_error;
		double a, b, disc;

		double diff_to_origin = 1.0;
		cont->subtract_delta(colm);

		collocation_fill<true, false>();
		fc = 0.5*colm.rhs.squaredNorm();

		while(1){
			if(fc <= (1.0 - 2.0*ALPHA) * f){ //Checking if the step with current lambda is adequatly small
				eqn_error = fc;
				return 0;
			}
	//		std::cout << "Continuing Lambda search\n";

			alam2 = alam1;
			f2 = f1;
			alam1 = alamc;
			f1 = fc;

			if(alamc == 1.0){ //Initial lambda step
				alamc = f/(f + fc);
			}
			else{ //Secondary lambda step
				a = ((f1 + (2.0 * alam1 - 1.0) * f)/(alam1 * alam1) - (f2 + (2.0 * alam2 - 1.0) * f)/(alam2 * alam2))/(alam1-alam2);
				b = (alam1 * (f2 + (2.0 * alam2 - 1.0) * f)/(alam2 * alam2) - alam2 * (f1 + (2.0 * alam1 - 1.0) * f)/(alam1 * alam1))/(alam1-alam2);
				if(a == 0.0)
					alamc = f/b;
				else{
					disc = b * b + 6.0 * a * f;
					if(disc < 0.0)
						alamc = 0.5 * alam1;
					else if(b <= 0.0)
						alamc = (-b + std::sqrt(disc))/(3.0 * a);
					else
						alamc = 2.0 * f/(b + std::sqrt(disc));
				}
				if(alamc > 0.5 * alam1)
					alamc = 0.5 * alam1;
			}
			alamc = FMAX(alamc, 0.1 * alam1);

			if(alamc < alamin){
				std::cout << "Reached minimal lambda " << alamc << '\n';
				cont->subtract_delta(colm, -diff_to_origin);
				collocation_fill<true, false>();
				return 1;
			}

			cont->subtract_delta(colm, alamc-diff_to_origin);
			diff_to_origin = alamc;

			collocation_fill<true, false>();
			fc = 0.5*colm.rhs.squaredNorm();
		}
	}

	int iteration(double stepmax){
		int code;//._.

		collocation_fill<false, true>();

		colm.compute_delta();

		// cont->subtract_delta(colm);
		// collocation_fill<true, false>();
		// eqn_error = 0.5*colm.rhs.squaredNorm();

		// std::cout << "Searching for the next step\n";

		code = lnsearch(stepmax);

		iterations++;

		return 0;//code;
	}

public:
	int solve(void){
		verify_container_solver_correspondence();
		int code;
		double err, stepmax;
		std::ios oldState(nullptr);
		oldState.copyfmt(std::cout);

		collocation_fill<true, false>();
		eqn_error = 0.5*colm.rhs.squaredNorm();

		stepmax = STPMX*FMAX(std::sqrt(cont->norm()), (double) cont->size());
		iterations = 0;

		for(int i = 0; i < MAX_ITER; i++){
			code = iteration(stepmax);
			err = std::sqrt(2.0*eqn_error)/(1.0*colm.size());
			std::cout << "Step is done, the error is: " << std::scientific << err << std::endl;
			if(err < COLLOC_THRSHLD){
				std::cout << "Calculation is done, the error is: " << std::scientific << err << std::endl;
				std::cout << "Number of iterations is: " << iterations << "\n" << std::endl;
				std::cout.copyfmt(oldState);
				return 0;
			}
			if(code != 0){
				std::cout << "Reached the minimal newton step, the error is: " << std::scientific << err << std::endl;
				std::cout << "Number of iterations is: " << iterations << "\n" << std::endl;
				std::cout.copyfmt(oldState);
				return 1;
			}
		}

		std::cout << "Solution did not converge, no data is recorded\n" << std::endl;
		std::cout.copyfmt(oldState);
		return 1;
	}
};

#endif