#include <container.h>

#ifndef POINT_CLASS_H
#define POINT_CLASS_H

class Point{
	int n_dim;
	int nfunc;

	std::vector<unsigned int> ind;
	std::vector<std::pair<bool,bool>> doubling;

	const Bulk_Domain *curr_domain;

	void get_derivatives(bool compute_mixed_derivatives){
		int indx, n_i_tmp;
		double ch_coef_1, ch_coef_2, rescale;
		double ch_coef_1_alt, ch_coef_2_alt, tmp1, tmp2;
		//Computing non-mixed derivatives
		for(int i = 0; i < n_dim; i++){
			int shift = 1;
			for(int j = i+1; j < n_dim; j++) shift *= curr_domain->dim_num(j);
			indx = cont->index_st[curr_domain->bulk_domain_number] + (posit - ind[i]*shift)*nfunc;
			shift = (shift-1)*nfunc;
			n_i_tmp = curr_domain->dim_num(i);

			for(int ii = 0; ii < n_i_tmp; ii++){
				if(doubling[i].first){
					if(doubling[i].second){
						tmp1 = cheb_Dn1_el(2*n_i_tmp-1, n_i_tmp-ind[i]-1, n_i_tmp-ii-1);
						tmp2 = cheb_Dn1_el(2*n_i_tmp-1, n_i_tmp-ind[i]-1, n_i_tmp+ii);
						ch_coef_1 = tmp1 + tmp2;
						ch_coef_1_alt = tmp1 - tmp2;

						tmp1 = cheb_Dn2_el(2*n_i_tmp-1, n_i_tmp-ind[i]-1, n_i_tmp-ii-1);
						tmp2 = cheb_Dn2_el(2*n_i_tmp-1, n_i_tmp-ind[i]-1, n_i_tmp+ii);
						ch_coef_2 = tmp1 + tmp2;
						ch_coef_2_alt = tmp1 - tmp2;
					}
					else{
						tmp1 = cheb_Dn1_el(2*n_i_tmp-1, 2*n_i_tmp-ind[i]-1, 2*n_i_tmp-ii-1);
						tmp2 = cheb_Dn1_el(2*n_i_tmp-1, 2*n_i_tmp-ind[i]-1, ii);
						ch_coef_1 = tmp1 + tmp2;
						ch_coef_1_alt = tmp1 - tmp2;

						tmp1 = cheb_Dn2_el(2*n_i_tmp-1, 2*n_i_tmp-ind[i]-1, 2*n_i_tmp-ii-1);
						tmp2 = cheb_Dn2_el(2*n_i_tmp-1, 2*n_i_tmp-ind[i]-1, ii);
						ch_coef_2 = tmp1 + tmp2;
						ch_coef_2_alt = tmp1 - tmp2;
					}
				}
				else{
					ch_coef_1 = cheb_Dn1_el(n_i_tmp-1, n_i_tmp-ind[i]-1, n_i_tmp-ii-1);
					ch_coef_2 = cheb_Dn2_el(n_i_tmp-1, n_i_tmp-ind[i]-1, n_i_tmp-ii-1);
				}
				for(int fn = 0; fn < nfunc; fn++){
					if(doubling[i].first && ((!cont->parity_st[i].empty() && cont->parity_st[i][fn]) || 
											 (!cont->parity_end[i].empty() && cont->parity_end[i][fn])))
					{
						dphi[i][fn] += ch_coef_1_alt * cont->data[indx];
						ddphi[i][fn] += ch_coef_2_alt * cont->data[indx++];
					}
					else{
						dphi[i][fn] += ch_coef_1 * cont->data[indx];
						ddphi[i][fn] += ch_coef_2 * cont->data[indx++];
					}
				}
				indx += shift;
			}

			if(doubling[i].first)
				rescale = 1.0/(curr_domain->domain_end(i) - curr_domain->domain_start(i));
			else
				rescale = 2.0/(curr_domain->domain_end(i) - curr_domain->domain_start(i));
			for(int fn = 0; fn < nfunc; fn++){
				dphi[i][fn] *= rescale;
				ddphi[i][fn] *= rescale*rescale;
			}
		}

		//Computing mixed derivatives
		if(compute_mixed_derivatives){
			int dpos = n_dim;
			double *tmp = new double[nfunc];
			for(int i = 0; i < n_dim; i++){
				n_i_tmp = curr_domain->dim_num(i);
				for(int j = i+1; j < n_dim; j++){
					int shift_j = 1, n_j_tmp = curr_domain->dim_num(j);
					for(int ii = j+1; ii < n_dim; ii++) shift_j *= curr_domain->dim_num(ii);
					int shift_i = shift_j;
					for(int ii = i+1; ii <= j; ii++) shift_i *= curr_domain->dim_num(ii);
					indx = cont->index_st[curr_domain->bulk_domain_number] + (posit - ind[i]*shift_i - ind[j]*shift_j)*nfunc;
					shift_i = (shift_i - curr_domain->dim_num(j)*shift_j)*nfunc;
					shift_j = (shift_j - 1)*nfunc;

					for(int ii = 0; ii < n_i_tmp; ii++){
						for(int fn = 0; fn < nfunc; fn++) tmp[fn] = 0.0;

						if(doubling[i].first){
							if(doubling[i].second){
								tmp1 = cheb_Dn1_el(2*n_i_tmp-1, n_i_tmp-ind[i]-1, n_i_tmp-ii-1);
								tmp2 = cheb_Dn1_el(2*n_i_tmp-1, n_i_tmp-ind[i]-1, n_i_tmp+ii);
								ch_coef_1 = tmp1 + tmp2;
								ch_coef_1_alt = tmp1 - tmp2;
							}
							else{
								tmp1 = cheb_Dn1_el(2*n_i_tmp-1, 2*n_i_tmp-ind[i]-1, 2*n_i_tmp-ii-1);
								tmp2 = cheb_Dn1_el(2*n_i_tmp-1, 2*n_i_tmp-ind[i]-1, ii);
								ch_coef_1 = tmp1 + tmp2;
								ch_coef_1_alt = tmp1 - tmp2;
							}
						}
						else
							ch_coef_1 = cheb_Dn1_el(n_i_tmp-1, n_i_tmp-ind[i]-1, n_i_tmp-ii-1);

						for(int jj = 0; jj < n_j_tmp; jj++){
							if(doubling[j].first){
								if(doubling[j].second){
									tmp1 = cheb_Dn1_el(2*n_j_tmp-1, n_j_tmp-ind[j]-1, n_j_tmp-jj-1);
									tmp2 = cheb_Dn1_el(2*n_j_tmp-1, n_j_tmp-ind[j]-1, n_j_tmp+jj);
									ch_coef_2 = tmp1 + tmp2;
									ch_coef_2_alt = tmp1 - tmp2;
								}
								else{
									tmp1 = cheb_Dn1_el(2*n_j_tmp-1, 2*n_j_tmp-ind[j]-1, 2*n_j_tmp-jj-1);
									tmp2 = cheb_Dn1_el(2*n_j_tmp-1, 2*n_j_tmp-ind[j]-1, jj);
									ch_coef_2 = tmp1 + tmp2;
									ch_coef_2_alt = tmp1 - tmp2;
								}
							}
							else
								ch_coef_2 = cheb_Dn1_el(n_j_tmp-1, n_j_tmp-ind[j]-1, n_j_tmp-jj-1);

							for(int fn = 0; fn < nfunc; fn++){
								if(doubling[i].first && ((!cont->parity_st[i].empty() && cont->parity_st[i][fn]) || 
														 (!cont->parity_end[i].empty() && cont->parity_end[i][fn])))
									tmp[fn] += ch_coef_2_alt * cont->data[indx++];
								else
									tmp[fn] += ch_coef_2 * cont->data[indx++];
							}
							indx += shift_j;
						}
						for(int fn = 0; fn < nfunc; fn++){
							if(doubling[i].first && ((!cont->parity_st[i].empty() && cont->parity_st[i][fn]) || 
													 (!cont->parity_end[i].empty() && cont->parity_end[i][fn])))
								ddphi[dpos][fn] += ch_coef_1_alt * tmp[fn];
							else
								ddphi[dpos][fn] += ch_coef_1 * tmp[fn];
						}
						indx += shift_i;
					}

					rescale = 4.0/(curr_domain->domain_end(i) - curr_domain->domain_start(i))
								/(curr_domain->domain_end(j) - curr_domain->domain_start(j));
					if(doubling[i].first) rescale *= 0.5;
					if(doubling[j].first) rescale *= 0.5;

					for(int fn = 0; fn < nfunc; fn++)
						ddphi[dpos][fn] *= rescale;
					dpos++;
				}
			}
			delete [] tmp;
		}
	}

public:
	const Container *cont;

	int posit;
	double *x;

	double *phi;

	double **dphi;
	double **ddphi;

	Point(const Container *conti): cont(conti), n_dim(conti->n_dim), nfunc(conti->nfunc), doubling(conti->n_dim), ind(conti->n_dim, 0){
		int part_size = n_dim*(n_dim+1)/2;
		phi = new double[nfunc];

		x = new double[n_dim];
		dphi = new double*[n_dim];
		ddphi = new double*[part_size];

		for(int i = 0; i < n_dim; i++)
			dphi[i] = new double[nfunc];
		for(int i = 0; i < part_size; i++)
			ddphi[i] = new double[nfunc];
	}
	~Point(){
		for(int i = 0; i < n_dim*(n_dim+1)/2; i++)
			delete [] ddphi[i];
		for(int i = 0; i < n_dim; i++)
			delete [] dphi[i];
		delete [] phi; delete [] dphi; delete [] ddphi;
		delete [] x;
	}

	void change_container(const Container *conti){
		if(n_dim != conti->n_dim || nfunc != conti->nfunc){
			std::cout << "Number of dimensions or functions cannot be change in a Point class" << std::endl;
			throw("Wrong point class change");
		}
		cont = conti;
	}

	void get_point_fields_only(int dm, const std::vector<unsigned int> &indi){
		ind = indi;
		curr_domain = cont->first_bulk_domain;
		while(curr_domain->bulk_domain_number != dm){
			curr_domain = curr_domain->next_bulk_domain();
			if(curr_domain == NULL){
				std::cout << "Domain index " << dm << " is out of range" << std::endl;
				throw("Domain index is out of range");
			}
		}
		for(int dim = 0; dim < n_dim; ++dim){
			if(!cont->parity_st[dim].empty() && std::fabs(curr_domain->domain_start(dim) - cont->partition_start(dim)) < COMP_TOL)
				doubling[dim] = {true, true};
			else if(!cont->parity_end[dim].empty() && std::fabs(curr_domain->domain_end(dim) - cont->partition_end(dim)) < COMP_TOL)
				doubling[dim] = {true, false};
			else
				doubling[dim] = {false, false};
		}

		for(int i = 0; i < n_dim; i++){
			if(ind[i] >= curr_domain->dim_num(i)){
				std::cout << "Point " << ind[i] << " along the dimension " << i << " is out of range " << curr_domain->dim_num(i)-1 << "\n";
				throw("Point i is out of range");
			}
		}

		for(int i = 0; i < n_dim; i++){
			if(doubling[i].first){
				if(doubling[i].second)
					x[i] = ch_x(2.0*curr_domain->domain_start(i)-curr_domain->domain_end(i), curr_domain->domain_end(i),
								2*curr_domain->dim_num(i)-1, curr_domain->dim_num(i)+ind[i]);
				else
					x[i] = ch_x(curr_domain->domain_start(i), 2.0*curr_domain->domain_end(i)-curr_domain->domain_start(i),
								2*curr_domain->dim_num(i)-1, ind[i]);
			}
			else
				x[i] = ch_x(curr_domain->domain_start(i), curr_domain->domain_end(i), curr_domain->dim_num(i)-1, ind[i]);
		}

		posit = ind[0];
		for(int i = 1; i < n_dim; i++){
			posit = curr_domain->dim_num(i)*posit + ind[i];
		}

		for(int fn = 0; fn < nfunc; fn++){
			phi[fn] = cont->data[cont->index_st[dm] + nfunc*posit + fn];
		}
		for(int i = 0; i < n_dim; i++)
			for(int fn = 0; fn < nfunc; fn++)
				dphi[i][fn] = 0.0;
		for(int i = 0; i < n_dim*(n_dim+1)/2; i++)
			for(int fn = 0; fn < nfunc; fn++)
				ddphi[i][fn] = 0.0;
	}

	void get_point(unsigned int dm, const std::vector<unsigned int> &indi, bool compute_mixed_derivatives){
		get_point_fields_only(dm, indi);
		get_derivatives(compute_mixed_derivatives);
	}

	template<class>
	friend class Point_Compactified;
};

template<class Compactification>
class Point_Compactified: public Point{
public:
	const Compactification &compact;

	double *x_c;
	double *dx;
	double *ddx;

	Point_Compactified(const Compactification &compi, const Container *conti): Point(conti), compact(compi) {
		x_c = new double[n_dim];
		dx = new double[n_dim];
		ddx = new double[n_dim];
	};

	~Point_Compactified(){
		delete [] x_c;
		delete [] dx;
		delete [] ddx;
	}

	void get_point_fields_only(int dm, const std::vector<unsigned int> &indi){
		Point::get_point_fields_only(dm, indi);

		for(int i = 0; i < n_dim; i++){
			x_c[i] = x[i];
			x[i] = compact.xc(i, x_c);
		}
	}

	//Returns 1 if the point lies at compactified infinity
	bool get_point(int dm, const std::vector<unsigned int> &indi, bool compute_mixed_derivatives){
		Point::get_point_fields_only(dm, indi);

		for(int i = 0; i < n_dim; i++){
			x_c[i] = x[i];
			x[i] = compact.xc(i, x_c);
		}
		for(int i = 0; i < n_dim; i++)
			if(std::isinf(x[i])) return true;
		for(int i = 0; i < n_dim; i++){
			dx[i] = compact.dxc(i, x_c);
			ddx[i] = compact.ddxc(i, x_c);
		}

		Point::get_derivatives(compute_mixed_derivatives);

		for(int i = 0; i < n_dim; i++){
			for(int fn = 0; fn < nfunc; fn++){
				dphi[i][fn] /= dx[i];
				ddphi[i][fn] = (ddphi[i][fn] - ddx[i]*dphi[i][fn])/(dx[i]*dx[i]);
			}
		}

		if(compute_mixed_derivatives){
			int dpos = n_dim;
			for(int i = 0; i < n_dim; i++){
				for(int j = i+1; j < n_dim; j++){
					for(int fn = 0; fn < nfunc; fn++)
						ddphi[dpos][fn] /= dx[i]*dx[j];
					++dpos;
				}
			}
		}

		return false;
	}
};

#endif