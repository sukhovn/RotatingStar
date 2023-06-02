#include <collocation_solution.h>
#include <nelder-mead.h>

#ifndef COLLOCATION_ROUTINED_H
#define COLLOCATION_ROUTINED_H

class Collocation_Nelder_Mead: public Nelder_Mead{
	const Evaluator &eval;
	std::vector<double> spec_min_x;
	double spec_min_val;

public:
	const std::vector<double> &min_x;
	const double &min_val;

	Collocation_Nelder_Mead(const Evaluator &evali): Nelder_Mead(evali.n_dim), eval(evali),
		spec_min_x(evali.n_dim, 0.0), spec_min_val(0.0),
		min_x(spec_min_x), min_val(spec_min_val) {};

	double function(const std::vector<double> &x){
		std::vector<double> tmp(x);
		box(tmp, eval.partition_start(), eval.partition_end());
		return eval(tmp, 0);
	}

	void set_simplex(const std::vector<double> &posit, double scale = 0.1){
		x_array[0] = posit;
		for(int i = 0; i < nvar; i++){
			x_array[i+1] = posit;
			x_array[i+1][i] += scale*(eval.partition_end(i) - eval.partition_start(i));
		}

		for(int i = 0; i <= nvar; i++)
			vals[i] = function(x_array[i]);

		sort_points();
	}

	double minimize(bool messages = false){
		int i;

		for(i = 0; i < max_iter; i++){
			nm_step();
			if(largest_side() < max_side) break;
		}

		for(int j = 0; j <= nvar; j++)
			box(x_array[j], eval.partition_start(), eval.partition_end());

		if(i < max_iter && messages){
			std::cout << "Minimization is done\n";
			std::cout << "Number of iterations: " << i << std::endl;
			std::cout << "Minimum is reached at:\n";
			for(int i = 0; i < nvar; i++){
				std::cout << x_array[0][i];
				if(i < nvar-1) std::cout << " ";
			}
			std::cout << std::endl;
			std::cout << "Minimal value reached is: " << vals[0] << std::endl;
		}
		else if(messages){
			std::cout << "Maximal number of iterations reached\n";
			std::cout << "Number of iterations: " << i << std::endl;
			std::cout << "Best minimum reached at:\n";
			for(int i = 0; i < nvar; i++){
				std::cout << x_array[0][i];
				if(i < nvar-1) std::cout << " ";
			}
			std::cout << std::endl;
			std::cout << "Best minimal value reached: " << vals[0] << std::endl;
		}

		spec_min_x = x_array[0];
		spec_min_val = vals[0];
		return vals[0];
	}

	//Diagnostic routine that keeps track of all the steps
	double minimize_track_csv(const std::string &file_name, bool messages = true){
		double res = Nelder_Mead::minimize_track_csv(file_name, messages);
		spec_min_x = x_array[0];
		spec_min_val = vals[0];
		return res;
	}

	//Repeats minimization a few times starting at random points
	double safe_minimization(bool messages = false, int tries = 5, double scale = 0.1){
		srand (time(NULL));
		for(int i = 0; i < tries; i++){
			random_x(eval.partition_start(), eval.partition_end(), scale);
			minimize(messages);
			if(i == 0 || vals[0] < spec_min_val){
				spec_min_val = vals[0];
				spec_min_x = x_array[0];
			}
		}
		return spec_min_val;
	}
};

#endif