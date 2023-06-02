#include <boson_star.h>
#include <solution_parameters.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void Boson_Star::read_parameters(const Parameter_File &reader){
	int shift;

	poly_coef.resize(4);

	if(reader.name() == "Solitonic" || reader.name() == "Initial"){
		shift = 2;
		double mass = reader[0];
		double sigma = reader[1];

		if(sigma > 0.0){
			poly_coef[0] = mass*mass;
			poly_coef[1] = -4.0*mass*mass/(sigma*sigma);
			poly_coef[2] = 4.0*mass*mass/std::pow(sigma, 4.0);
		}
		else{
			poly_coef[0] = mass*mass;
		}
	}
	else{
		shift = 4;
		for(int i = 0; i < 5; i++)
			poly_coef[i] = reader[i];
	}

	kappa = reader[shift];
	params["omega"] = omega = reader[shift+1];
	m_ang = reader[shift+2];

	if(reader.name() != "Initial"){
		params["komar_mass"] = reader[shift+3];
		params["komar_momentum"] = reader[shift+4];
		params["radius_charge"] = reader[shift+5];
		params["radius_komar"] = reader[shift+6];
		params["phi_max"] = reader[shift+7];
		params["min_gtt"] = reader[shift+8];
	}
}

void Boson_Star::record_parameters(Parameter_File &reader){
	int shift;

	if(reader.name() == "Solitonic" || reader.name() == "Initial"){
		shift = 2;
		if(std::fabs(poly_coef[1]) < 1.0e-14 && std::fabs(poly_coef[2]) < 1.0e-14){
			reader[0] = std::sqrt(poly_coef[0]);
			reader[1] = 0.0;
		}
		else{
			bool coef_3_cond = (std::fabs(poly_coef[2] - 0.25*poly_coef[1]*poly_coef[1]/poly_coef[0])/poly_coef[2] < 1.0e-14); 
			if(coef_3_cond){
				reader[0] = std::sqrt(poly_coef[0]);
				reader[1] = std::sqrt(-4.0*poly_coef[0]/poly_coef[1]);
			}
			else{
				reader[0] = nan("");
				reader[1] = nan("");
			}
		}
	}
	else{
		shift = 4;
		for(int i = 0; i < 4; i++)
			reader[i] = poly_coef[i];
	}

	reader[shift+0] = kappa;
	reader[shift+1] = params["omega"];
	reader[shift+2] = m_ang;

	reader[shift+3] = params["komar_mass"];
	reader[shift+4] = params["komar_momentum"];
	reader[shift+5] = params["radius_charge"];
	reader[shift+6] = params["radius_komar"];
	reader[shift+7] = params["phi_max"];
	reader[shift+8] = params["min_gtt"];
}

Boson_Star::Boson_Star(const Domain_Partitioning &part):
	Rotating_Boson_Star_Equations(), Container(5, part),
	pfile((const char*) __Star_solution_parameters_dat, __Star_solution_parameters_dat_len),
	params({{"omega", 1.0}, {"komar_mass", 0.0}, {"komar_momentum", 0.0}, {"radius_charge", 0.0}, {"radius_komar", 0.0}, {"phi_max", 0.0}, {"min_gtt", 0.0}}),
	ps_length(0.0)
{
	set_doubling_begin(0, {false, false, false, false, true});
	set_doubling_begin(1, {false, false, false, false, true});
	if(partition.n_dim != 2){
		std::cout << "Wrong number of dimensions in partitioning" << std::endl;
		throw("Wrong number of dimensions in partitioning");
	}
};

Boson_Star::Boson_Star(const std::string &folder_name):
	Rotating_Boson_Star_Equations(), Container(5, get_partitioning_from_data(folder_name)),
	pfile((const char*) __Star_solution_parameters_dat, __Star_solution_parameters_dat_len),
	params({{"omega", 1.0}, {"komar_mass", 0.0}, {"komar_momentum", 0.0}, {"radius_charge", 0.0}, {"radius_komar", 0.0}, {"phi_max", 0.0}, {"min_gtt", 0.0}}),
	ps_length(0.0)
{
	if(partition.n_dim != 2){
		std::cout << "Wrong number of dimensions in partitioning" << std::endl;
		throw("Wrong number of dimensions in partitioning");
	}
	load(folder_name);
};

Boson_Star::Boson_Star(const Boson_Star &star):
	Rotating_Boson_Star_Equations(star), Container(star),
	pfile((const char*) __Star_solution_parameters_dat, __Star_solution_parameters_dat_len),
	params(star.params), residual_error(star.residual_error), ps_length(star.ps_length) {};

Boson_Star &Boson_Star::operator=(const Boson_Star &star){
	this->Rotating_Boson_Star_Equations::operator=(star);
	this->Container::operator=(star);
	params = star.params;
	residual_error = star.residual_error;
	ps_length = star.ps_length;
	return *this;
}

Boson_Star gaussian_seed_initial_configuration(const std::string &parameter_file, const Domain_Partitioning &part, bool verbose){
	#include <init_setup_parameters.h>
	std::string str((const char*) __Star_init_setup_parameters_dat, __Star_init_setup_parameters_dat_len);
	Parameter_File pfile(str);

	if(pfile.name() != "Initial"){
		std::cout << "Incompatible setup parameters in initial\n";
		throw("Incompatible setup parameters in initial");
	}

	pfile.read_file(parameter_file);
	
	Boson_Star star(part);
	star.read_parameters(pfile);
	
	struct{
		double sx;
		double sz;
		double f0;
		
		double operator()(int fn, double *x){
			double rho = x[0]*sin(x[1]);
			double z = x[0]*cos(x[1]);

			switch(fn){
				case 0:
					return 0.0;
					break;
				case 1:
					return 0.0;
					break;
				case 2:
					return 0.0;
					break;
				case 3:
					return 0.0;
					break;
				case 4:
					return f0 * rho * rho * exp(-rho*rho/(sx*sx)) * exp(-z*z/(sz*sz));
					break;
			}
			std::cout << "Wrong equation number" << std::endl;
			throw("Wrong equation number");
		}
	} init;
	init.sx = pfile[5]; init.sz = pfile[6]; init.f0 = pfile[7];

	if(verbose){
		std::cout << "Initial potential parameters:" << std::endl;
		std::cout << "Mass: " << pfile[0] << std::endl;
		std::cout << "Sigma: " << pfile[1] << std::endl;
		std::cout << std::endl;

		std::cout << "Initial star parameters:" << std::endl;
		std::cout << "Frequency w: " << pfile[3] << std::endl;
		std::cout << "Angular harmonic number m: " << pfile[4] << std::endl;
		std::cout << std::endl;

		std::cout << "Initial partition:" << std::endl;
		print_partitioning(part);
	}
	Radial_Compactification compact(part, verbose);
	if(verbose){
		std::cout << "\nGaussian seed parameters:" << std::endl;
		std::cout << "sx value: " << pfile[5] << std::endl;
		std::cout << "sz value: " << pfile[6] << std::endl;
		std::cout << "f0 value: " << pfile[7] << std::endl;
		std::cout << std::endl;
	}

	Collocation_Solution<Radial_Compactification> coll_sol(compact, &star);
	coll_sol.fill_seed(init);

	return star;
}

void Boson_Star_Solver::find_star_parameters(void){
	star->params["omega"] = star->omega;
	std::vector<double> approx_min;

	auto factor_function = [](int fn, double *x, double *phi, double **dphi, double **ddphi)
	{
		return std::exp(0.5*(phi[0] - phi[1]));
	};

	Evaluator metric_prefactor(evaluate_function(factor_function, 1, {1.0}));

	//Calculating density functions
	auto density_function = [this](int fn, double *x, double *phi, double **dphi, double **ddphi)
	{
		if(x[0] == 0.0 || x[1] == 0) return 0.0;
		double factor = 1.0;
		if(fn == 1) factor = star->m_ang;
		if(fn == 2) factor = std::exp(2.0 * phi[2] + 0.5 * (phi[0] - phi[1]));
		return star->kappa*factor*star->densities(fn, x, phi, dphi, ddphi)*star->area_element(x)*dxc_factor(xc_inv(0, x));
	};

	Container densities = evaluate_function(density_function, 3, {0.0, 0.0, 0.0});
	densities.set_doubling_begin(0, {false, false, false});
	densities.set_doubling_begin(1, {true, true, true});
	Evaluator densities_eval(densities);

	//Calculating mass, momentum and energy
	std::vector<double> charges = densities_eval.integrate();
	star->params["komar_mass"] = charges[0];
	star->params["komar_momentum"] = charges[1];
	star->params["adm_energy"] = charges[2];

	//Finding partial integral of densities to find radii
	int int_size = density_extrapolation_resolution;
	double xst = partition_start(0), xend = partition_end(0);
	double *xrange = new double[int_size];
	ch_x_array(xst, xend, int_size-1, xrange);
	double *density_integrals = new double[3*int_size];
	densities_eval.integrate_keep_first_dimension(int_size, xrange, density_integrals);

	//Performing FFT to find spectra
	fftw_r2r_kind kind = FFTW_REDFT00;
	fftw_plan planCOS = fftw_plan_many_r2r(1, &int_size, 3,
								density_integrals, NULL, 1, int_size,
								density_integrals, NULL, 1, int_size, &kind, FFTW_ESTIMATE);
	fftw_execute_r2r(planCOS, density_integrals, density_integrals);
	fftw_destroy_plan(planCOS);

	for(int i = 0; i < 3*int_size; ++i)
		density_integrals[i] = density_integrals[i]/(1.0*(int_size-1));
	for(int i = 0; i < 3; ++i){
		double *spect = density_integrals + int_size*i;
		cheb_int(xst, xend, spect, int_size-1, spect);
		spect[0] -= 2.0*chebev(xst, xend, spect, int_size, xst);
	}

	//Evaluating Komar radius
	double *spec = density_integrals, rad;
	secant_solve([&xst, &xend, &spec, &int_size, &charges](double x){return chebev(xst, xend, spec, int_size, x) - RAD_FRACTION*charges[0];}, xst, xend, &rad);
	star->params["radius_komar"] = metric_prefactor({rad, M_PI_2}, 0)*xc(0, &rad);

	//Evaluating charge radius
	spec = density_integrals+int_size;
	secant_solve([&xst, &xend, &spec, &int_size, &charges](double x){return chebev(xst, xend, spec, int_size, x) - RAD_FRACTION*charges[1];}, xst, xend, &rad);
	star->params["radius_charge"] = metric_prefactor({rad, M_PI_2}, 0)*xc(0, &rad);

	//Evaluating energy radius
	spec = density_integrals+2*int_size;
	secant_solve([&xst, &xend, &spec, &int_size, &charges](double x){return chebev(xst, xend, spec, int_size, x) - RAD_FRACTION*charges[2];}, xst, xend, &rad);
	star->params["radius_adm"] = metric_prefactor({rad, M_PI_2}, 0)*xc(0, &rad);

	delete [] xrange;
	delete [] density_integrals;

	//Printing out masses and radii
	std::cout << "Komar mass is: " << star->params["komar_mass"] << std::endl;
	std::cout << "Komar angular momentum is: " << star->params["komar_momentum"] << std::endl;
	std::cout << "ADM energy is: " << star->params["adm_energy"] << std::endl;

	std::cout << "Komar radius estimate is: " << star->params["radius_komar"] << std::endl;
	std::cout << "Charge radius estimate is: " << star->params["radius_charge"] << std::endl;
	std::cout << "ADM energy radius estimate is: " << star->params["radius_adm"] << std::endl;

	//Calculating maximal scalar field value
	struct{
		double xmax;
		double ymax;
		double phi_max = 0.0;

		double operator()(int fn, double *x, double *phi, double **dphi, double **ddphi){
			if(phi[4] > phi_max){
				xmax = x[0];
				ymax = x[1];
				phi_max = phi[4];
			}
			return -phi[4];
		}
	} neg_phi;

	Container sc_field = evaluate_function(neg_phi, 1, {0.0});
	sc_field.set_doubling_begin(0, {true});
	sc_field.set_doubling_begin(1, {true});
	Evaluator sc_field_eval(sc_field);
	Collocation_Nelder_Mead sc_field_minimizer(sc_field_eval);
	sc_field_minimizer.set_simplex({neg_phi.xmax, neg_phi.ymax});
	star->params["phi_max"] = -sc_field_minimizer.minimize();

	std::cout << "Maximal value of the scalar field is: " << star->params["phi_max"] << std::endl;

	//Calculating minimal g_tt value
	struct{
		double xmin;
		double ymin;
		double gtt_min = 1.0;

		double operator()(int fn, double *x, double *phi, double **dphi, double **ddphi){
			double tmp = x[0]*std::sin(x[1])*phi[3];
			tmp = std::exp(phi[0]+phi[1]) - std::exp(phi[0]-phi[1])*tmp*tmp;
			if(tmp < gtt_min){
				xmin = x[0];
				ymin = x[1];
				gtt_min = tmp;
			}
			return tmp;
		}
	} g_tt_function;

	Container g_tt = evaluate_function(g_tt_function, 1, {1.0});
	g_tt.set_doubling_begin(0, {false});
	g_tt.set_doubling_begin(1, {false});
	Evaluator g_tt_eval(g_tt);
	Collocation_Nelder_Mead g_tt_minimizer(g_tt_eval);
	g_tt_minimizer.set_simplex({g_tt_function.xmin, g_tt_function.ymin});
	star->params["min_gtt"] = g_tt_minimizer.minimize();
	double gtt_min_x = g_tt_minimizer.min_x[0];

	std::cout << "Minimal value of g_tt component " << star->params["min_gtt"] << " reached at r = " << xc(0, &gtt_min_x) << std::endl;
	if(star->params["min_gtt"] < 0)
		std::cout << "The star has an ergoregion" << std::endl;
	else
		std::cout << "The star does not have an ergoregion" << std::endl;
}

void Boson_Star_Solver::plot_polar_cut(const std::string &folder_name, int nr, int nz){
	struct stat sttus = {0};
	if(stat(folder_name.c_str(), &sttus) == -1){
		mkdir(folder_name.c_str(), 0700);
	}

	struct{
		double operator()(int fn, double *x, double *phi, double **dphi, double **ddphi){
			if(fn == 0)
				return phi[4];
			double tmp = x[0]*std::sin(x[1])*phi[3];
			return std::exp(phi[0]+phi[1]) - std::exp(phi[0]-phi[1])*tmp*tmp;
		}
	} functions_plot;

	std::ofstream out_coords(folder_name + "/coords.dat");
	int nfplot = 1;
	std::vector<std::ofstream> out_f(1);
	out_f[0].open(folder_name + "/field.dat");
	std::vector<double> inf_val = {0.0};
	if(star->params["min_gtt"] < 0){
		nfplot = 2;
		inf_val = {0.0, 1.0};
		out_f.push_back(std::ofstream(folder_name + "/g_tt.dat"));
	}

	Container plot_functions = evaluate_function(functions_plot, nfplot, inf_val);
	if(nfplot == 1){
		plot_functions.set_doubling_begin(0, {true});
		plot_functions.set_doubling_begin(1, {true});
	}
	else{
		plot_functions.set_doubling_begin(0, {true, false});
		plot_functions.set_doubling_begin(1, {true, false});
	}
	// plot_functions.save(folder_name);
	Evaluator eval(plot_functions);

	double rad = 1.5*star->params["radius_komar"];
	double dr = 2.0*rad/(1.0*(nr-1));
	double dz = 2.0*rad/(1.0*(nz-1));

	for(std::ofstream &ofs : out_f)
		ofs << std::setprecision(15) << std::fixed;

	double rc, zc;
	std::vector<double> coord(2);

	for(int i = 0; i < nr-1; i++)
		out_coords << -rad + i*dr << ",";
	out_coords << rad << std::endl;
	for(int i = 0; i < nz-1; i++)
		out_coords << - rad + i*dz << ",";
	out_coords << rad << std::endl;

	for(int i = 0; i < nz; i++){
		zc = - rad + i*dz;
		for(int j = 0; j < nr; j++){
			rc = -rad + j*dr;
			coord[0] = std::sqrt(rc*rc + zc*zc);
			coord[1] = std::asin(std::fabs(zc/coord[0]));
			coord[0] = xc_inv(0, &(coord[0]));
			for(int fn = 0; fn < nfplot; fn++){
				out_f[fn] << eval(coord, fn);
				if(j != nr-1) out_f[fn] << ",";
				else out_f[fn] << std::endl;
			}
		}
	}
}