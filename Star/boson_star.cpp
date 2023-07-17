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
		params["M2"] = reader[shift+9];
		params["S3"] = reader[shift+10];
		// params["M4"] = reader[shift+11];
		// params["S5"] = reader[shift+12];
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
	reader[shift+9] = params["M2"];
	reader[shift+10] = params["S3"];
	// reader[shift+11] = params["M4"];
	// reader[shift+12] = params["S5"];
}

Boson_Star::Boson_Star(const Domain_Partitioning &part):
	Rotating_Boson_Star_Equations(), Container(5, part),
	pfile((const char*) __Star_solution_parameters_dat, __Star_solution_parameters_dat_len),
	params({{"omega", 1.0}, {"komar_mass", 0.0}, {"komar_momentum", 0.0}, {"radius_charge", 0.0}, {"radius_komar", 0.0}, {"phi_max", 0.0}, {"min_gtt", 0.0},
			{"M2", 0.0}, {"S3", 0.0}}),
	ps_length(0.0)
{
	if(std::fabs(part.partition_start(0) - 0.0) < 1.0e-10)
		set_doubling_begin(0, {false, false, false, false, true});
	if(std::fabs(part.partition_start(1) - 0.0) < 1.0e-10)
		set_doubling_begin(1, {false, false, false, false, true});
	if(partition.n_dim != 2){
		std::cout << "Wrong number of dimensions in partitioning" << std::endl;
		throw("Wrong number of dimensions in partitioning");
	}
};

Boson_Star::Boson_Star(const std::string &folder_name):
	Rotating_Boson_Star_Equations(), Container(5, get_partitioning_from_data(folder_name)),
	pfile((const char*) __Star_solution_parameters_dat, __Star_solution_parameters_dat_len),
	params({{"omega", 1.0}, {"komar_mass", 0.0}, {"komar_momentum", 0.0}, {"radius_charge", 0.0}, {"radius_komar", 0.0}, {"phi_max", 0.0}, {"min_gtt", 0.0},
			{"M2", 0.0}, {"S3", 0.0}}),
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

void Boson_Star_Solver::find_gh_moments(bool find_sub, const std::string &save_folder){
	double xend = 10.0*star->params["radius_komar"];
	double xst = 0.9*xend;

	auto multipole_function = [this](int fn, double *x, double *phi, double **dphi, double **ddphi)
	{
		if(x[0] == 0.0 || x[1] == 0) return 0.0;
		switch(fn){
			case 0:
				return star->b_mode(0, x, phi);
				break;
			case 1:
				return star->rho_mode(2, x, phi);
				break;
			case 2:
				return star->omega_mode(2, x, phi);
				break;
			case 3:
				return star->b_mode(2, x, phi);
				break;
			case 4:
				return star->rho_mode(4, x, phi);
				break;
			case 5:
				return star->omega_mode(4, x, phi);
				break;
			default:
				std::cout << "Wrong function number" << std::endl;
				throw("wrong function number");
				break;
		}
	};

	auto estimate_zero = [](int n, double *x, double *y){
		double xs = 0.0, ys = 0.0, x2 = 0.0, y2 = 0.0, xy = 0.0;
		for(int i = 0; i < n; ++i){
			xs += x[i];
			ys += y[i];
			x2 += x[i]*x[i];
			y2 += y[i]*y[i];
			xy += x[i]*y[i];
		}
		double beta_est = (ys*x2 - xs*xy)/(n*x2 - xs*xs);
		double sigma_est = (n*xy*xy - 2*xs*ys*xy + x2*ys*ys + y2*xs*xs - n*x2*y2)/(n*xs*xs - n*n*x2);
		double beta_err = sigma_est*x2/(n*x2 - xs*xs);
		double err = std::sqrt(std::fabs(beta_err/(beta_est*beta_est)));
		return std::make_pair(beta_est, err);
	};

	int n_asymptotics = 3;
	if(find_sub){
		n_asymptotics = 6;
	}

	Container multipoles = evaluate_function(multipole_function, n_asymptotics);
	if(find_sub){
		multipoles.set_doubling_begin(0, {false, false, false, false, false, false});
		multipoles.set_doubling_begin(1, {false, true, true, false, true, true});
	}
	else{
		multipoles.set_doubling_begin(0, {false, false, false});
		multipoles.set_doubling_begin(1, {false, true, true});
	}
	Evaluator multipoles_eval(multipoles);

	int int_size = 11;
	double *xrange = new double[int_size];
	for(int i = 0; i < int_size; ++i){
		double tmp = xst + (xend - xst)/(1.0*(int_size-1))*i;
		xrange[i] = xc_inv(0, &tmp);
	}
	double *density_integrals = new double[n_asymptotics*int_size];
	multipoles_eval.integrate_keep_first_dimension(int_size, xrange, density_integrals);

	//Rescaling functions
	for(int i = 0; i < int_size; ++i){
		xrange[i] = xst + (xend - xst)/(1.0*(int_size-1))*i;
		density_integrals[i] *= std::pow(xrange[i], 2);
		density_integrals[int_size+i] *= std::pow(xrange[i], 3);
		density_integrals[2*int_size+i] *= std::pow(xrange[i], 5);
		if(find_sub){
			density_integrals[3*int_size+i] *= std::pow(xrange[i], 4);
			density_integrals[4*int_size+i] *= std::pow(xrange[i], 5);
			density_integrals[5*int_size+i] *= std::pow(xrange[i], 7);
		}
	}

	if(!save_folder.empty()){
		std::ofstream file(save_folder + "/moment_functions.csv");
		if(find_sub)
			file << "r,b0,rho2,omega2,b2,rho4,omega4" << std::endl;
		else
			file << "r,b0,rho2,omega2" << std::endl;
		for(int i = 0; i < int_size; ++i){
			file << xrange[i];
			for(int fn = 0; fn < n_asymptotics; ++fn)
				file << ',' << density_integrals[fn*int_size+i];
			file << std::endl;
		}
	}

	for(int i = 0; i < int_size; ++i) xrange[i] = 1.0/xrange[i];

	double M = star->params["komar_mass"];
	double J = -star->params["komar_momentum"];
	std::pair<double, double> asym_coef;
	double B0, rho2, omega2;

	asym_coef = estimate_zero(int_size, xrange, density_integrals);
	if(asym_coef.second > 0.001){
		std::cout << "B0 coefficient estimation error is " << asym_coef.second << std::endl;
		std::cout << "Couldn't estimate B0 coefficient\nCouldn't estimate leading GH moments" << std::endl;
		return;
	}
	else B0 = asym_coef.first;

	asym_coef = estimate_zero(int_size, xrange, density_integrals+int_size);
	if(asym_coef.second > 0.001){
		std::cout << "rho2 coefficient estimation error is " << asym_coef.second << std::endl;
		std::cout << "Couldn't estimate rho2 coefficient\nCouldn't estimate leading GH moments" << std::endl;
		return;
	}
	else rho2 = asym_coef.first;

	asym_coef = estimate_zero(int_size, xrange, density_integrals+2*int_size);
	if(asym_coef.second > 0.001){
		std::cout << "omega2 coefficient estimation error is " << asym_coef.second << std::endl;
		std::cout << "Couldn't estimate omega2 coefficient\nCouldn't estimate leading GH moments" << std::endl;
		return;
	}
	else omega2 = asym_coef.first;

	double M2 = 0.5*rho2 + M*(M*M + 4.0*B0)/3.0;
	double S3 = -0.75*omega2 + 3.0*J*(M*M + 4.0*B0)/5.0;
	M2 /= J*J/M;
	S3 /= J*J*J/(M*M);

	star->params["M2"] = M2;
	star->params["S3"] = S3;

	std::cout << "Reduced multipole moment M2: " << star->params["M2"] << std::endl;
	std::cout << "Reduced multipole moment S3: " << star->params["S3"] << std::endl;

	if(find_sub){
		double B2, rho4, omega4;

		asym_coef = estimate_zero(int_size, xrange, density_integrals+3*int_size);
		if(std::fabs(asym_coef.first) > 1.0e-10 && asym_coef.second > 0.001){
			std::cout << "B2 coefficient estimation error is " << asym_coef.second << std::endl;
			std::cout << "Couldn't estimate B2 coefficient\nCouldn't estimate subleading GH moments" << std::endl;
			star->params["M4"] = 0.0;
			star->params["S5"] = 0.0;
			return;
		}
		else B2 = asym_coef.first;

		asym_coef = estimate_zero(int_size, xrange, density_integrals+4*int_size);
		if(asym_coef.second > 0.001){
			std::cout << "rho4 coefficient estimation error is " << asym_coef.second << std::endl;
			std::cout << "Couldn't estimate rho4 coefficient\nCouldn't estimate subleading GH moments" << std::endl;
			star->params["M4"] = 0.0;
			star->params["S5"] = 0.0;
			return;
		}
		else rho4 = asym_coef.first;

		asym_coef = estimate_zero(int_size, xrange, density_integrals+5*int_size);
		if(asym_coef.second > 0.001){
			std::cout << "omega4 coefficient estimation error is " << asym_coef.second << std::endl;
			std::cout << "Couldn't estimate omega4 coefficient\nCouldn't estimate subleading GH moments" << std::endl;
			star->params["M4"] = 0.0;
			star->params["S5"] = 0.0;
			return;
		}
		else omega4 = asym_coef.first;

		double M4 = 16.0*B0*B0*M/5.0 - 8.0*B2*M/7.0 - 6.0*J*J*M/35.0 + 32.0*B0*M*M*M/21.0 + 19.0*std::pow(M, 5)/105.0 + 12.0*B0*rho2/7.0 + 4.0*M*M*rho2/7.0 - rho4/2.0;
		double S5 = 48.0*B0*B0*J/7.0 - 40.0*B2*J/21.0 - 2.0*J*J*J/7.0 + 208.0/63.0*B0*J*M*M + 25.0*J*std::pow(M, 4)/63.0 + 5.0*J*M*rho2/21.0 - 10.0*B0*omega2/3.0 - 5.0*M*M*omega2/6.0 + 5.0*omega4/6.0;

		M4 /= J*J*J*J/(M*M*M);
		S5 /= J*J*J*J*J/(M*M*M*M);

		star->params["M4"] = M4;
		star->params["S5"] = S5;

		std::cout << "Reduced multipole moment M4: " << star->params["M4"] << std::endl;
		std::cout << "Reduced multipole moment S5: " << star->params["S5"] << std::endl;
	}

	delete [] xrange;
	delete [] density_integrals;
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