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
	coordinate_komar_radius(0.5*partition.partition_end(0)), ps_length(0.0)
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
	coordinate_komar_radius(0.5*partition.partition_end(0)), ps_length(0.0)
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
	params(star.params), coordinate_komar_radius(star.coordinate_komar_radius),
	residual_error(star.residual_error), ps_length(star.ps_length) {};

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

// void Boson_Star_Solver::plot_polar_cut(const char *folder_name, int nr, int nz){
// 	struct stat sttus = {0};
// 	if (stat(folder_name, &sttus) == -1) {
// 		mkdir(folder_name, 0700);
// 	}

// 	Evaluator eval(*star);
// 	int n_func = star->number_of_functions();
// 	double rad = 1.5*star->radius_komar;
// 	double dr = 2.0*rad/(1.0*(nr-1));
// 	double dz = 2.0*rad/(1.0*(nz-1));

// 	std::ofstream out_coords;
// 	std::ofstream *out_f = new std::ofstream[n_func];
// 	std::string file_name;

// 	file_name = std::string(folder_name) + "/coords.dat";
// 	out_coords.open(file_name);
// 	for(int fn = 0; fn < n_func; fn++){
// 		file_name = std::string(folder_name) + "/f" + std::to_string(fn) + ".dat";
// 		out_f[fn].open(file_name);
// 		out_f[fn] << std::setprecision(15) << std::fixed;
// 	}

// 	double rc, zc;
// 	std::vector<double> coord(2);

// 	for(int i = 0; i < nr-1; i++)
// 		out_coords << -rad + i*dr << ",";
// 	out_coords << rad << std::endl;
// 	for(int i = 0; i < nz-1; i++)
// 		out_coords << - rad + i*dz << ",";
// 	out_coords << rad << std::endl;

// 	for(int i = 0; i < nz; i++){
// 		zc = - rad + i*dz;
// 		for(int j = 0; j < nr; j++){
// 			rc = -rad + j*dr;
// 			coord[0] = std::sqrt(rc*rc + zc*zc);
// 			coord[1] = std::asin(std::fabs(zc/coord[0]));
// 			coord[0] = xc_inv(0, &(coord[0]));
// 			for(int fn = 0; fn < n_func; fn++){
// 				out_f[fn] << eval(coord, fn);
// 				if(j != nr-1) out_f[fn] << ",";
// 				else out_f[fn] << std::endl;
// 			}
// 		}
// 	}

// 	for(int fn = 0; fn < n_func; fn++)
// 		out_f[fn].close();
// 	delete [] out_f;
// 	out_coords.close();
// }