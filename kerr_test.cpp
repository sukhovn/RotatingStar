#include <boson_star.h>
#include <supplement.h>

int main(int argc, char const *argv[])
{
	double m = 1.0;
	double a = 0.5;

	Domain_Partitioning part({2.0, 0.0}, {38.0, M_PI_2}, {40, 20});
	part.domain_split({10.0, M_PI_4}, 0);

	std::cout << "Choosing partition:\n" << part << std::endl;

	Boson_Star star(part);
	star.poly_coef.resize(4);
	star.poly_coef = {0.0, 0.0, 0.0, 0.0};
	star.kappa = 0.0;
	star.m_ang = 0.0;
	star.omega = 0.0;
	star.params["omega"] = 0.0;
	star.params["komar_mass"] = m;
	star.params["komar_momentum"] = -m*a;
	star.params["radius_komar"] = 10*m;

	struct{
		double m;
		double a;
		
		double operator()(int fn, double *x){
			double sin_th = std::sin(x[1]), cos_th = std::cos(x[1]);
			double rf = 1.0 - (m*m - a*a)/(4.0*x[0]*x[0]);
			double r2 = 1.0 + m/x[0] + (m*m - a*a)/(4.0*x[0]*x[0]);
			double num = r2*r2 + a*a/(x[0]*x[0])*cos_th*cos_th;
			double denom = std::pow(r2*r2 + a*a/(x[0]*x[0]), 2.0) - a*a/(x[0]*x[0])*rf*rf*sin_th*sin_th;
			switch(fn){
				case 0:
					return std::log(rf);
					break;
				case 1:
					return std::log(rf*num/denom);
					break;
				case 2:
					return 0.5*std::log(num);
					break;
				case 3:
					return 2*m*a/(x[0]*x[0]*x[0])*r2/denom;
					break;
				case 4:
					return 0.0;
					break;
			}
			std::cout << "Wrong equation number" << std::endl;
			throw("Wrong equation number");
		}
	} init;
	init.m = m; init.a = a;
	Radial_Compactification compact(part, true);

	Boson_Star_Solver solver(&star);
	solver.fill_seed(init);

	std::vector<double> equation_error = solver.find_function_maximum(&Rotating_Boson_Star_Equations::equations, 4);
	std::cout << "First equation error maximum is: " << equation_error[0] << std::endl;
	std::cout << "Second equation error maximum is: " << equation_error[1] << std::endl;
	std::cout << "Third equation error maximum is: " << equation_error[2] << std::endl;
	std::cout << "Fourth equation error maximum is: " << equation_error[3] << std::endl;
	std::cout << std::endl;

	star.residual_error = solver.find_function_maximum(&Rotating_Boson_Star_Equations::residuals, 2);
	std::cout << "First residual maximum is: " << star.residual_error[0] << std::endl;
	std::cout << "Second residual maximum is: " << star.residual_error[1] << std::endl;
	std::cout << std::endl;

	solver.find_gh_moments(true);
	
	std::cout << std::endl;
	star.save("../kerr_test");

	return 0;
}