#include <feynumeric/feynumeric.hpp>

#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
#include <effective_lagrangian_model.hpp>

double f_a(double x){ return x*x*x*x;}
double f_b(double x){ return x*x*x;}
double f_c(double x){ return x*x;}
double f_d(double x){ return x;}
double f_e(double x){ return 1.;}

double target(double x){ return 0.69 * std::pow(x, 4) + 0.42 * std::pow(x, 3) + 0.1337 * std::pow(x, 2.) + 0.666 * std::pow(x, 1) + 1; }

double sigmoid(double x){ return 1./(1. + std::exp(-x));}

struct fit_pair{
	double param;
	std::function<double(double)> func;
};

int main(int argc, char** argv){
	using namespace Feynumeric;

	Command_Line_Manager cmd(argc, argv);
	cmd.register_command("particle_file", true, "file with particle parameters");
	cmd.register_command("coupling_constants", true, "file with coupling constants");
	cmd.crash_on_missing_mandatory_command();

	Particle_Manager P(cmd.as_string("particle_file"));
	auto const& Proton = P.get("proton");
	auto const& Pi_Plus = P.get("pi+");
	init_vertices(P, cmd.as_string("coupling_constants"));

	std::vector<Polynomial> s_poly;
	std::vector<Polynomial> c_poly;

	std::vector<std::string> resonances = {"D1232","D1600", "D1620", "D1700", "D1750", "D1900", "D1905", "D1910", "D1920", "D1930", "D1940", "D1950"};

	std::vector<std::vector<FPolynomial<1>>> fpolynomials(resonances.size());
	std::vector<std::vector<FPolynomial<1>>> fpolynomials_conjugated(resonances.size());

	for( std::size_t k = 0; k < resonances.size(); ++k){
		auto particle = P.get(FORMAT("{}pp", resonances[k]));

		/* width */
		std::vector<std::vector<Polynomial>> width_polynomials(1);
		auto n_spin_states = particle->spin().n_states() * Proton->spin().n_states();

		for( std::size_t i = 0; i < n_spin_states; ++i ){
			Polynomial temp;
			temp.load(FORMAT("data/polynomials/polynomial_decay_{}_0_{}.txt", particle->name(), i));
			width_polynomials[0].push_back(temp);
		}

		Amplitude<0> M_width(width_polynomials, {particle}, {Proton, Pi_Plus});

		/* scattering */
		n_spin_states = Proton->spin().n_states() * Proton->spin().n_states();
		std::vector<std::vector<Polynomial>> scattering_polynomials(2);

		for( std::size_t i = 0; i < n_spin_states; ++i ){
			Polynomial temp_s, temp_c;
			temp_s.load(FORMAT("data/polynomials/polynomial_scattering_{}_0_{}.txt", particle->name(), i));
			temp_c.load(FORMAT("data/polynomials/polynomial_scattering_{}_1_{}.txt", particle->name(), i));
			scattering_polynomials[0].push_back(temp_s);
			scattering_polynomials[1].push_back(temp_c);
		}

		func_t<1> form_factor = [=](double sqrt_s) mutable{
			auto const lambda = 0.8;
			auto const l4 = std::pow(lambda, 4);
			auto const delta = sqrt_s * sqrt_s - particle->mass() * particle->mass();
			return l4 / ( delta * delta + l4 );
		};

		func_t<1> breit_wigner = [=](double sqrt_s) mutable {
			return 1. / ( sqrt_s * sqrt_s - particle->mass() * particle->mass() +
			              1.i * sqrt_s * M_width.width(sqrt_s) );
		};

		fpolynomials[k].reserve(scattering_polynomials[0].size());
		for( std::size_t i = 0; i < scattering_polynomials[0].size(); ++i ){
			func_t<1> f = scattering_polynomials[0][i];
			f = f * breit_wigner;
			f = f * form_factor;
			fpolynomials[k].emplace_back(f, scattering_polynomials[1][i]);
			auto f_con = f >> [](Complex x){ return std::conj(x); };
			fpolynomials_conjugated[k].emplace_back(f_con, scattering_polynomials[1][i].conjugate());
		}
	}

	/** t_channel **/
	auto particle = P.get("rho0");

	auto n_spin_states = Proton->spin().n_states() * Proton->spin().n_states();
	std::vector<std::vector<Polynomial>> scattering_polynomials(2);
	for( std::size_t i = 0; i < n_spin_states; ++i ){
		Polynomial temp_s, temp_c;
		temp_s.load(FORMAT("data/polynomials/polynomial_scattering_{}_0_{}.txt", particle->name(), i));
		temp_c.load(FORMAT("data/polynomials/polynomial_scattering_{}_1_{}.txt", particle->name(), i));
		scattering_polynomials[0].push_back(temp_s);
		scattering_polynomials[1].push_back(temp_c);
	}

	func_t<1> form_factor = [=](double sqrt_s) mutable{
		auto const lambda = 0.8;
		auto const l4 = std::pow(lambda, 4);
		auto const delta = sqrt_s * sqrt_s - particle->mass() * particle->mass();
		return l4 / ( delta * delta + l4 );
	};

	func_t<1> breit_wigner = [=](double sqrt_s) mutable {
		return 1. / ( sqrt_s * sqrt_s - particle->mass() * particle->mass());
	};

	for( std::size_t i = 0; i < scattering_polynomials[0].size(); ++i ){
		func_t<1> f = scattering_polynomials[0][i];
		f = f * breit_wigner;
		f = f * form_factor;
		fpolynomials.back().emplace_back(f, scattering_polynomials[1][i]);
		auto f_con = f >> [](Complex x){ return std::conj(x); };
		fpolynomials_conjugated.back().emplace_back(f_con, scattering_polynomials[1][i].conjugate());
	}

	Timer stopwatch;
	stopwatch.start();

	FPolynomial<1> sum;
	for( std::size_t i = 0; i < fpolynomials.size(); ++i ){
		FPolynomial<1> temp = fpolynomials[i][0];
		for( std::size_t j = 1; j < n_spin_states; ++j ){
//			temp = temp + fpolynomials[i][j];
		}
	}
/*
	func_t<1> scattering(){
		using namespace std::placeholders;
		std::function<Complex(Complex)> conjugate = [](Complex c){ return std::conj(c); };

		if( N == 1 ){
			auto qout = std::bind(momentum, _1, _outgoing[0]->mass(), _outgoing[1]->mass());
			auto qin = std::bind(momentum, _1, _incoming[0]->mass(), _incoming[1]->mass());
			func_t<1> phase_space = [=, this](double s) mutable{
				return phase_space2(n_pol, s, qout(s), qin(s));
			};
			func_t<1> sum = [](double){return 0.;};
			for( auto const& f : _fpolynomials ){
				auto fc = f >> conjugate;
				auto full = fc * f;
				auto integrated = full.integrate(-1, 1);
				sum = sum + integrated;
			}
			return phase_space * sum;
		}
		critical_error("scattering only supports 2 final state particles (i.e. N=1)");
	}

	/*

	double const x_min = -1.;
	double const x_max = 1.;

	std::default_random_engine generator;
	generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_real_distribution<double> distribution(-2,2);
	auto r = [&](){return distribution(generator);};

	std::vector<fit_pair> fit_pairs{{r(), f_a}, {r(),f_b}, {r(), f_c}, {r(), f_d}, {r(), f_e}};

	for( auto [x, y] : fit_pairs ){
		std::cout << "p: " << sigmoid(x) << "\n";
	}
	std::cout << "\n";

	std::size_t const n_epochs = 2000;
	double rate = 0.04;

	std::vector<double> x_values(20, 0.);
	for( std::size_t i = 0; i < x_values.size(); ++i ){
		x_values[i] = x_min + i * (x_max-x_min) / x_values.size();
	}

	for( std::size_t i = 0; i < n_epochs; ++i ){
		if( i > 0 && i%200 == 0 ){
			rate *= 1.25;
		}
		double loss = 0.;
		std::vector<double> derivatives(fit_pairs.size(), 0.);
		for( std::size_t j = 0; j < x_values.size(); ++j ){
			double const x = x_values[j];
			std::vector<double> sigmoids(fit_pairs.size());
			std::vector<double> sigmoids_left(fit_pairs.size());
			std::vector<double> sigmoids_right(fit_pairs.size());
			double const epsilon = 0.001;
			for( std::size_t k = 0; k < sigmoids.size(); ++k ){
				sigmoids[k] = sigmoid(fit_pairs[k].param);
				sigmoids_left[k] = sigmoid(fit_pairs[k].param - epsilon);
				sigmoids_right[k] = sigmoid(fit_pairs[k].param + epsilon);
			}
			double y_hat = 0.;
			std::vector<double> y_hat0(fit_pairs.size(), 0.);
			std::vector<double> y_hat1(fit_pairs.size(), 0.);
			for( std::size_t k = 0; k < fit_pairs.size(); ++k ){
				double f = fit_pairs[k].func(x);
//				std::cout << "f: " << f << "\n";
//				std::cout << FORMAT("sigmoids[{}]: {}\n", k, sigmoids[k]);
//				std::cout << "y_hat: += " << sigmoids[k] * f << "\n";
				y_hat += sigmoids[k] * f;
//				std::cout << "y_hat: " << y_hat << "\n";
				for( std::size_t l = 0; l < fit_pairs.size(); ++l ){
					double ff = fit_pairs[l].func(x);
//					std::cout << "ff: " << ff << "\n";
//					std::cout << FORMAT("sigmoids[{}]: {}\n", l, sigmoids[l]);
//					std::cout << FORMAT("sigmoids_L[{}]: {}\n", l, sigmoids_left[l]);
//					std::cout << FORMAT("sigmoids_R[{}]: {}\n", l, sigmoids_right[l]);
//					std::cout << "y_hat: += " << sigmoids[l] * ff << "\n";
					if( k ==l  ){
						y_hat0[k] += sigmoids_left[l] * ff;
						y_hat1[k] += sigmoids_right[l] * ff;
					}else{
						y_hat0[k] += sigmoids[l] * ff;
						y_hat1[k] += sigmoids[l] * ff;
					}
//					std::cout << "y_hat0: " << y_hat0[k] << "\n";
//					std::cout << "y_hat1: " << y_hat1[k] << "\n";
				}

			}
			double const diff = target(x) - y_hat;
			loss += diff * diff;
			for( std::size_t k = 0; k < fit_pairs.size(); ++k ){
//				derivatives[k] += -2. * diff * fit_pairs[k].func(x) * sigmoids[k] * (1-sigmoids[k]);
				double loss0 = std::pow(target(x) - y_hat0[k], 2);
				double loss1 = std::pow(target(x) - y_hat1[k], 2);
//				std::cout << "loss:  " << diff*diff  << "\n";
//				std::cout << "loss0: " << loss0 << "\n";
//				std::cout << "loss1: " << loss1 << "\n";
//				std::cout << "numeric: " << (loss1 - loss0) / (2* epsilon) << "\n";
//				std::cout << "exact: " << -2. * diff * fit_pairs[k].func(x) * sigmoids[k] * (1-sigmoids[k]) << "\n";
				derivatives[k] += (loss1 - loss0) / (2* epsilon);
			}
		}
		for( std::size_t j = 0; j < fit_pairs.size(); ++j ){
			fit_pairs[j].param -= rate * derivatives[j];
		}
		if( loss < 1.e-4 ){
			break;
		}
		if( i > 0 && i%100 == 0 ){
			std::cout << "rate: " << rate << " ls: " << loss << "\n";
		}
	}
	for( auto [x, y] : fit_pairs ){
		std::cout << "p: " << sigmoid(x) << "\n";

	}
	*/
}