

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

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

int main(){
	double const x_min = -1.;
	double const x_max = 1.;

	std::vector<fit_pair> fit_pairs{{0., f_a}, {0.,f_b}, {0., f_c}, {0., f_d}, {0., f_e}};

	std::size_t const n_epochs = 500;
	double rate = 1.2;

	std::vector<double> x_values(20, 0.);
	for( std::size_t i = 0; i < x_values.size(); ++i ){
		x_values[i] = x_min + i * (x_max-x_min) / x_values.size();
	}

	for( std::size_t i = 0; i < n_epochs; ++i ){
		if( i > 0 && i%20 == 0 ){
			rate /= 1.1;
		}
		double loss = 0.;
		std::vector<double> derivatives(fit_pairs.size(), 0.);
		for( std::size_t j = 0; j < x_values.size(); ++j ){
			double const x = x_values[j];
			std::vector<double> sigmoids(fit_pairs.size());
			for( std::size_t k = 0; k < sigmoids.size(); ++k ){
				sigmoids[k] = sigmoid(fit_pairs[k].param);
			}
			double y_hat = 0.;
			for( std::size_t k = 0; k < fit_pairs.size(); ++k ){
				y_hat += sigmoids[k] * fit_pairs[k].func(x);
			}
			double const diff = target(x) - y_hat;
			loss += diff * diff;
			for( std::size_t k = 0; k < fit_pairs.size(); ++k ){
				derivatives[k] += -2. * diff * fit_pairs[k].func(x) * sigmoids[k] * (1-sigmoids[k]);
			}
		}
		for( std::size_t j = 0; j < fit_pairs.size(); ++j ){
			fit_pairs[j].param -= rate * derivatives[j];
		}
		if( loss < 1.e-4 ){
			break;
		}
		std::cout << "ls: " << loss << "\n";
	}
	for( auto [x, y] : fit_pairs ){
		std::cout << "p: " << sigmoid(x) << "\n";
	}
}