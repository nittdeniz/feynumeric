#include "utility.hpp"

namespace Feynumeric
{
    bool almost_identical(double a, double b, double epsilon, double cutoff)
    {
    	if( std::abs(a) < cutoff && std::abs(b) < cutoff ){
    		return true;
    	}
        return std::abs(a-b) < std::abs(epsilon * std::min(a, b));
    }

    bool almost_identical(Complex aa, Complex bb, double epsilon, double cutoff)
    {
        return almost_identical(aa.real(), bb.real(), epsilon, cutoff) && almost_identical(aa.imag(), bb.imag(), epsilon, cutoff);
    }

    Complex dot3(Matrix const& a, Matrix const& b)
    {
        return a.at(0) * b.at(0) + a.at(1) * b.at(1) + a.at(2) * b.at(2);
    }

    Complex dot4(Matrix const& a, Matrix const& b)
    {
        return a.at(0) * b.at(0) - a.at(1) * b.at(1) - a.at(2) * b.at(2) - a.at(3) * b.at(3);
    }

	double int_pow(double c, std::size_t n){
		double result{1};
		for( std::size_t i = 0; i < n; ++i ){
			result *= c;
		}
		return result;
	}

	Complex int_pow(Complex c, std::size_t n){
		Complex result{1};
		for( std::size_t i = 0; i < n; ++i ){
			result *= c;
		}
		return result;
	}

	std::vector<double> lin_space(double from, double to, std::size_t steps){
		double delta = (to-from) / steps;
		std::vector<double> result;
		result.reserve(steps+1);
		for( std::size_t i = 0; i < steps; ++i ){
			result.push_back(from + i * delta);
		}
		result.push_back(to);
		return result;
	}

	std::vector<double> weighted_space(double from, double a, double b, double to, std::size_t N){
		std::vector<double> values1 = lin_space(from, a, N/2-1);
		std::vector<double> values2 = lin_space(a, b, N);
		std::vector<double> values3 = lin_space(b, to, N/2-1);

		std::vector<double> values;
		values.insert(values.end(), values1.begin(), values1.end());
		values.insert(values.end(), values2.begin(), values2.end());
		values.insert(values.end(), values3.begin(), values3.end());
		return values;
	}
}

