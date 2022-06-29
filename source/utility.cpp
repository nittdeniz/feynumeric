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
}

