#include "feynumeric/format.hpp"
#include "feynumeric/matrix.hpp"
#include "feynumeric/polynomial.hpp"
#include "feynumeric/utility.hpp"

#include <iostream>
#include <sstream>

namespace Feynumeric{
	Polynomial::Polynomial(std::size_t order)
	: n(order+1)
	{
	}

	Polynomial::Polynomial(std::vector<Complex> const& coefficients)
	: _coefficients(coefficients)
	, n(coefficients.size())
	{

	}


	void Polynomial::fit(std::vector<Point> const& data){
		std::vector<Complex> y_data;
		y_data.reserve(n);
		Matrix cramer(n, n, 0);
		for( std::size_t i = 0; i < n; ++i ){
			Complex temp_b{0.};
			for( std::size_t j = i; j < n; ++j ){
				double temp{0.};
				for( std::size_t x = 0; x < data.size(); ++x ){
					temp += std::pow(data[x].x, i+j);
					if( j == i )
					temp_b += std::pow(data[x].x, i) * data[x].y;
				}
				cramer[i * n + j] = temp;
				cramer[j * n + i] = temp;
			}
			y_data.push_back(temp_b);
		}
		auto detM = cramer.det();
		if( almost_identical(detM, 0.) ){
			_coefficients = std::vector<Complex>(n, 0);
			return;
		}
		_coefficients.resize(n);
		for( std::size_t i = 0; i < n; ++i ){
			cramer.swap_col(i, y_data);
			auto temp = cramer.det();
						auto div = temp / detM;
			_coefficients[i] = div;
			cramer.swap_col(i, y_data);
		}
	}

	std::string Polynomial::to_string(char x) const{
		std::stringstream result;
		for( std::size_t i = 0; i < _coefficients.size(); ++i ){
			if( i > 0 ){
				result << "+";
			}
			result << _coefficients[i];
//			result << FORMAT("({:f}", _coefficients[i].real());
//			if( _coefficients[i].imag() < 0 ){
//				result << FORMAT("{:f}I", _coefficients[i].imag());
//			}else{
//				result << FORMAT("+{:f}I", _coefficients[i].imag());
//			}
			result << "*" << x << "^" << i;
		}
		return result.str();
	}

	Complex Polynomial::integrate(double a, double b){
		Complex result{0.};
		for( std::size_t i = 0; i < _coefficients.size(); ++i ){
			result += _coefficients[i]/(i+1.) * (int_pow(b, i+1) - int_pow(a, i+1));
		}
		return result;
	}
}
