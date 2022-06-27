#include "feynumeric/format.hpp"
#include "feynumeric/matrix.hpp"
#include "feynumeric/polynomial.hpp"

#include <iostream>
#include <sstream>

namespace Feynumeric{
	Polynomial::Polynomial(std::size_t order)
	: n(order+1)
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
		_coefficients.resize(n);
		for( std::size_t i = 0; i < n; ++i ){
			cramer.swap_col(i, y_data);
			_coefficients[i] = cramer.det() / detM;
			cramer.swap_col(i, y_data);
		}
	}

	std::string Polynomial::to_string(char x) const{
		std::stringstream result;
		for( std::size_t i = 0; i < _coefficients.size(); ++i ){
			result << FORMAT("({:f}", _coefficients[i].real());
			if( _coefficients[i].imag() < 0 ){
				result << FORMAT("{:f}I", _coefficients[i].imag());
			}else{
				result << FORMAT("+{:f}I", _coefficients[i].imag());
			}
			result << ")*" << x << "^" << i << " + ";
		}
		return result.str();
	}
}
