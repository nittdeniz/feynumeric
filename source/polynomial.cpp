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
		//std::cout << cramer << "\n";
//		std::cout << "detM: " << detM << "\n";
		if( almost_identical(detM, 0., 1.e-16, 1.e-16) ){
			_coefficients = std::vector<Complex>(n, 0);
			return;
		}
		_coefficients.resize(n);
		for( std::size_t i = 0; i < n; ++i ){
			cramer.swap_col(i, y_data);
			auto temp = cramer.det();
			auto div = temp / detM;
//			std::cout << "div: " << temp << " /  " << detM << " = " << div << "\n";
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

	Polynomial operator+(Polynomial const& lhs, Polynomial const& rhs){
		Polynomial result(std::max(lhs.n, rhs.n)-1);
		for( std::size_t i = 0; i < result.n; ++i ){
			if( i < lhs.n ){
				result._coefficients[i] += lhs._coefficients.at(i);
			}
			if( i < rhs.n ){
				result._coefficients[i] += rhs._coefficients.at(i);
			}
		}
		return result;
	}

	Polynomial operator-(Polynomial const& lhs, Polynomial const& rhs){
		Polynomial result(std::max(lhs.n, rhs.n)-1);
		for( std::size_t i = 0; i < result.n; ++i ){
			if( i < lhs.n ){
				result._coefficients[i] += lhs._coefficients.at(i);
			}
			if( i < rhs.n ){
				result._coefficients[i] -= rhs._coefficients.at(i);
			}
		}
		return result;
	}

	Polynomial operator*(Polynomial const& lhs, Polynomial const& rhs){
		Polynomial result(lhs.n * rhs.n - 1);
		for( std::size_t i = 0; i < lhs.n; ++i ){
			for( std::size_t j = 0; j < rhs.n; ++j ){
				result._coefficients[i+j] += lhs._coefficients[i] * rhs._coefficients[j];
			}
		}
		return result;
	}

	Polynomial operator*(Polynomial const& lhs, Complex const& rhs){
		Polynomial copy(lhs);
		for( auto& c : copy._coefficients ){
			c *= rhs;
		}
		return copy;
	}

	Polynomial operator*(Complex const& lhs, Polynomial const& rhs){
		return rhs * lhs;
	}

	Polynomial operator/(Polynomial const& lhs, Complex const& rhs){
		return lhs * (1./rhs);
	}

	Polynomial Polynomial::conjugate() const{
		Polynomial result(*this);
		for( auto& coef : result._coefficients ){
			coef = std::conj(coef);
		}
		return result;
	}

	Polynomial::Polynomial(Polynomial const& other)
	: _coefficients(other._coefficients)
	, n(other.n)
	{

	}

	Polynomial& Polynomial::operator=(Polynomial const& other){
		n = other.n;
		_coefficients = other._coefficients;
		return *this;
	}

	Complex Polynomial::operator()(double x) const{
		Complex result{0.};
		for( std::size_t i = 0; i < _coefficients.size(); ++i ){
			result += _coefficients[i] * int_pow(x, i);
		}
		return result;
	}
}
