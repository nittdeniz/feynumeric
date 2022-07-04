#ifndef FEYNUMERIC_AMPLITUDE_HPP
#define FEYNUMERIC_AMPLITUDE_HPP

#include "functions.hpp"
#include "polynomial.hpp"

#include <vector>

namespace Feynumeric
{
	template<std::size_t N>
	class Amplitude
	{
	private:
		std::vector<FPolynomial<N>> _fpolynomials;

	public:
		Amplitude(std::vector<std::vector<Polynomial>> const& polynomials);

		void scale(func_t<N> const& func){
			for( auto& fp : _fpolynomials ){
				for( auto& f : fp._coefficients ){
					f = f * func;
				}
			}
		}

		Complex operator()(double x, auto&& args...){
			return _fpolynomials[1](x, args);
		}
	};

	template<>
	Amplitude<0>::Amplitude(std::vector<std::vector<Polynomial>> const& polynomials){
		_fpolynomials.reserve(polynomials[0].size());
		for( std::size_t i = 0; i < polynomials[0].size(); ++i ){
			func_t<0> f = [&](){return Complex(1.);};
			_fpolynomials.emplace_back(f, polynomials[0][i]);
		}
	}

	template<>
	Amplitude<1>::Amplitude(std::vector<std::vector<Polynomial>> const& polynomials){
		_fpolynomials.reserve(polynomials[0].size());
		for( std::size_t i = 0; i < polynomials[0].size(); ++i ){
			func_t<1> f = polynomials[0][i];
			FPolynomial<1> fp(f, polynomials[1][i]);
			_fpolynomials.push_back(fp);
//			_fpolynomials.emplace_back(f, polynomials[1][i]);
		}
	}
}


#endif