#ifndef FEYNUMERIC_POLYNOMIAL_HPP
#define FEYNUMERIC_POLYNOMIAL_HPP

#include <feynumeric/complex.hpp>

#include <array>
#include <vector>


namespace Feynumeric
{
	class Polynomial
	{
		std::vector<Complex> _coefficients;
		std::size_t n;
	public:
		Polynomial(std::size_t order);


		void fit(std::vector<std::pair<double, Complex>> const& data);
	};
}

#endif