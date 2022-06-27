#ifndef FEYNUMERIC_POLYNOMIAL_HPP
#define FEYNUMERIC_POLYNOMIAL_HPP

#include <feynumeric/complex.hpp>

#include <array>
#include <vector>


namespace Feynumeric
{
	struct Point{
		double x;
		Complex y;
	};
	class Polynomial
	{
		std::vector<Complex> _coefficients;
		std::size_t n;
	public:
		Polynomial(std::size_t order);


		void fit(std::vector<Point> const& data);
		std::string to_string(char x) const;
	};
}

#endif