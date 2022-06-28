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
		Polynomial(std::vector<Complex> const& coefficients);
		void fit(std::vector<Point> const& data);
		Complex integrate(double a, double b);
		std::string to_string(char x) const;
	};
}

#endif