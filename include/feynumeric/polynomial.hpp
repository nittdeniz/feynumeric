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
		Polynomial(Polynomial const& other);
		Polynomial& operator=(Polynomial const& other);
		void fit(std::vector<Point> const& data);
		Complex integrate(double a, double b);
		std::string to_string(char x) const;

		Polynomial conjugate() const;

		Complex operator()(double x) const;

		friend Polynomial operator+(Polynomial const& lhs, Polynomial const& rhs);
		friend Polynomial operator-(Polynomial const& lhs, Polynomial const& rhs);
		friend Polynomial operator*(Polynomial const& lhs, Polynomial const& rhs);
		friend Polynomial operator*(Polynomial const& lhs, Complex const& rhs);
		friend Polynomial operator*(Complex const& lhs, Polynomial const& rhs);
		friend Polynomial operator/(Polynomial const& lhs, Complex const& rhs);
	};
}

#endif