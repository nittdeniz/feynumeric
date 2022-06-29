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
	class Polynomial2D;

	class Polynomial
	{
		std::size_t n;
		std::size_t order;
		std::vector<Complex> _coefficients;
	public:
		Polynomial(std::size_t order = 0);
		Polynomial(std::vector<Complex> const& coefficients);
		Polynomial(Polynomial const& other);
		void fit(std::vector<Point> const& data);
		Complex integrate(double a, double b);
		std::string to_string(char x) const;

		Polynomial conjugate() const;

		Complex operator()(double x) const;

		Polynomial& operator=(Polynomial const& other);
		Polynomial& operator+=(Polynomial const& other);
		Polynomial& operator-=(Polynomial const& other);
		Polynomial& operator*=(Polynomial const& other);
		Polynomial& operator*=(Complex const& scale);
		Polynomial& operator/=(Complex const& scale);

		friend Polynomial operator+(Polynomial const& lhs, Polynomial const& rhs);
		friend Polynomial operator-(Polynomial const& lhs, Polynomial const& rhs);
		friend Polynomial operator*(Polynomial const& lhs, Polynomial const& rhs);
		friend Polynomial operator*(Polynomial const& lhs, Complex const& rhs);
		friend Polynomial operator*(Complex const& lhs, Polynomial const& rhs);
		friend Polynomial operator/(Polynomial const& lhs, Complex const& rhs);

		friend class Polynomial2D;
	};

	class Polynomial2D{
		std::array<std::vector<Polynomial>, 2> _coefficients;
		std::array<Polynomial, 2> _polynomials;

		Polynomial2D(Polynomial const& first, Polynomial const& second);
	};
}

#endif