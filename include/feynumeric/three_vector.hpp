#ifndef FEYNUMERIC_THREE_VECTOR_HPP
#define FEYNUMERIC_THREE_VECTOR_HPP

#include <array>

#include "complex.hpp"

namespace Feynumeric
{
	class Three_Vector
	{
		std::array<Complex, 3> _data;
	public:
		Three_Vector(Complex x = Complex{}, Complex y = Complex{}, Complex z = Complex{});
		Three_Vector(Three_Vector const& copy);
		Three_Vector& operator=(Three_Vector const& copy);

		Three_Vector& operator+=(Three_Vector const& other);
		Three_Vector& operator-=(Three_Vector const& other);
		template<typename T>
		Three_Vector& operator*=(T const& other);
		template<typename T>
		Three_Vector& operator/=(T const& other);

		static Three_Vector from_spherical(double radius, double cos_theta, double cos_phi);

		Complex x() const;
		Complex y() const;
		Complex z() const;

		void x(Complex const& x);
		void y(Complex const& y);
		void z(Complex const& z);

		Complex at(std::size_t i) const;

		double square() const;
		double theta() const;
		double phi() const;

		friend Complex dot(Three_Vector const& lhs, Three_Vector const& rhs);
		friend Three_Vector operator+(Three_Vector const& lhs, Three_Vector const& rhs);
		friend Three_Vector operator-(Three_Vector const& lhs, Three_Vector const& rhs);
		template<typename T>
		friend Three_Vector operator*(Three_Vector const& lhs, T const& rhs);
		template<typename T>
		friend Three_Vector operator*(T const& lhs, Three_Vector const& rhs);
		template<typename T>
		friend Three_Vector operator/(Three_Vector const& lhs, T const& rhs);
	};
}
#endif // FEYNUMERIC_THREE_VECTOR_HPP