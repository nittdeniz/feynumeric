#ifndef FEYNUMERIC_THREE_VECTOR_HPP
#define FEYNUMERIC_THREE_VECTOR_HPP

#include <array>

#include "complex.hpp"
#include "matrix.hpp"

namespace Feynumeric
{
	class Three_Vector : public Matrix
	{
	public:
		Three_Vector(Complex x = Complex{}, Complex y = Complex{}, Complex z = Complex{});
		Three_Vector(Three_Vector const& copy);
		Three_Vector& operator=(Three_Vector const& copy);

		Three_Vector& operator+=(Three_Vector const& other);
		Three_Vector& operator-=(Three_Vector const& other);

		static Three_Vector from_spherical(double radius, double cos_theta, double cos_phi);

		Complex x() const;
		Complex y() const;
		Complex z() const;

		void x(Complex const& x);
		void y(Complex const& y);
		void z(Complex const& z);

		Three_Vector Rx(double phi) const;
		Three_Vector Ry(double phi) const;
		Three_Vector Rz(double phi) const;

		double squared() const;
		double magnitude() const;
		double theta() const;
		double phi() const;

		double cos_theta() const;
		double cos_phi() const;

		Three_Vector align(Three_Vector const& other) const;

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

	Complex dot(Three_Vector const& lhs, Three_Vector const& rhs);
	Three_Vector operator+(Three_Vector const& lhs, Three_Vector const& rhs);
	Three_Vector operator-(Three_Vector const& lhs, Three_Vector const& rhs);
	template<typename T>
	Three_Vector operator*(Three_Vector const& lhs, T const& rhs)
	{
		return Three_Vector(lhs._data[0] * rhs, lhs._data[1] * rhs, lhs._data[2] * rhs);
	}
	template<typename T>
	Three_Vector operator*(T const& lhs, Three_Vector const& rhs){
		return rhs * lhs;
	}
	template<typename T>
	Three_Vector operator/(Three_Vector const& lhs, T const& rhs){
		auto inverse = 1./rhs;
		return lhs * inverse;
	}
}
#endif // FEYNUMERIC_THREE_VECTOR_HPP