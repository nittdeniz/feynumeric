#ifndef FEYNUMERIC_FOUR_VECTOR_HPP
#define FEYNUMERIC_FOUR_VECTOR_HPP

#include <array>

#include "complex.hpp"
#include "three_vector.hpp"

namespace Feynumeric
{
	class Four_Vector
	{
	private:
		Complex _temporal;
		Three_Vector _spatial;
		bool _contra_variant;
	public:
		Four_Vector(Complex t = Complex{}, Complex x = Complex{}, Complex y = Complex{}, Complex z = Complex{});
		Four_Vector(Four_Vector const& copy);
		Four_Vector& operator=(Four_Vector const& copy);

		Four_Vector& operator+=(Four_Vector const& rhs);
		Four_Vector& operator-=(Four_Vector const& rhs);
		template<typename T>
		Four_Vector& operator*=(T const& rhs);
		template<typename T>
		Four_Vector& operator/=(T const& rhs);

		double theta();
		double phi();
		double spatial_square();
		double square();

		Complex co(std::size_t i) const;
		Complex contra(std::size_t i) const;

		friend Four_Vector operator+(Four_Vector const& lhs, Four_Vector const& rhs);
		friend Four_Vector operator-(Four_Vector const& lhs, Four_Vector const& rhs);
		template<typename T>
		friend Four_Vector operator*(Four_Vector const& lhs, T const& rhs);
		template<typename T>
		friend Four_Vector operator*(T const& lhs, Four_Vector const& rhs);
		template<typename T>
		friend Four_Vector operator/(Four_Vector const& lhs, T const& rhs);

		friend Complex dot(Four_Vector const& lhs, Four_Vector const& rhs);
	};
}
#endif // FEYNUMERIC_FOUR_VECTOR_HPP