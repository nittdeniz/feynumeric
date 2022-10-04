#ifndef FEYNUMERIC_FOUR_VECTOR_HPP
#define FEYNUMERIC_FOUR_VECTOR_HPP

#include <array>

#include "complex.hpp"
#include "lorentz_index.hpp"
#include "three_vector.hpp"

namespace Feynumeric
{
	class Four_Vector
	{
	private:
		Complex _temporal;
		Three_Vector _spatial;
	public:
		Four_Vector(Complex t = Complex{}, Complex x = Complex{}, Complex y = Complex{}, Complex z = Complex{});
		Four_Vector(Complex t, Three_Vector const& r);
		Four_Vector(Four_Vector const& copy);
		Four_Vector& operator=(Four_Vector const& copy);

		Four_Vector& operator+=(Four_Vector const& rhs);
		Four_Vector& operator-=(Four_Vector const& rhs);
		template<typename T>
		Four_Vector& operator*=(T const& rhs);
		template<typename T>
		Four_Vector& operator/=(T const& rhs);

		double theta() const;
		double phi() const;
		double spatial_squared() const;
		double squared() const;

		Three_Vector const& spatial() const;

		Complex t() const;
		Complex E() const;
		Complex x() const;
		Complex y() const;
		Complex z() const;

		Complex co(Lorentz_Index const& mu) const;
		Complex co(Lorentz_Index_Ptr const& mu) const;
		Complex contra(Lorentz_Index const& mu) const;
		Complex contra(Lorentz_Index_Ptr const& mu) const;

		Four_Vector boost(Four_Vector const& p) const;
		Four_Vector align(Four_Vector const& p) const;
		Four_Vector align(Three_Vector const& p) const;

		Four_Vector Rx(double phi) const;
		Four_Vector Ry(double phi) const;
		Four_Vector Rz(double phi) const;

		double cos_theta() const;
		double cos_phi() const;

		friend Four_Vector operator+(Four_Vector const& lhs, Four_Vector const& rhs);
		friend Four_Vector operator-(Four_Vector const& lhs);
		friend Four_Vector operator-(Four_Vector const& lhs, Four_Vector const& rhs);
		template<typename T>
		friend Four_Vector operator*(Four_Vector const& lhs, T const& rhs);
		template<typename T>
		friend Four_Vector operator*(T const& lhs, Four_Vector const& rhs);
		template<typename T>
		friend Four_Vector operator/(Four_Vector const& lhs, T const& rhs);

		friend Complex dot(Four_Vector const& lhs, Four_Vector const& rhs);

        friend std::ostream& operator<<(std::ostream& out, Four_Vector const& v);
	};

	Four_Vector operator+(Four_Vector const& lhs, Four_Vector const& rhs);
	Four_Vector operator-(Four_Vector const& lhs);
	Four_Vector operator-(Four_Vector const& lhs, Four_Vector const& rhs);
	template<typename T>
	Four_Vector operator*(Four_Vector const& lhs, T const& rhs);
	template<typename T>
	Four_Vector operator*(T const& lhs, Four_Vector const& rhs);
	template<typename T>
	Four_Vector operator/(Four_Vector const& lhs, T const& rhs);
    std::ostream& operator<<(std::ostream& out, Four_Vector const& v);

	Complex dot(Four_Vector const& lhs, Four_Vector const& rhs);
}
#endif // FEYNUMERIC_FOUR_VECTOR_HPP