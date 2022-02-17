#include "three_vector.hpp"

namespace Feynumeric
{
	Three_Vector::Three_Vector(Complex x, Complex y, Complex z)
	: _data({x, y, z})
	{

	}

	Three_Vector::Three_Vector(Three_Vector const& copy)
	: _data(copy._data)
	{

	}

	Three_Vector& Three_Vector::operator=(Three_Vector const& copy)
	{
		_data = copy._data;
		return *this;
	}

	Three_Vector& Three_Vector::operator+=(Three_Vector const& other)
	{
		_data[0] += other._data[0];
		_data[1] += other._data[1];
		_data[2] += other._data[2];
		return *this;
	}

	Three_Vector& Three_Vector::operator-=(Three_Vector const& other)
	{
		_data[0] -= other._data[0];
		_data[1] -= other._data[1];
		_data[2] -= other._data[2];
		return *this;
	}

	template<typename T>
	Three_Vector& Three_Vector::operator*=(T const& other)
	{
		_data[0] *= other;
		_data[1] *= other;
		_data[2] *= other;
		return *this;
	}

	template<typename T>
	Three_Vector& Three_Vector::operator/=(T const& other)
	{
		_data[0] /= other;
		_data[1] /= other;
		_data[2] /= other;
		return *this;
	}

	Three_Vector Three_Vector::from_spherical(double radius, double cos_theta, double cos_phi)
	{
		auto const sin_theta = std::sqrt(1 - cos_theta * cos_theta);
		auto const sin_phi = std::sqrt(1 - cos_phi * cos_phi);
		return Three_Vector(radius * sin_theta * cos_phi, radius * sin_theta * sin_phi, radius * cos_theta);
	}

	Complex Three_Vector::x() const
	{
		return _data[0];
	}

	Complex Three_Vector::y() const
	{
		return _data[1];
	}

	Complex Three_Vector::z() const
	{
		return _data[2];
	}

	void Three_Vector::x(Complex const& x)
	{
		_data[0] = x;
	}

	void Three_Vector::y(Complex const& y)
	{
		_data[1] = y;
	}

	void Three_Vector::z(Complex const& z)
	{
		_data[2] = z;
	}

	double Three_Vector::square() const
	{
		return (_data[0] * std::conj(_data[0]) + _data[1] * std::conj(_data[1]) + _data[2] * std::conj(_data[2])).real();
	}

	double Three_Vector::theta() const
	{
		auto const& x = _data[0].real();
		auto const& y = _data[1].real();
		auto const& z = _data[2].real();
		double const rho = std::sqrt(x*x + y*y);
		if( rho == 0 && z == 0 )
		{
			return 0.;
		}
		return std::atan2(z, rho);
	}

	double Three_Vector::phi() const
	{
		auto const& x = _data[0].real();
		auto const& y = _data[1].real();
		if( x == 0 && y == 0 )
		{
			return 0.;
		}
		return std::atan2(x, y);
	}

	Complex Three_Vector::at(std::size_t i) const
	{
		return _data[i];
	}

	Complex dot(Three_Vector const& lhs, Three_Vector const& rhs)
	{
		return lhs.x() * rhs.x() + lhs.y() * rhs.y() + lhs.z() * rhs.z();
	}

	Three_Vector operator+(Three_Vector const& lhs, Three_Vector const& rhs)
	{
		Three_Vector copy(lhs);
		copy += rhs;
		return copy;
	}

	Three_Vector operator-(Three_Vector const& lhs, Three_Vector const& rhs)
	{
		Three_Vector copy(lhs);
		copy -= rhs;
		return copy;
	}

	template<typename T>
	Three_Vector operator*(Three_Vector const& lhs, T const& rhs)
	{
		Three_Vector copy(lhs);
		copy *= rhs;
		return copy;
	}

	template<typename T>
	Three_Vector operator*(T const& lhs, Three_Vector const& rhs)
	{
		return rhs * lhs;
	}

	template<typename T>
	Three_Vector operator/(Three_Vector const& lhs, T const& rhs)
	{
		Three_Vector copy(lhs);
		copy /= rhs;
		return copy;
	}
}