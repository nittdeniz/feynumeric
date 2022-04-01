#include "three_vector.hpp"

namespace Feynumeric
{
	Three_Vector::Three_Vector(Complex x, Complex y, Complex z)
	: Matrix(3, 1, {x, y, z})
	{

	}

	Three_Vector::Three_Vector(Three_Vector const& copy)
	: Matrix(copy)
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

	Three_Vector Three_Vector::Rx(double phi) const
	{
		auto const c = std::cos(phi);
		auto const s = std::sqrt(1 - c*c);
		auto const& x = _data[0];
		auto const& y = _data[1];
		auto const& z = _data[2];
		return Three_Vector(x, y * c - z * s, z * c + y * s);
	}

	Three_Vector Three_Vector::Ry(double phi) const
	{
		auto const c = std::cos(phi);
		auto const s = std::sqrt(1 - c*c);
		auto const& x = _data[0];
		auto const& y = _data[1];
		auto const& z = _data[2];
		return Three_Vector(x * c + z * s, y, z * c - x * s);
	}

	Three_Vector Three_Vector::Rz(double phi) const
	{
		auto const c = std::cos(phi);
		auto const s = std::sqrt(1 - c*c);
		auto const& x = _data[0];
		auto const& y = _data[1];
		auto const& z = _data[2];
		return Three_Vector(x * c - y * s, y * c + x * s, z);
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

	Three_Vector Three_Vector::align(Three_Vector const& other) const
	{
		return Rz(-phi()).Ry(other.theta()-theta()).Rz(other.phi());
	}

	double Three_Vector::squared() const
	{
		return (_data[0] * std::conj(_data[0]) + _data[1] * std::conj(_data[1]) + _data[2] * std::conj(_data[2])).real();
	}

	double Three_Vector::magnitude() const
	{
		return std::sqrt(squared());
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
		auto a = std::atan2(rho, z);
		return a;
	}

	double Three_Vector::phi() const
	{
		auto const& x = _data[0].real();
		auto const& y = _data[1].real();
		if( x == 0 && y == 0 )
		{
			return 0.;
		}
		auto a = std::atan2(y, x);
		return a;
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

	double Three_Vector::cos_theta() const
	{
		if( magnitude() == 0 ){
			return 1.;
		}
		return _data[2].real() / magnitude();
	}

	double Three_Vector::cos_phi() const
	{
		if( _data[0].real() == 0 && _data[1].real() == 0 )
		{
			return 1;
		}
		Three_Vector projection(_data[0], _data[1], 0);
		return projection.x().real() / projection.magnitude();
	}
}