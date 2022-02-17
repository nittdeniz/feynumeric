#include "four_vector.hpp"

namespace Feynumeric
{
	Four_Vector::Four_Vector(Complex t, Complex x, Complex y, Complex z)
	: _temporal(t)
	, _spatial(x, y, z)
	{

	}

	Four_Vector::Four_Vector(Four_Vector const& copy)
	: _temporal(copy._temporal)
	  , _spatial(copy._spatial)
	{

	}

	Four_Vector& Four_Vector::operator=(Four_Vector const& copy)
	{
		_temporal = copy._temporal;
		_spatial = copy._spatial;
		return *this;
	}

	Four_Vector& Four_Vector::operator+=(Four_Vector const& rhs)
	{
		_temporal += rhs._temporal;
		_spatial += rhs._spatial;
		return *this;
	}

	Four_Vector& Four_Vector::operator-=(Four_Vector const& rhs)
	{
		_temporal -= rhs._temporal;
		_spatial -= rhs._spatial;
		return *this;
	}

	template<typename T>
	Four_Vector& Four_Vector::operator*=(T const& rhs)
	{
		_temporal *= rhs;
		_spatial *= rhs;
		return *this;
	}

	template<typename T>
	Four_Vector& Four_Vector::operator/=(T const& rhs)
	{
		_temporal /= rhs;
		_spatial /= rhs;
		return *this;
	}

	double Four_Vector::theta()
	{
		return _spatial.theta();
	}

	double Four_Vector::phi()
	{
		return _spatial.phi();
	}

	double Four_Vector::spatial_square()
	{
		return spatial_square();
	}

	double Four_Vector::square()
	{
		return (_temporal * std::conj(_temporal) - spatial_square()).real();
	}

	Complex Four_Vector::co(std::size_t i) const
	{
		return i == 0 ? _temporal : -_spatial.at(i);
	}

	Complex Four_Vector::contra(std::size_t i) const
	{
		return i == 0 ? _temporal : _spatial.at(i);
	}

	Four_Vector operator+(Four_Vector const& lhs, Four_Vector const& rhs)
	{
		Four_Vector copy(lhs);
		copy._temporal += rhs._temporal;
		copy._spatial += rhs._spatial;
		return copy;
	}

	Four_Vector operator-(Four_Vector const& lhs, Four_Vector const& rhs)
	{
		Four_Vector copy(lhs);
		copy._temporal -= rhs._temporal;
		copy._spatial -= rhs._spatial;
		return copy;
	}

	template<typename T>
	Four_Vector operator*(Four_Vector const& lhs, T const& rhs)
	{
		Four_Vector copy(lhs);
		copy._temporal *= rhs;
		copy._spatial *= rhs;
		return copy;
	}

	template<typename T>
	Four_Vector operator*(T const& lhs, Four_Vector const& rhs)
	{
		return rhs * lhs;
	}

	template<typename T>
	Four_Vector operator/(Four_Vector const& lhs, T const& rhs)
	{
		Four_Vector copy(lhs);
		copy._temporal /= rhs;
		copy._spatial /= rhs;
		return copy;
	}

	Complex dot(Four_Vector const& lhs, Four_Vector const& rhs)
	{
		return lhs._temporal * rhs._temporal - dot(lhs._spatial, rhs._spatial);
	}
}