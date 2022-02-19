#include "four_vector.hpp"
#include "utility.hpp"

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

	Four_Vector::Four_Vector(Complex t, Three_Vector const& r)
	: _temporal(t)
	, _spatial(r)
	{

	}

	Three_Vector const& Four_Vector::spatial() const
	{
		return _spatial;
	}

	Four_Vector Four_Vector::Rx(double phi) const
	{
		return Four_Vector(_temporal, _spatial.Rx(phi));
	}

	Four_Vector Four_Vector::Ry(double phi) const
	{
		return Four_Vector(_temporal, _spatial.Ry(phi));
	}

	Four_Vector Four_Vector::Rz(double phi) const
	{
		return Four_Vector(_temporal, _spatial.Rz(phi));
	}

	Four_Vector Four_Vector::boost(Four_Vector const& p) const
	{
		auto pp = p.squared();
		if( almost_identical(pp, 0) ){
			critical_error("Boosting with light-like momentum.");
		}
		if( pp < 0 )
		{
			return boost(-p);
		}
		auto const& t = _temporal;
		auto const& x = _spatial.at(0);
		auto const& y = _spatial.at(1);
		auto const& z = _spatial.at(2);

		auto const beta = p.spatial()/p.E();
		auto const beta2 = beta.squared();
		auto const gamma = 1./std::sqrt(1-beta.squared());

		auto const Ltt = gamma;
		auto const Lxx = (gamma - 1) * beta.x() * beta.x() / beta2 + 1.;
		auto const Lyy = (gamma - 1) * beta.y() * beta.y() / beta2 + 1.;
		auto const Lzz = (gamma - 1) * beta.z() * beta.z() / beta2 + 1.;

		auto const Ltx = beta.x() * gamma;
		auto const Lty = beta.y() * gamma;
		auto const Ltz = beta.z() * gamma;

		auto const Lxy = (gamma - 1) * beta.x() * beta.y() / beta2;
		auto const Lyz = (gamma - 1) * beta.y() * beta.z() / beta2;
		auto const Lzx = (gamma - 1) * beta.z() * beta.x() / beta2;

		auto const Lxt = Ltx;
		auto const Lyt = Lty;
		auto const Lzt = Ltz;

		auto const Lyx = Lxy;
		auto const Lzy = Lyz;
		auto const Lxz = Lzx;

		return Four_Vector( Ltt * t + Ltx * x + Lty * y + Ltz * z,
					        Lxt * t + Lxx * x + Lxy * y + Lxz * z,
					        Lyt * t + Lyx * x + Lyy * y + Lyz * z,
					        Lzt * t + Lzx * x + Lzy * y + Lzz * z);
	}

	Four_Vector Four_Vector::align(Four_Vector const& p) const
	{
		Four_Vector copy(*this);
		copy._spatial.align(p.spatial());
		return copy;
	}

	Four_Vector Four_Vector::align(Three_Vector const& p) const
	{
		Four_Vector copy(*this);
		copy._spatial.align(p);
		return copy;
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

	double Four_Vector::theta() const
	{
		return _spatial.theta();
	}

	double Four_Vector::phi() const
	{
		return _spatial.phi();
	}

	double Four_Vector::spatial_squared() const
	{
		return _spatial.squared();
	}

	double Four_Vector::squared() const
	{
		return (_temporal * std::conj(_temporal) - spatial_squared()).real();
	}

	Complex Four_Vector::co(Lorentz_Index mu) const
	{
		return mu == 0 ? _temporal : -_spatial.at(static_cast<int>(mu) - 1);
	}

	Complex Four_Vector::co(Lorentz_Index_Ptr mu) const
	{
		return *mu == 0 ? _temporal : -_spatial.at(static_cast<int>(*mu) - 1);
	}

	Complex Four_Vector::contra(Lorentz_Index mu) const
	{
		return mu == 0 ? _temporal : _spatial.at(static_cast<int>(mu) - 1);
	}

	Complex Four_Vector::contra(Lorentz_Index_Ptr mu) const
	{
		return *mu == 0 ? _temporal : _spatial.at(static_cast<int>(*mu) - 1);
	}

	Four_Vector operator+(Four_Vector const& lhs, Four_Vector const& rhs)
	{
		Four_Vector copy(lhs);
		copy._temporal += rhs._temporal;
		copy._spatial += rhs._spatial;
		return copy;
	}

	Four_Vector operator-(Four_Vector const& lhs)
	{
		Four_Vector copy(lhs);
		copy._temporal *= -1;
		copy._spatial *= -1;
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

	template Four_Vector operator*(int const& lhs, Four_Vector const& rhs);
	template Four_Vector operator*(double const& lhs, Four_Vector const& rhs);
	template Four_Vector operator*(Complex const& lhs, Four_Vector const& rhs);

	template Four_Vector operator*(Four_Vector const& rhs, int const& lhs);
	template Four_Vector operator*(Four_Vector const& rhs, double const& lhs);
	template Four_Vector operator*(Four_Vector const& rhs, Complex const& lhs);

	Complex dot(Four_Vector const& lhs, Four_Vector const& rhs)
	{
		return lhs._temporal * rhs._temporal - dot(lhs._spatial, rhs._spatial);
	}

	Complex Four_Vector::t() const
	{
		return _temporal;
	}

	Complex Four_Vector::E() const
	{
		return _temporal;
	}

	Complex Four_Vector::x() const
	{
		return _spatial.x();
	}

	Complex Four_Vector::y() const
	{
		return _spatial.y();
	}

	Complex Four_Vector::z() const
	{
		return _spatial.z();
	}
}