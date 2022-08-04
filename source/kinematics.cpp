#include "feynumeric/kinematics.hpp"

namespace Feynumeric
{

	Kinematics::Kinematics(double sqrt_s, std::size_t n_in, std::size_t n_out)
	: _sqrt_s(sqrt_s)
	, _n_in(n_in)
	, _n_out(n_out)
	{
		_momenta.resize(n_in + n_out);
        _angles.resize(n_out-1);
	}

	double Kinematics::sqrt_s() const
	{
		return _sqrt_s;
	}

	Four_Vector const& Kinematics::incoming(std::size_t i) const
	{
		return _momenta[i];
	}

	Four_Vector const& Kinematics::outgoing(std::size_t i) const
	{
		return _momenta[_n_in + i];
	}

	Four_Vector const& Kinematics::momentum(std::size_t i) const
	{
		return _momenta[i];
	}

	void Kinematics::incoming(std::size_t i, Four_Vector const& p)
	{
		_momenta[i] = p;
	}

	void Kinematics::outgoing(std::size_t i, Four_Vector const& p)
	{
		_momenta[_n_in+i] = p;
	}

    void Kinematics::angle(std::size_t i, double a)
    {
        _angles[i] = a;
    }

    double Kinematics::angle(std::size_t i) const
    {
        return _angles[i];
    }


}