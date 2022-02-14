#include "feynumeric/kinematics.hpp"

namespace Feynumeric
{

	Kinematics::Kinematics(double sqrt_s, std::size_t n_in, std::size_t n_out)
	: _sqrt_s(sqrt_s)
	, _n_in(n_in)
	, _n_out(n_out)
	{
		_momenta.resize(n_in + n_out);
	}

	double Kinematics::sqrt_s() const
	{
		return _sqrt_s;
	}

	Four_Momentum const& Kinematics::incoming(std::size_t i) const
	{
		return _momenta[i];
	}

	Four_Momentum const& Kinematics::outgoing(std::size_t i) const
	{
		return _momenta[_n_in + i];
	}

	Four_Momentum const& Kinematics::momentum(std::size_t i) const
	{
		return _momenta[i];
	}

	void Kinematics::incoming(std::size_t i, Four_Momentum const& p)
	{
		_momenta[i] = p;
	}

	void Kinematics::outgoing(std::size_t i, Four_Momentum const& p)
	{
		_momenta[_n_in+i] = p;
	}


}