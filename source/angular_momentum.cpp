#include <iostream>

#include "feynumeric/angular_momentum.hpp"
#include "feynumeric/messages.hpp"

namespace Feynumeric
{
    Angular_Momentum::Angular_Momentum(double j, double m, bool massless)
    : _j(j)
    , _m(m)
    , _massless(massless)
    {
        if( !is_valid_spin(j) || !is_valid_spin(std::abs(m)) )
        {
        	critical_error("Angular Momentum is not of the form n/2 (int n >= 0)");
        }
        if( std::abs(m) > j )
        {
            critical_error("|m| > j in Angular Momentum.");
        }
        if( massless && std::abs(m) != j )
        {
        	critical_error("Massless particle must have |m| = j.");
        }
    }

	Angular_Momentum::Angular_Momentum(const Angular_Momentum& J)
	: _j(J.j())
	, _m(J.m())
	, _massless(massless)
	{

	}

	Angular_Momentum& Angular_Momentum::operator=(const Angular_Momentum& J){
		_j = J._j;
		_m = J._m;
		_massless = J._massless;
		return *this;
	}

	bool Angular_Momentum::is_valid_spin(double value)
    {
        return value >= 0 && (static_cast<int>(2*value) - 2*value) == 0;
    }

    bool Angular_Momentum::is_half_odd_integer() const
    {
        return (static_cast<int>(_j) - _j) != 0;
    }

    double Angular_Momentum::j() const
    {
        return _j;
    }

    double Angular_Momentum::m() const
    {
        return _m;
    }

	bool Angular_Momentum::massless() const
	{
		return _massless;
	}

	void Angular_Momentum::m(double new_m)
    {
        if( std::abs(new_m) > _j || (static_cast<int>(_j-new_m) != _j-new_m) )
        {
            std::cerr << "Invalid value for AngularMomentum::m: (" << _j << " " << new_m << ")\n";
            std::abort();
        }
        _m = new_m;
    }

	void Angular_Momentum::reset()
	{
		_m = _j;
	}

	std::size_t Angular_Momentum::n_states() const
	{
		return (_massless && _j > 0 )? 2 : 2 * _j + 1;
	}

	Angular_Momentum& operator++(Angular_Momentum& J)
	{
		J._m++;
		if( J._m > J._j )
		{
			J._m = -J._j;
		}
		return *this;
	}
}
