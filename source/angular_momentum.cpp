#include <iostream>

#include "Feynumeric/angular_momentum.hpp"

namespace Feynumeric
{
    Angular_Momentum::Angular_Momentum(double j, double m)
    : _j(j)
    , _m(m)
    {
        if( !is_valid_spin(j) || !is_valid_spin(std::abs(m)) )
        {
            std::cerr << "Angular Momentum is not of the form n/2 (int n >= 0)\n";
            std::abort();
        }
        if( std::abs(m) > j )
        {
            std::cerr << "m > j in Angular Momentum.\n";
            std::abort();
        }
    }

	Angular_Momentum::Angular_Momentum(const Angular_Momentum& J)
	: _j(J.j())
	, _m(J.m())
	{

	}

	Angular_Momentum& Angular_Momentum::operator=(const Angular_Momentum& J){
		_j = J.j();
		_m = J.m();
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

    void Angular_Momentum::m(double new_m)
    {
        if( std::abs(new_m) > _j || (static_cast<int>(_j-new_m) != _j-new_m) )
        {
            std::cerr << "Invalid value for AngularMomentum::m: (" << _j << " " << new_m << ")\n";
            std::abort();
        }
        _m = new_m;
    }
}
