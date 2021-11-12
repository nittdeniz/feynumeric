#include <iostream>

#include "feyncalc/angular_momentum.hpp"

namespace Feyncalc
{
    Angular_Momentum::operator double() const
    {
        return 0;
    }

    Angular_Momentum::Angular_Momentum(double value)
    {
        if( is_valid_spin(value) )
        {
            _value = value;
        }
        else
        {
            std::cerr << "Angular Momentum is not of the form n/2 (int n >= 0)\n";
            std::abort();
        }
    }

    bool Angular_Momentum::is_valid_spin(double value)
    {
        return value >= 0 && (static_cast<int>(2*value) - 2*value) == 0;
    }

}
