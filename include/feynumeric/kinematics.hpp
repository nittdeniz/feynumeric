#ifndef Feynumeric_KINEMATICS_HPP
#define Feynumeric_KINEMATICS_HPP

#include <vector>

namespace Feynumeric
{
    struct Kinematics
    {
        double sqrt_s;
        std::vector<double> momenta;
        std::vector<double> invariant_masses;
        std::vector<double> cosines;
    };
}

#endif // Feynumeric_KINEMATICS_HPP