#include "utility.hpp"

namespace Feynumeric
{
    bool almost_identical(double a, double b, double epsilon)
    {
        return std::abs(a-b) < std::abs(epsilon * std::min(a, b));
    }
}

