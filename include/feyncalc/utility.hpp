#ifndef FEYNCALC_UTILITY_HPP
#define FEYNCALC_UTILITY_HPP

#include <algorithm>
#include <iterator>

namespace Feyncalc
{
    template<typename T, typename U>
    bool contains(T const& container, U const& value)
    {
        return std::find(container.begin(), container.end(), value) != container.end();
    }

    bool almost_identical(double a, double b, double epsilon=1.e-9);
}
#endif // FEYNCALC_UTILITY_HPP