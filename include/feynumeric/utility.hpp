#ifndef Feynumeric_UTILITY_HPP
#define Feynumeric_UTILITY_HPP

#include <algorithm>
#include <iterator>
#include "feynumeric/complex.hpp"
#include "feynumeric/matrix.hpp"
#include "feynumeric/messages.hpp"

namespace Feynumeric
{
    template<typename T, typename U>
    bool contains(T const& container, U const& value)
    {
        return std::find(container.begin(), container.end(), value) != container.end();
    }

    bool almost_identical(double a, double b, double epsilon=1.e-9);

    template <auto Start, auto End, auto Inc, class F>
    constexpr void constexpr_for(F&& f)
    {
        if constexpr (Start < End)
        {
            f(std::integral_constant<decltype(Start), Start>());
            constexpr_for<Start + Inc, End, Inc>(f);
        }
    }

    Complex dot3(Matrix const& a, Matrix const& b)
    {
        return a.at(0) * b.at(0) + a.at(1) * b.at(1) + a.at(2) * b.at(2);
    }

    Complex dot4(Matrix const& a, Matrix const& b)
    {
        return a.at(0) * b.at(0) - a.at(1) * b.at(1) - a.at(2) * b.at(2) - a.at(3) * b.at(3);
    }
}
#endif // Feynumeric_UTILITY_HPP