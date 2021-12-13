#ifndef FEYNUMERIC_CONSTEXPR_MATH_HPP
#define FEYNUMERIC_CONSTEXPR_MATH_HPP

#include <array>

namespace Feynumeric
{
    constexpr double constexpr_sqrt(double x, double curr, double prev);
    constexpr double constexpr_sqrt(double x);

    constexpr std::size_t constexpr_binomial(std::size_t n, std::size_t k);
    constexpr std::size_t n_total_partitions(std::size_t j, std::size_t m);
    constexpr std::size_t constexpr_pow(std::size_t a, std::size_t n);
    constexpr std::size_t constexpr_factorial(std::size_t n);
    constexpr double abs(double x);
    constexpr bool is_angular_momentum(double x);
    constexpr bool is_half_integer(double x);
    constexpr bool is_integer(double x);

    constexpr double clebsch_gordan(std::size_t j1, int m1, std::size_t j2, int m2, std::size_t J, int M);

    template<std::size_t n>
    constexpr std::array<int, n> cummulated_sum(std::array<int, n> const& array)
    {
        std::array<int, n> accumulated = array;
        constexpr_for<1, n, 1>([&](auto k){
            accumulated[k] = accumulated[k] + accumulated[k-1];
        });
        return accumulated;
    }

    template<std::size_t n>
    constexpr double cummulated_clebsch_gordan(std::array<int, n> const& list)
    {
        double result = 1;
        auto accum = cummulated_sum(list);
        for( std::size_t k = 0; k < n-1; ++k)
        {
            result *= clebsch_gordan(k+1, accum[k], 1, list[k+1], k+2, accum[k+1]);
        }
        return result;
    }
}

#endif // FEYNCALC_CONSTEXPR_MATH_HPP