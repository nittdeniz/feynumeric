#ifndef FEYNUMERIC_CONSTEXPR_MATH_HPP
#define FEYNUMERIC_CONSTEXPR_MATH_HPP

#include <array>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace Feynumeric
{
    template<typename T>
    constexpr T abs(T x)
    {
        return x > 0
               ? x
               : -x;
    }

    template<typename T>
    constexpr T min(T a, T b)
    {
        return a < b ? a : b;
    }
    template<typename T>
    constexpr T max(T a, T b)
    {
        return a > b? a : b;
    }

    template<typename T>
    constexpr bool is_almost_equal(T a, T b, T rel_epsilon = 1.e-15)
    {
        if( std::isnan(a) || std::isnan(b) )
        {
            return false;
        }
        return abs(a-b) < std::abs(max(a, b)) * rel_epsilon;
    }

    constexpr double constexpr_sqrt(double x, double curr, double prev)
    {
        return curr == prev
               ? curr
               : constexpr_sqrt(x, 0.5 * (curr + x / curr), curr);
    }



    constexpr double constexpr_sqrt(double x)
    {
        return x >= 0 && x < std::numeric_limits<double>::infinity()
               ? constexpr_sqrt(x, x, 0)
               : std::numeric_limits<double>::quiet_NaN();
    }

    constexpr std::size_t constexpr_binomial(std::size_t n, std::size_t k)
    {
        if (k == 0 || n == k) {
            return 1;
        }
        else {
            return constexpr_binomial(n - 1, k - 1) + constexpr_binomial(n - 1, k);
        }
    }
//    constexpr std::size_t n_total_partitions(std::size_t j, std::size_t m);

    template<typename T>
    constexpr T constexpr_pow(T a, std::size_t n) {
        T r = 1;
        while( n-- ){
            r *= a;
        }
        return r;
    }

    constexpr std::size_t constexpr_factorial(std::size_t n)
    {
        return n == 0? 1 : n * constexpr_factorial(n-1);
    }


    constexpr bool is_angular_momentum(double x)
    {
        return static_cast<int>(2*x) - 2*x == 0;
    }

    constexpr bool is_integer(double x)
    {
        return static_cast<int>(x) - x == 0;
    }

    constexpr bool is_half_integer(double x)
    {
        return !is_integer(x) && is_integer(2*x);
    }


    constexpr double clebsch_gordan(double j1, double m1, double j2, double m2, double J, double M)
    {
        if( !is_angular_momentum(j1) || !is_angular_momentum(m1) || !is_angular_momentum(j2) || !is_angular_momentum(m2) || !is_angular_momentum(J) )
        {
            throw std::domain_error("clebsch_gordan: all parameters must be multiples of 0.5.");
        }
        if( (is_integer(j1) && is_half_integer(m1)) || (is_half_integer(j1) && is_integer(m1)) )
        {
            throw std::domain_error("clebsch_gordan: j1 and m1 must both be integral or half integral.");
        }
        if( (is_integer(j2) && is_half_integer(m2)) || (is_half_integer(j2) && is_integer(m2)) )
        {
            throw std::domain_error("clebsch_gordan: j2 and m2 must both be integral or half integral.");
        }
        if( (is_integer(j1+j2) && !is_integer(J)) || (is_half_integer(j1+j2) && !is_half_integer(J)) )
        {
            throw std::domain_error("clebsch_gordan: J is no valid value for j1 and j2 combination.");
        }
        if( j1 < 0 || j2 < 0 )
        {
            throw std::domain_error("clebsch_gordan: j1 and j2 must be non-negative.");
        }
        if( J > j1+j2 || J < abs(j1-j2) )
        {
            throw std::domain_error("clebsch_gordan: J must be |j1-j2| <= J <= j1+j2.");
        }
        if( abs(m1) > j1 || abs(m2) > j2 || abs(M) > J )
        {
            throw std::domain_error("clebsch_gordan: m1 and m2 must be |m1| <= j1 and |m2| <= j2 and |m1+m2| <= J");
        }
        if( m1 + m2 != M )
        {
            throw std::domain_error("clebsch_gordan: m1 + m2 != M.");
        }

        auto const& f = constexpr_factorial;
        double numerator = (2*J+1) * f(J+j1-j2) * f(J-j1+j2) * f(j1+j2-J);
        numerator *= f(J+M) * f(J-M) * f(j1-m1) * f(j1+m1) * f(j2-m2) * f(j2+m2);
        double denominator = f(j1+j2+J+1);
        int const min = std::max(0., std::max(j2-J-m1, j1+m2-J));
        int const max = std::min(j2+m2, std::min(j1-m1, j1+j2-J));
        double sum = min > max? 1 : 0;
        for( int k = min; k <= max; ++k )
        {
            sum += constexpr_pow(-1., k) / ( f(k) * f(j1+j2-J-k) * f(j1-m1-k) * f(j2+m2-k) * f(J-j2+m1+k) * f(J-j1-m2+k) );
        }
        return constexpr_sqrt(numerator / denominator) * sum;
    }

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