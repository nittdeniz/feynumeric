#include "constexpr_math.hpp"
#include <array>
#include <stdexcept>
#include <limits>


namespace Feynumeric
{
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

    constexpr std::size_t constexpr_binomial(std::size_t n, std::size_t k) {
        if (k == 0 || n == k) {
            return 1;
        }
        else {
            return constexpr_binomial(n - 1, k - 1) + constexpr_binomial(n - 1, k);
        }
    }

    constexpr std::size_t n_total_partitions(std::size_t j, std::size_t m)
    {

        std::size_t result = 0;
        for (std::size_t k = m; k <= (j + m) / 2; ++k) {
            result += constexpr_binomial(j, k) * constexpr_binomial(j - k, k - m);
        }
        return result;
    }

    constexpr std::size_t constexpr_pow(std::size_t a, std::size_t n) {
        std::size_t r = 1;
        while (n--) {
            r *= a;
        }
        return r;
    }

    constexpr int abs(int m)
    {
        return m < 0? -m : m;
    }

    template<std::size_t j, int m>
    constexpr std::array<int, j*n_total_partitions(j, abs(m))> integer_partition()
    {
        std::array<int, j*n_total_partitions(j, abs(m))> result{};
        std::size_t i = 0;
        for (std::size_t k = 0; k < constexpr_pow(3, j); ++k) {
            std::array<std::size_t, j> partition{};
            std::size_t sum = 0;
            for (std::size_t l = 0; l < j; ++l) {
                partition[l] = -((k/static_cast<std::size_t>(constexpr_pow(3, l)))%3-1);
                sum += partition[l];
            }
            if (sum == m)
            {
                for (std::size_t l = 0; l < j; ++l) {
                    result[j*i + l] = partition[l];
                }
                ++i;
            }
        }
        return result;
    }

    constexpr std::size_t constexpr_factorial(std::size_t n)
    {
        return n == 0? 1 : n * constexpr_factorial(n-1);
    }

    constexpr bool is_integer(double x)
    {
        return static_cast<int>(x) - x == 0;
    }

    constexpr bool is_half_integer(double x)
    {
        return !is_integer(x) && is_integer(2*x);
    }

    constexpr bool is_angular_momentum(double x)
    {
        return static_cast<int>(2*x) - 2*x == 0;
    }

    constexpr double abs(double x)
    {
        return x > 0
               ? x
               : -x;
    }


    double constexpr clebsch_gordan(double const j1, double const m1, double const j2, double const m2, double const J)
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
        double const  M = m1+m2;
        if( abs(m1) > j1 || abs(m2) > j2 || abs(M) > J )
        {
            throw std::domain_error("clebsch_gordan: m1 and m2 must be |m1| <= j1 and |m2| <= j2 and |m1+m2| <= J");
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
            sum += constexpr_pow(-1, k)*1. / ( f(k) * f(j1+j2-J-k) * f(j2+m2-k) * f(J-j2+m1+k) * f(J-j1-m2+k) );
        }
        return constexpr_sqrt(numerator / denominator) * sum;
    }


}