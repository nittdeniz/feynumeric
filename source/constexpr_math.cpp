#include "constexpr_math.hpp"
#include <array>
#include <stdexcept>
#include <limits>


namespace Feynumeric
{


//    constexpr std::size_t constexpr_binomial(std::size_t n, std::size_t k)

    constexpr std::size_t n_total_partitions(std::size_t j, std::size_t m)
    {

        std::size_t result = 0;
        for (std::size_t k = m; k <= (j + m) / 2; ++k) {
            result += constexpr_binomial(j, k) * constexpr_binomial(j - k, k - m);
        }
        return result;
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
}