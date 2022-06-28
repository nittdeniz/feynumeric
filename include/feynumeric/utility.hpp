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
    bool almost_identical(Complex a, Complex b, double epsilon=1.e-9);

	double int_pow(double c, std::size_t n);
    Complex int_pow(Complex c, std::size_t n);

//    constexpr std::size_t binomial(std::size_t n, std::size_t k) {
//        if (k == 0 || n == k) {
//            return 1;
//        }
//        else {
//            return binomial(n - 1, k - 1) + binomial(n - 1, k);
//        }
//    }

/** formula from the picture **/

//    constexpr std::size_t n(std::size_t j, std::size_t m)
//    {
//
//        std::size_t result = 0;
//        for (std::size_t k = m; k <= (j + m) / 2; ++k) {
//            result += binomial(j, k) * binomial(j - k, k - m);
//        }
//        return result;
//    }
/** constexpr power function **/

//    constexpr std::size_t pow(std::size_t a, std::size_t n) {
//        std::size_t r = 1;
//        while (n--) {
//            r *= a;
//        }
//        return r;
//    }

/** actual function in question **/

//    constexpr int abs(int m)
//    {
//        return m < 0? -m : m;
//    }

//    std::vector<int> integer_partition(std::size_t j, int m)
//    {
//        std::vector<int> result;
//        result.resize(j * n(j, abs(m)));
//        std::size_t i = 0;
//        for (std::size_t k = 0; k < ::pow(3, j); ++k) {
//            std::vector<int> partition(j);
//            int sum = 0;
//            for (std::size_t l = 0; l < j; ++l) {
//                partition[l] = -((k/static_cast<std::size_t>(::pow(3, l)))%3-1);
//                sum += partition[l];
//            }
//            if (sum == m)
//            {
//                for (std::size_t l = 0; l < j; ++l) {
//                    result[j*i + l] = partition[l];
//                }
//                ++i;
//            }
//        }
//        return result;
//    }

//    constexpr int factorial(int n)
//    {
//        if( n < 0 ) return 0;
//        return n == 0? 1 : n * factorial(n-1);
//    }

//    template <auto Start, auto End, auto Inc, class F>
//    constexpr void constexpr_for(F&& f)
//    {
//        if constexpr (Start < End)
//        {
//            f(std::integral_constant<decltype(Start), Start>());
//            constexpr_for<Start + Inc, End, Inc>(f);
//        }
//    }
//
//    constexpr double constexpr_sqrt(double x, double curr, double prev);
//    constexpr double constexpr_sqrt(double x);
//
///*
//    template <auto Start, auto End, auto Inc, class F>
//    constexpr void constexpr_for(F&& f)
//    {
//        if constexpr (Start < End)
//        {
//            f(std::integral_constant<decltype(Start), Start>());
//            constexpr_for<Start + Inc, End, Inc>(f);
//        }
//    }
//*/
//    /** constexpr binomials **/
//
//    constexpr std::size_t binomial(std::size_t n, std::size_t k) {
//        if (k == 0 || n == k) {
//            return 1;
//        }
//        else {
//            return binomial(n - 1, k - 1) + binomial(n - 1, k);
//        }
//    }
//
///** formula from the picture **/
//
//    constexpr std::size_t n(std::size_t j, std::size_t m)
//    {
//
//        std::size_t result = 0;
//        for (std::size_t k = m; k <= (j + m) / 2; ++k) {
//            result += binomial(j, k) * binomial(j - k, k - m);
//        }
//        return result;
//    }
//
///** constexpr power function **/
//
//    constexpr std::size_t pow(std::size_t a, std::size_t n) {
//        std::size_t r = 1;
//        while (n--) {
//            r *= a;
//        }
//        return r;
//    }
//
///** actual function in question **/
//
//    constexpr int abs(int m)
//    {
//        return m < 0? -m : m;
//    }
//
//    template<std::size_t j, int m>
//    constexpr std::array<int, j*n(j, abs(m))> integer_partition()
//    {
//        std::array<int, j*n(j, abs(m))> result{};
//        std::size_t i = 0;
//        for (std::size_t k = 0; k < ::pow(3, j); ++k) {
//            std::array<std::size_t, j> partition{};
//            std::size_t sum = 0;
//            for (std::size_t l = 0; l < j; ++l) {
//                partition[l] = -((k/static_cast<std::size_t>(::pow(3, l)))%3-1);
//                sum += partition[l];
//            }
//            if (sum == m)
//            {
//                for (std::size_t l = 0; l < j; ++l) {
//                    result[j*i + l] = partition[l];
//                }
//                ++i;
//            }
//        }
//        return result;
//    }
//
//    constexpr int factorial(int n)
//    {
//        if( n < 0 ) return 0;
//        return n == 0? 1 : n * factorial(n-1);
//    }
//
//    template<std::size_t n>
//    constexpr std::array<int, n> accumulate(std::array<int, n> const& array)
//    {
//        std::array<int, n> accumulated = array;
//        constexpr_for<1, n, 1>([&](auto k){
//            accumulated[k] = accumulated[k] + accumulated[k-1];
//        });
//        return accumulated;
//    }
//
//    template<std::size_t n>
//    constexpr double cummulated_clebsch_gordan(std::array<int, n> const& list)
//    {
//        double result = 1;
//        auto accum = accumulate(list);
//        for( std::size_t k = 0; k < n-1; ++k)
//        {
//            result *= clebsch_gordan(k+1, accum[k], list[k+1]);
//        }
//        return result;
//    }

    Complex dot3(Matrix const& a, Matrix const& b);
    Complex dot4(Matrix const& a, Matrix const& b);
}
#endif // Feynumeric_UTILITY_HPP