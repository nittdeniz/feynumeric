#ifndef Feynumeric_DIRAC_HPP
#define Feynumeric_DIRAC_HPP

#include <array>

#include "feynumeric/angular_momentum.hpp"
#include "feynumeric/matrix.hpp"
#include "feynumeric/particle.hpp"
#include <functional>

namespace Feynumeric
{
    extern std::array<Matrix, 4> GA;
    extern std::array<Matrix, 4> GAC;

    constexpr double constexpr_sqrt(double x, double curr, double prev);
    constexpr double constexpr_sqrt(double x);


    template <auto Start, auto End, auto Inc, class F>
    constexpr void constexpr_for(F&& f)
    {
        if constexpr (Start < End)
        {
            f(std::integral_constant<decltype(Start), Start>());
            constexpr_for<Start + Inc, End, Inc>(f);
        }
    }

    /** constexpr binomials **/

    constexpr std::size_t binomial(std::size_t n, std::size_t k) {
        if (k == 0 || n == k) {
            return 1;
        }
        else {
            return binomial(n - 1, k - 1) + binomial(n - 1, k);
        }
    }

/** formula from the picture **/

    constexpr std::size_t n(std::size_t j, std::size_t m)
    {

        std::size_t result = 0;
        for (std::size_t k = m; k <= (j + m) / 2; ++k) {
            result += binomial(j, k) * binomial(j - k, k - m);
        }
        return result;
    }

/** constexpr power function **/

    constexpr std::size_t pow(std::size_t a, std::size_t n) {
        std::size_t r = 1;
        while (n--) {
            r *= a;
        }
        return r;
    }

/** actual function in question **/

    constexpr int abs(int m)
    {
        return m < 0? -m : m;
    }

    template<std::size_t j, int m>
    constexpr std::array<int, j*n(j, abs(m))> integer_partition()
    {
        std::array<int, j*n(j, abs(m))> result{};
        std::size_t i = 0;
        for (std::size_t k = 0; k < ::pow(3, j); ++k) {
            std::array<std::size_t, j> partition{};
            std::size_t sum = 0;
            for (std::size_t l = 0; l < j; ++l) {
                partition[l] = -((k/static_cast<std::size_t>(::pow(3, l)))%3-1);
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

    constexpr int factorial(int n)
    {
        if( n < 0 ) return 0;
        return n == 0? 1 : n * factorial(n-1);
    }

    template<std::size_t n>
    constexpr std::array<int, n> accumulate(std::array<int, n> const& array)
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
        auto accum = accumulate(list);
        for( std::size_t k = 0; k < n-1; ++k)
        {
            result *= clebsch_gordan(k+1, accum[k], list[k+1]);
        }
        return result;
    }

    template<std::size_t j, int lambda>
    std::complex<double> epsilon(double p, std::array<int, j> const& lorentz)
    {
        constexpr auto partition = integer_partition<j, lambda>();
        std::complex<double> result = 0;
        constexpr_for<0, partition.size(), j>([&](auto i){
            std::complex<double> val = 1;
            std::array<int, j> part{};
            constexpr_for<0, j, 1>([&](auto k){
                val *= epsilon<1, partition[k+i]>(p, {lorentz[k]});
                part[k] = partition[k+i];
            });
            result += cummulated_clebsch_gordan<j>(part) * val;
        });
        return result;
    }

    std::complex<double> epsilon(Angular_Momentum const& A, Four_Momentum const& p, std::vector<int> lorentz_indices)
    {

    }

    template<>
    std::complex<double> epsilon<1, 1>(double p, std::array<int, 1> const& lorentz)
    {
        using namespace std::complex_literals;
        constexpr std::array<std::complex<double>, 4> eps{0,-1,-1i,0};
        return 1./std::sqrt(2) * eps[lorentz[0]];
    }

    template<>
    std::complex<double> epsilon<1, 0>(double p, std::array<int, 1> const& lorentz)
    {
        constexpr std::array<std::complex<double>, 4> eps{0,0,0,1};
        return eps[lorentz[0]];
    }

    template<>
    std::complex<double> epsilon<1, -1>(double p, std::array<int, 1> const& lorentz)
    {
        using namespace std::complex_literals;
        constexpr std::array<std::complex<double>, 4> eps{0,1,-1i,0};
        return 1./std::sqrt(2) * eps[lorentz[0]];
        return -1.;
    }

    Matrix GS(Matrix const& a);

    Matrix dirac_sigma(Matrix const& a, Matrix const& b);

    Matrix epsilon1(Matrix const& p, Angular_Momentum const& lambda);
    Matrix epsilon2(Matrix const& p, Angular_Momentum const& lambda);

    Matrix u12(Matrix const& p, Angular_Momentum const& s);
    Matrix ubar12(Matrix const& p, Angular_Momentum const& s);

    Matrix u32(Matrix const& p, Angular_Momentum const& s);
    Matrix ubar32(Matrix const& p, Angular_Momentum const& s);

    Matrix u52(Matrix const& p, Angular_Momentum const& s);
    Matrix ubar52(Matrix const& p, Angular_Momentum const& s);

    Matrix v12(Matrix const& p, Angular_Momentum const& s);
    Matrix vbar12(Matrix const& p, Angular_Momentum const& s);

    Matrix v32(Matrix const& p, Angular_Momentum const& s);
    Matrix vbar32(Matrix const& p, Angular_Momentum const& s);

    Matrix v52(Matrix const& p, Angular_Momentum const& s);
    Matrix vbar52(Matrix const& p, Angular_Momentum const& s);

    Matrix Projector_12p(Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
    Matrix Projector_1(Matrix const& p, std::vector<std::size_t> const& lorentz_indices);

    Complex Breit_Wigner(Matrix const& p, double mass, std::function<double(Matrix const& p)> width);

    Matrix Propagator_0(Particle_Ptr const& particle, Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
    Matrix Propagator_12(Particle_Ptr const& particle, Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
    Matrix Propagator_1(Particle_Ptr const& particle, Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
    Matrix Propagator_32(Particle_Ptr const& particle, Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
    Matrix Propagator_2(Particle_Ptr const& particle, Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
    Matrix Propagator_52(Particle_Ptr const& particle, Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
}

#endif // Feynumeric_DIRAC_HPP


