#ifndef Feynumeric_DIRAC_HPP
#define Feynumeric_DIRAC_HPP

#include <array>

#include "feynumeric/angular_momentum.hpp"
#include "feynumeric/boost.hpp"
#include "feynumeric/matrix.hpp"
#include "feynumeric/momentum.hpp"
#include "feynumeric/particle.hpp"
#include <functional>

namespace Feynumeric
{
    extern std::array<Matrix, 4> GA;
    extern std::array<Matrix, 4> GAC;


    Matrix u(Four_Momentum const& p, Angular_Momentum const& s, std::vector<Lorentz_Index*> const& lorentz_indices);
    Matrix ubar(Four_Momentum const& p, Angular_Momentum const& s, std::vector<Lorentz_Index*> const& lorentz_indices);

    Matrix epsilon(Four_Momentum const& p, Angular_Momentum const& s, std::vector<Lorentz_Index*> const& lorentz_indices);
    Matrix epsilon_star(Four_Momentum const& p, Angular_Momentum const& s, std::vector<Lorentz_Index*> const& lorentz_indices);

    Matrix Projector(Particle_Ptr const& P, const Four_Momentum &p, const vector<Lorentz_Index*> &lorentz_indices, bool ignore_momentum = false);
    Matrix Propagator(Particle_Ptr const& P, const Four_Momentum &p, const vector<Lorentz_Index*> &lorentz_indices, bool ignore_momentum = false);



//
//    template<std::size_t j, int lambda>
//    inline std::complex<double> epsilon(double p, std::array<int, j> const& lorentz)
//    {
//        constexpr auto partition = integer_partition<j, lambda>();
//        std::complex<double> result = 0;
//        constexpr_for<0, partition.size(), j>([&](auto i){
//            std::complex<double> val = 1;
//            std::array<int, j> part{};
//            constexpr_for<0, j, 1>([&](auto k){
//                val *= epsilon<1, partition[k+i]>(p, {lorentz[k]});
//                part[k] = partition[k+i];
//            });
//            result += cummulated_clebsch_gordan<j>(part) * val;
//        });
//        return result;
//    }

//    std::complex<double> epsilon(Angular_Momentum const& A, Four_Momentum const& p, std::vector<int> lorentz_indices);
//
//    template<>
//    inline std::complex<double> epsilon<1, 1>(double p, std::array<int, 1> const& lorentz)
//    {
//        using namespace std::complex_literals;
//        Matrix epsilon(4, 1, {0, -1, -1i, 0});
//        auto boosted = boost(p, epsilon);
//        return 1./std::sqrt(2) * epsilon[lorentz[0]];
//    }
//
//    template<>
//    inline std::complex<double> epsilon<1, 0>(double p, std::array<int, 1> const& lorentz)
//    {
//        Matrix epsilon(4, 1, {0,0,0,1});
//        auto boosted = boost(p, epsilon);
//        return  epsilon[lorentz[0]];
//    }
//
//    template<>
//    inline std::complex<double> epsilon<1, -1>(double p, std::array<int, 1> const& lorentz)
//    {
//        using namespace std::complex_literals;
//        Matrix epsilon(4, 1, {0, 1, -1i, 0});
//        auto boosted = boost(p, epsilon);
//        return 1./std::sqrt(2) * epsilon[lorentz[0]];
//    }
//
//    Matrix GS(Matrix const& a);
//
//    Matrix dirac_sigma(Matrix const& a, Matrix const& b);
//
//    Matrix epsilon1(Matrix const& p, Angular_Momentum const& lambda);
//    Matrix epsilon2(Matrix const& p, Angular_Momentum const& lambda);
//
//    Matrix u12(Matrix const& p, Angular_Momentum const& s);
//    Matrix ubar12(Matrix const& p, Angular_Momentum const& s);
//
//    Matrix u32(Matrix const& p, Angular_Momentum const& s);
//    Matrix ubar32(Matrix const& p, Angular_Momentum const& s);
//
//    Matrix u52(Matrix const& p, Angular_Momentum const& s);
//    Matrix ubar52(Matrix const& p, Angular_Momentum const& s);
//
//    Matrix v12(Matrix const& p, Angular_Momentum const& s);
//    Matrix vbar12(Matrix const& p, Angular_Momentum const& s);
//
//    Matrix v32(Matrix const& p, Angular_Momentum const& s);
//    Matrix vbar32(Matrix const& p, Angular_Momentum const& s);
//
//    Matrix v52(Matrix const& p, Angular_Momentum const& s);
//    Matrix vbar52(Matrix const& p, Angular_Momentum const& s);
//
//    Matrix Projector_12p(Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
//    Matrix Projector_1(Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
//
//    Complex Breit_Wigner(Matrix const& p, double mass, std::function<double(Matrix const& p)> width);
//
//    Matrix Propagator_0(Particle_Ptr const& particle, Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
//    Matrix Propagator_12(Particle_Ptr const& particle, Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
//    Matrix Propagator_1(Particle_Ptr const& particle, Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
//    Matrix Propagator_32(Particle_Ptr const& particle, Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
//    Matrix Propagator_2(Particle_Ptr const& particle, Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
//    Matrix Propagator_52(Particle_Ptr const& particle, Matrix const& p, std::vector<std::size_t> const& lorentz_indices);
}

#endif // Feynumeric_DIRAC_HPP


