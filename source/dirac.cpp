#include "dirac.hpp"
#include "complex.hpp"
#include "constexpr_math.hpp"
#include "lorentz_index.hpp"
#include "product.hpp"
#include "sum.hpp"
namespace Feynumeric
{

    std::array<Matrix, 4> GA =
            {
                Matrix(4, 4, {1,0,0,0, 0,1,0,0, 0,0,-1,0, 0,0,0,-1}),
                Matrix(4, 4, {0,0,0,1, 0,0,1,0, 0,-1,0,0, -1,0,0,0}),
                Matrix(4, 4, {0,0,0,-1.i, 0,0,1.i,0, 0,1.i,0,0, -1.i,0,0,0}),
                Matrix(4, 4, {0,0,1,0, 0,0,0,-1, -1,0,0,0, 0,1,0,0})
            };

    std::array<Matrix, 4> GAC =
            {
                    Matrix(4, 4, {1,0,0,0, 0,1,0,0, 0,0,-1,0, 0,0,0,-1}),
                    -Matrix(4, 4, {0,0,0,1, 0,0,1,0, 0,-1,0,0, -1,0,0,0}),
                    -Matrix(4, 4, {0,0,0,-1.i, 0,0,1.i,0, 0,1.i,0,0, -1.i,0,0,0}),
                    -Matrix(4, 4, {0,0,1,0, 0,0,0,-1, -1,0,0,0, 0,1,0,0})
            };

    [[maybe_unused]] Matrix GS(const Matrix &matrix)
    {
        Complex const& a = matrix.at(0);
        Complex const& b = matrix.at(1);
        Complex const& c = matrix.at(2);
        Complex const& d = matrix.at(3);
        return Matrix(4, 4,
                      std::vector<Complex>(std::initializer_list<Complex>{
                        a, 0., -d, -b + 1.i * c,
                        0., a, -b - 1.i * c, d,
                        d, b - 1.i * c, -a, 0.,
                        b + 1.i * c, -d, 0., -a
                      }));
    }

    [[maybe_unused]] Matrix dirac_sigma(const Matrix &a, const Matrix &b)
    {
        return Complex(0, 1)/2. * (a*b - b*a);
    }

    Matrix u(const Four_Momentum &p, const Angular_Momentum &s, const vector<Lorentz_Index*> &lorentz_indices)
    {
        if( s.j() == 1/2 )
        {
            if( s.m() == 1/2 )
            {
                return Matrix(1,1, {1.i});
            }
            if( s.m() == -1/2 )
            {
                return Matrix(1,1, {1.i});
            }
            critical_error("Invalid state in polarisation vector.\n");
        }
        Matrix result;
        for( double k = -1; k <= 1; k+= 2 )
        {
            double const n = k/2.;
            if( abs(s.m()) <= s.j() && abs(s.m() - n) <= s.j()-0.5 )
            {
                Angular_Momentum const s1 = Angular_Momentum(s.j()-0.5, s.m()-n);
                Angular_Momentum const s2 = Angular_Momentum(0.5, n);
                result +=
                        clebsch_gordan(s1.j(), s.m()-n, 1, n, s.j(), s.m())
                        * epsilon(p, s1, lorentz_indices)
                        * u(p, s2, {});
            }
        }
        return result;
    }

    Matrix ubar(const Four_Momentum &p, const Angular_Momentum &s, const vector<Lorentz_Index*> &lorentz_indices)
    {
        return u(p, s, lorentz_indices).T().apply([](Complex const& z){return std::conj(z);}) * GA[0];
    }

    Matrix epsilon(const Four_Momentum &p, const Angular_Momentum &s, const vector<Lorentz_Index*> &lorentz_indices)
    {
        if( s.j() == 0 )
        {
            return Matrix(1,1,1);
        }
        if( s.j() == 1 )
        {
            if( s.m() == 1 )
            {
                return Matrix(1,1, {1.i});
            }
            if( s.m() == 0 )
            {
                return Matrix(1,1, {1.i});
            }
            if( s.m() == -1 )
            {
                return Matrix(1,1, {1.i});
            }
            critical_error("Invalid state in polarisation vector.\n");
        }
        Matrix result;
        for( int n = -1; n <= 1; n++ )
        {
            if( abs(s.m()) <= s.j() && abs(s.m() - n) <= s.j()-1 )
            {
                Angular_Momentum const s1 = Angular_Momentum(s.j()-1, s.m()-n);
                Angular_Momentum const s2 = Angular_Momentum(1, n);
                auto indices1 = std::vector<Lorentz_Index*>(lorentz_indices.begin(), lorentz_indices.end()-1);
                auto indices2 = std::vector<Lorentz_Index*>(lorentz_indices.end()-1, lorentz_indices.end());
                result +=
                        clebsch_gordan(s1.j(), s1.m(), s2.j(), s2.m(), s.j(), s.m())
                        * epsilon(p, s1, indices1)
                        * epsilon(p, s2, indices2);
            }
        }
        return result;
    }

    Matrix epsilon_star(const Four_Momentum &p, const Angular_Momentum &s, const vector<Lorentz_Index*> &lorentz_indices)
    {
        return epsilon(p, s, lorentz_indices).apply([](Complex const& z){return std::conj(z);});
    }

    Matrix Projector(Particle_Ptr const& P, const Four_Momentum &p, const vector<Lorentz_Index*> &lorentz_indices, bool ignore_momentum)
    {
        static Matrix const metric_tensor = Matrix(4, 4, {1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,-1});

        auto const& f = constexpr_factorial;

        // algorithm from 10.1140/epjc/s2005-02299-4
        auto A = [&](int r, int n)
        {
            if( r == 0 )
            {
                return 1.;
            }
            return std::pow(-0.5, r) * f(n)/(1. * f(r) * f(n-2*r)) *
            product<double>([&](int k){ return 1./(2*n-(2*k+1));}, 0, r-1);
        };

        std::size_t const half_size = lorentz_indices.size() / 2;
        std::vector<Lorentz_Index*> indices_left(lorentz_indices.begin(), lorentz_indices.begin() + half_size);
        std::vector<Lorentz_Index*> indices_right(lorentz_indices.begin() + half_size, lorentz_indices.end());

        // since the whole projector is symmetric, we sort the indices so we can use std::next_permutation
        std::sort(indices_left.begin(), indices_left.end());
        std::sort(indices_right.begin(), indices_right.end());

        auto all_permutations = [&](std::vector<Lorentz_Index*>& v)
        {
            std::vector<Lorentz_Index*> result;
            result.reserve(v.size() * f(v.size()));
            do
            {
                result.insert(result.end(), v.begin(), v.end());
            } while( std::next_permutation(v.begin(), v.end()) );
            return result;
        };

        std::vector<Lorentz_Index*> mu = all_permutations(indices_left);
        std::vector<Lorentz_Index*> nu = all_permutations(indices_right);

        int const n = static_cast<int>(P->spin().j());

        auto T = [&](Lorentz_Index* mu, Lorentz_Index* nu){
            return -metric_tensor.at(*mu, *nu) + (ignore_momentum? 0. : 1./p.squared() * p.at(*mu) * p.at(*nu));
        };

        int const fac_n = f(n);
        return Matrix(1, 1, 1) * std::pow((1./fac_n), 2) *
          sum<Complex>(
            [&](int j)
            {
                return sum<Complex>([&](int i)
                {
                    return sum<Complex>([&](int r)
                    {
                        return A(r, n) *
                        product<Complex>([&](int l){
                            return T(mu[i*n + l], nu[j*n+l]);
                            }, 2*r + 1, n)
                        * product<Complex>([&](int k){
                            return T(mu[i*n + 2*k - 1], mu[i*n + 2*k]) * T(nu[j*n+2*k-1], nu[j*n+2*k]);
                        }, 1, 2*r-1);
                    }, 0, n%2 == 0? n/2 : (n-1)/2);
                }, 1, fac_n);
            }, 1, fac_n
        );
    }

    Matrix Propagator(const Particle_Ptr &P, const Four_Momentum &p, const vector<Lorentz_Index*> &lorentz_indices,
                      bool ignore_momentum)
    {
        return Projector(P, p, lorentz_indices, ignore_momentum);
    }
//
//    [[maybe_unused]] Matrix epsilon(const Matrix& p, const Angular_Momentum &lambda)
//    {
//        if( almost_identical(dot4(p, p), 0) )
//        {
//            if( lambda.j() == 1 )
//            {
//
//            }
//        }
//        return p;
//    }
//
//    [[maybe_unused]] Matrix epsilon1(const Matrix &p, const Angular_Momentum &lambda)
//    {
//        switch( static_cast<int>(lambda.m()) )
//        {
//            case 1:
//            case 0:
//            case -1:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    [[maybe_unused]] Matrix epsilon2(const Matrix &p, const Angular_Momentum &lambda)
//    {
//        switch( static_cast<int>(lambda.m()) )
//        {
//            case 2:
//            case 1:
//            case 0:
//            case -1:
//            case -2:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    [[maybe_unused]] Matrix u12(const Matrix &p, const Angular_Momentum &s)
//    {
//        switch( static_cast<int>(2*s.m()) )
//        {
//            case 1:
//            case -1:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    [[maybe_unused]] Matrix ubar12(const Matrix &p, const Angular_Momentum &s)
//    {
//        switch( static_cast<int>(2*s.m()) )
//        {
//            case 1:
//            case -1:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    [[maybe_unused]] Matrix u32(const Matrix &p, const Angular_Momentum &s)
//    {
//        switch( static_cast<int>(2*s.m()) )
//        {
//            case 3:
//            case 1:
//            case -1:
//            case -3:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    [[maybe_unused]] Matrix ubar32(const Matrix &p, const Angular_Momentum &s)
//    {
//        switch( static_cast<int>(2*s.m()) )
//        {
//            case 3:
//            case 1:
//            case -1:
//            case -3:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    [[maybe_unused]] Matrix u52(const Matrix &p, const Angular_Momentum &s)
//    {
//        switch( static_cast<int>(2*s.m()) )
//        {
//            case 5:
//            case 3:
//            case 1:
//            case -1:
//            case -3:
//            case -5:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    [[maybe_unused]] Matrix ubar52(const Matrix &p, const Angular_Momentum &s)
//    {
//        switch( static_cast<int>(2*s.m()) )
//        {
//            case 5:
//            case 3:
//            case 1:
//            case -1:
//            case -3:
//            case -5:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    [[maybe_unused]] Matrix v12(const Matrix &p, const Angular_Momentum &s)
//    {
//        switch( static_cast<int>(2*s.m()) )
//        {
//            case 1:
//            case -1:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    [[maybe_unused]] Matrix vbar12(const Matrix &p, const Angular_Momentum &s)
//    {
//        switch( static_cast<int>(2*s.m()) )
//        {
//            case 1:
//            case -1:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    [[maybe_unused]] Matrix v32(const Matrix &p, const Angular_Momentum &s)
//    {
//        switch( static_cast<int>(2*s.m()) )
//        {
//            case 3:
//            case 1:
//            case -1:
//            case -3:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    [[maybe_unused]] Matrix vbar32(const Matrix &p, const Angular_Momentum &s)
//    {
//        switch( static_cast<int>(2*s.m()) )
//        {
//            case 3:
//            case 1:
//            case -1:
//            case -3:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    [[maybe_unused]] Matrix v52(const Matrix &p, const Angular_Momentum &s)
//    {
//        switch( static_cast<int>(2*s.m()) )
//        {
//            case 5:
//            case 3:
//            case 1:
//            case -1:
//            case -3:
//            case -5:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    [[maybe_unused]] Matrix vbar52(const Matrix &p, const Angular_Momentum &s)
//    {
//        switch( static_cast<int>(2*s.m()) )
//        {
//            case 5:
//            case 3:
//            case 1:
//            case -1:
//            case -3:
//            case -5:
//            default:
//                break;
//        }
//        return Matrix(p);
//    }
//
//    std::complex<double> epsilon(Angular_Momentum const& A, Four_Momentum const& p, std::vector<int> lorentz_indices)
//    {
//        return A.j() * p.E() * lorentz_indices.size();
//    }
//
//    Matrix Projector_12p(const Matrix &p, const vector<std::size_t> &lorentz_indices)
//    {
//        return Matrix();
//    }
//
//    Matrix Projector_1(const Matrix &p, const vector<std::size_t> &lorentz_indices)
//    {
//        return Matrix();
//    }
//
//    Complex Breit_Wigner(const Matrix &p, double mass, std::function<double(const Matrix &)> width)
//    {
//        return Feynumeric::Complex();
//    }
//
//    Matrix Propagator_0(const Particle_Ptr &particle, const Matrix &p, const vector<std::size_t> &lorentz_indices)
//    {
//        return Matrix();
//    }
//
//    Matrix Propagator_12(const Particle_Ptr &particle, const Matrix &p, const vector<std::size_t> &lorentz_indices)
//    {
//        return Matrix();
//    }
//
//    Matrix Propagator_1(const Particle_Ptr &particle, const Matrix &p, const vector<std::size_t> &lorentz_indices)
//    {
//        return Matrix();
//    }
//
//    Matrix Propagator_32(const Particle_Ptr &particle, const Matrix &p, const vector<std::size_t> &lorentz_indices)
//    {
//        return Matrix();
//    }
//
//    Matrix Propagator_2(const Particle_Ptr &particle, const Matrix &p, const vector<std::size_t> &lorentz_indices)
//    {
//        return Matrix();
//    }
//
//    Matrix Propagator_52(const Particle_Ptr &particle, const Matrix &p, const vector<std::size_t> &lorentz_indices)
//    {
//        return Matrix();
//    }
}