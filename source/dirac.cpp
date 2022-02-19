#include <iostream>
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

    std::array<std::array<double, 4>, 4> MT = {{{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}}};

	[[maybe_unused]] Matrix GS(Four_Vector const& p)
	{
		Complex const& a = p.E();
		Complex const& b = p.x();
		Complex const& c = p.y();
		Complex const& d = p.z();
		return Matrix(4, 4,
		              std::vector<Complex>(std::initializer_list<Complex>{
				              a, 0., -d, -b + 1.i * c,
				              0., a, -b - 1.i * c, d,
				              d, b - 1.i * c, -a, 0.,
				              b + 1.i * c, -d, 0., -a
		              }));
	}

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

    Matrix u(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin)
    {
		return u(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->spin(), edge_ptr->lorentz_indices());
    }

	Matrix v(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin)
	{
		return v(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->spin(), edge_ptr->lorentz_indices());
	}

	Matrix v(Particle_Ptr const& P, const Four_Vector &p, Angular_Momentum_Ptr s, const std::vector<Lorentz_Index_Ptr> &lorentz_indices)
	{
		if( p.E().real() < 0 )
		{
			Four_Vector p_new = -p;
			return u(P, p_new, s, lorentz_indices);
		}
		if( s->j() == 1./2 )
		{
			auto N = p.E()+P->mass();
			if( s->m() == 1./2 )
			{
				auto result = std::sqrt(N) * Matrix(4,1, {(p.x() - 1.i * p.y())/N, -p.z()/N, 0, 1});
				return result;
			}
			if( s->m() == -1./2 )
			{
				return std::sqrt(N) * Matrix(4,1, {p.z()/N, (p.x()+1.i * p.y())/N, 1, 0});
			}
			critical_error("Invalid state in polarisation vector.\n");
		}
		Matrix result(4, 1);
		for( double k = -1; k <= 1; k+= 2 )
		{
			double const n = k/2.;
			if( abs(s->m()) <= s->j() && abs(s->m() - n) <= s->j()-0.5 )
			{
				Angular_Momentum_Ptr s1 = std::make_shared<Angular_Momentum>(Angular_Momentum(s->j()-0.5, s->m()-n));
				Angular_Momentum_Ptr s2 = std::make_shared<Angular_Momentum>(Angular_Momentum(0.5, n));
				result +=
						clebsch_gordan(s1->j(), s->m()-n, 1, n, s->j(), s->m())
						* epsilon(P, p, s1, lorentz_indices)
						* v(P, p, s2, {});
			}
		}
		return result;
	}

    Matrix u(Particle_Ptr const& P, const Four_Vector &p, Angular_Momentum_Ptr s, const std::vector<Lorentz_Index_Ptr> &lorentz_indices)
    {
		if( p.E().real() < 0 )
		{
			Four_Vector p_new = -p;
			return u(P, p_new, s, lorentz_indices);
		}
        if( s->j() == 1./2 )
        {
        	auto N = p.E()+P->mass();
            if( s->m() == 1./2 )
            {
            	auto result = std::sqrt(N) * Matrix(4,1, {1, 0, p.z()/N, (p.x()+1.i * p.y())/N});
                return result;
            }
            if( s->m() == -1./2 )
            {
                return std::sqrt(N) * Matrix(4,1, {0, 1, (p.x()-1.i * p.y())/N, -p.z()/N});
            }
            critical_error("Invalid state in polarisation vector.\n");
        }
        Matrix result(4, 1);
        for( double k = -1; k <= 1; k+= 2 )
        {
            double const n = k/2.;
            if( abs(s->m()) <= s->j() && abs(s->m() - n) <= s->j()-0.5 )
            {
                Angular_Momentum_Ptr s1 = std::make_shared<Angular_Momentum>(Angular_Momentum(s->j()-0.5, s->m()-n));
                Angular_Momentum_Ptr s2 = std::make_shared<Angular_Momentum>(Angular_Momentum(0.5, n));
                result +=
                        clebsch_gordan(s1->j(), s->m()-n, 1, n, s->j(), s->m())
                        * epsilon(P, p, s1, lorentz_indices)
                        * u(P, p, s2, {});
            }
        }
        return result;
    }

	Matrix vbar(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin)
	{
		return vbar(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->spin(), edge_ptr->lorentz_indices());
	}

	Matrix vbar(Particle_Ptr const& P, const Four_Vector &p, Angular_Momentum_Ptr s, const std::vector<Lorentz_Index_Ptr> &lorentz_indices)
	{
		return v(P, p, s, lorentz_indices).T().apply([](Complex const& z){return std::conj(z);}) * GA[0];
	}

    Matrix ubar(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin)
	{
        return ubar(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->spin(), edge_ptr->lorentz_indices());
	}



    Matrix ubar(Particle_Ptr const& P, const Four_Vector &p, Angular_Momentum_Ptr s, const std::vector<Lorentz_Index_Ptr> &lorentz_indices)
    {
        return u(P, p, s, lorentz_indices).T().apply([](Complex const& z){return std::conj(z);}) * GA[0];
    }

    Matrix epsilon(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin)
    {
    	return epsilon(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->spin(), edge_ptr->lorentz_indices());
    }

    Matrix epsilon_star(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin)
    {
	    return epsilon_star(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->spin(), edge_ptr->lorentz_indices());
    }

    Matrix epsilon(Particle_Ptr const& P, const Four_Vector &pp, Angular_Momentum_Ptr s, const std::vector<Lorentz_Index_Ptr> &lorentz_indices)
    {
		Four_Vector const p = pp.E().real() > 0 ? pp : -pp;
        if( s->j() == 0 )
        {
            return Matrix(1,1,1);
        }

        auto align = [&p](Four_Vector const& epsilon)
        {
        	auto theta = p.theta();
        	auto phi = p.phi();
        	auto y = epsilon.Ry(theta);
        	auto z = y.Rz(phi);
        	return z;
        };

        if( P->mass() == 0 )
        {
	        if( s->j() == 1 )
	        {
		        double constexpr isqrt2 = 1./constexpr_sqrt(2.);
		        if( s->m() == 1 )
		        {
			        static auto const vector = Four_Vector(0, -isqrt2, -1.i * isqrt2, 0);
			        auto aligned = align(vector);
			        return aligned.contra(lorentz_indices[0]);
		        }
		        if( s->m() == -1 )
		        {
			        static auto const vector = Four_Vector(0, isqrt2, -1.i * isqrt2, 0);
			        auto aligned = align(vector);
			        return aligned.contra(lorentz_indices[0]);
		        }
		        critical_error("Invalid state in polarisation vector.\n");
	        }
        }
        else
        {
	        if( s->j() == 1 )
	        {
		        double constexpr isqrt2 = 1./constexpr_sqrt(2.);
		        if( s->m() == 1 )
		        {
			        static auto const vector = Four_Vector(0, -isqrt2, -1.i * isqrt2, 0);
			        auto aligned = align(vector);
			        return aligned.contra(lorentz_indices[0]);
		        }
		        if( s->m() == 0 )
		        {
			        auto const vector = Four_Vector(std::sqrt(p.spatial_squared())/P->mass(),0,0,p.E()/P->mass());
			        auto aligned = align(vector);
			        return aligned.contra(lorentz_indices[0]);
		        }
		        if( s->m() == -1 )
		        {
			        static auto const vector = Four_Vector(0, isqrt2, -1.i * isqrt2, 0);
			        auto aligned = align(vector);
			        return aligned.contra(lorentz_indices[0]);
		        }
		        critical_error("Invalid state in polarisation vector.\n");
	        }
        }
        Matrix result;
        for( int n = -1; n <= 1; n++ )
        {
            if( abs(s->m()) <= s->j() && abs(s->m() - n) <= s->j()-1 )
            {
                Angular_Momentum_Ptr s1 = std::make_shared<Angular_Momentum>(Angular_Momentum(s->j()-1, s->m()-n));
                Angular_Momentum_Ptr s2 = std::make_shared<Angular_Momentum>(Angular_Momentum(1, n));
                auto indices1 = std::vector<Lorentz_Index_Ptr>(lorentz_indices.begin(), lorentz_indices.end()-1);
                auto indices2 = std::vector<Lorentz_Index_Ptr>(lorentz_indices.end()-1, lorentz_indices.end());
                result +=
                        clebsch_gordan(s1->j(), s1->m(), s2->j(), s2->m(), s->j(), s->m())
                        * epsilon(P, p, s1, indices1)
                        * epsilon(P, p, s2, indices2);
            }
        }
        return result;
    }


    Matrix epsilon_star(Particle_Ptr const& P, const Four_Vector &p, Angular_Momentum_Ptr s, const std::vector<Lorentz_Index_Ptr> &lorentz_indices)
    {
        return epsilon(P, p, s, lorentz_indices).apply([](Complex const& z){return std::conj(z);});
    }



	Matrix Projector(Feynman_Graph::Edge_Ptr edge_ptr, const Kinematics& kin, bool ignore_momentum){
		return Projector(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->lorentz_indices(), ignore_momentum);
	}

	Matrix Propagator(Feynman_Graph::Edge_Ptr edge_ptr, const Kinematics& kin, bool ignore_momentum){
		return Propagator(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->lorentz_indices(), ignore_momentum);
	}

	Matrix Projector(Particle_Ptr const& P, const Four_Vector &p, const std::vector<Lorentz_Index_Ptr> &lorentz_indices, bool ignore_momentum)
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
        std::vector<Lorentz_Index_Ptr> indices_left(lorentz_indices.begin(), lorentz_indices.begin() + half_size);
        std::vector<Lorentz_Index_Ptr> indices_right(lorentz_indices.begin() + half_size, lorentz_indices.end());

        // since the whole projector is symmetric, we sort the indices so we can use std::next_permutation
        std::sort(indices_left.begin(), indices_left.end());
        std::sort(indices_right.begin(), indices_right.end());

        auto all_permutations = [&](std::vector<Lorentz_Index_Ptr>& v)
        {
            std::vector<Lorentz_Index_Ptr> result;
            result.reserve(v.size() * f(v.size()));
            do
            {
                result.insert(result.end(), v.begin(), v.end());
            } while( std::next_permutation(v.begin(), v.end()) );
            return result;
        };

        std::vector<Lorentz_Index_Ptr> mu = all_permutations(indices_left);
        std::vector<Lorentz_Index_Ptr> nu = all_permutations(indices_right);

        int const n = static_cast<int>(P->spin().j());

        auto T = [&](Lorentz_Index_Ptr const& mu, Lorentz_Index_Ptr const& nu){
        	if( mu == nullptr || nu == nullptr )
	        {
        		critical_error("Lorentz_Index_Ptr is nullptr.\n");
	        }
        	auto pmu = p.co(mu);
        	auto pnu = p.co(nu);
        	auto denominator = p.squared();
        	auto mt = metric_tensor.at(*mu, *nu);
        	auto numerator = pmu * pnu;
			auto momentum_contribution = (numerator.real() == 0 && numerator.imag() == 0 && denominator == 0) ? 1 : numerator/denominator;
        	auto result = -mt + (ignore_momentum? 0. : momentum_contribution);
            return result;
        };

        Matrix dirac_structure = Matrix(1,1,1);
        if( P->spin().is_half_odd_integer() )
        {
			dirac_structure = GS(p) + P->mass();
        }

        int const fac_n = f(n);
        return dirac_structure * Matrix(1, 1, 1) * std::pow((1./fac_n), 2) *
          sum<Complex>(
            [&](int j)
            {
                return sum<Complex>([&](int i)
                {
                    return sum<Complex>([&](int r)
                    {
                        return A(r, n) *
                        product<Complex>([&](int k){
                            return T(mu[i*n + k -1], nu[j*n+k -1]);
                            }, 2*r + 1, n)
                        * product<Complex>([&](int k){
                            return T(mu[i*n + 2*k - 1 -1], mu[i*n + 2*k -1]) * T(nu[j*n+2*k-1 -1], nu[j*n+2*k -1]);
                        }, 1, r);
                    }, 0, n%2 == 0? n/2 : (n-1)/2);
                }, 0, fac_n-1);
            }, 0, fac_n-1
        );
    }

	Matrix dirac_sigma(Lorentz_Index_Ptr mu, Lorentz_Index_Ptr nu)
	{
		return GA[*mu]*GA[*nu] - GA[*nu]*GA[*mu];
	}

	Matrix dirac_sigma(Lorentz_Index_Ptr mu, Four_Vector const& p)
	{
		return GA[*mu] * GS(p) - GS(p) * GA[*mu];
	}

	Matrix dirac_sigma(Four_Vector const& p, Lorentz_Index_Ptr nu)
	{
		return GS(p) * GA[*nu] - GA[*nu] * GS(p);
	}

	Matrix Propagator(const Particle_Ptr &P, const Four_Vector &p, const std::vector<Lorentz_Index_Ptr> &lorentz_indices,
                      bool ignore_momentum)
    {
        return Projector(P, p, lorentz_indices, ignore_momentum);
    }
}