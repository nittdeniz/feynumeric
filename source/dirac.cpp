#include <iostream>

#include "dirac.hpp"
#include "complex.hpp"
#include "constexpr_math.hpp"
#include "feynman_graph.hpp"
#include "lorentz_index.hpp"
#include "product.hpp"
#include "sum.hpp"
namespace Feynumeric
{

	std::array<Matrix, 4> GA =
			{
					Matrix(4, 4, {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1}),
					Matrix(4, 4, {0, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 0, 0}),
					Matrix(4, 4, {0, 0, 0, -1.i, 0, 0, 1.i, 0, 0, 1.i, 0, 0, -1.i, 0, 0, 0}),
					Matrix(4, 4, {0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 1, 0, 0})
			};

	std::array<Matrix, 4> GAC =
			{
					Matrix(4, 4, {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1}),
					-Matrix(4, 4, {0, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 0, 0}),
					-Matrix(4, 4, {0, 0, 0, -1.i, 0, 0, 1.i, 0, 0, 1.i, 0, 0, -1.i, 0, 0, 0}),
					-Matrix(4, 4, {0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 1, 0, 0})
			};

	Matrix GA5 = Matrix(4, 4, {0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0});

	std::array<std::array<double, 4>, 4> MT = {{{1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, -1, 0}, {0, 0, 0, -1}}};

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

	[[maybe_unused]] Matrix GS(const Matrix& matrix)
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

	[[maybe_unused]] Matrix dirac_sigma(const Matrix& a, const Matrix& b)
	{
		return Complex(0, 1) / 2. * ( a * b - b * a );
	}

	Matrix dirac_sigmac(Lorentz_Index_Ptr mu, Lorentz_Index_Ptr nu)
	{
		return dirac_sigma(GAC[*mu], GAC[*nu]);
	}

	Matrix dirac_sigmac(Lorentz_Index_Ptr mu, Four_Vector const& p)
	{
		return dirac_sigma(GAC[*mu], GS(p));
	}

	Matrix dirac_sigmac(Four_Vector const& p, Lorentz_Index_Ptr nu)
	{
		return dirac_sigma(GS(p), GAC[*nu]);
	}

	Matrix u(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin)
	{
		return u(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->spin(), edge_ptr->lorentz_indices());
	}

	Matrix v(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin)
	{
		return v(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->spin(), edge_ptr->lorentz_indices());
	}

	Matrix v(Particle_Ptr const& P, const Four_Vector& p, Angular_Momentum_Ptr s,
	         const std::vector<Lorentz_Index_Ptr>& lorentz_indices)
	{
		if( p.E().real() < 0 ){
			Four_Vector p_new = -p;
			return v(P, p_new, s, lorentz_indices);
		}
		if( s->j() == 1. / 2 ){
			auto N = p.E() + P->mass();
			if( s->m() == 1. / 2 ){
				auto result = std::sqrt(N) * Matrix(4, 1, {( p.x() - 1.i * p.y()) / N, -p.z() / N, 0, 1});
				return result;
			}
			if( s->m() == -1. / 2 ){
				return std::sqrt(N) * Matrix(4, 1, {p.z() / N, ( p.x() + 1.i * p.y()) / N, 1, 0});
			}
			critical_error("Invalid state in polarisation vector.\n");
		}
		Matrix result(4, 1);
		for( double k = -1; k <= 1; k += 2 ){
			double const n = k / 2.;
			if( abs(s->m()) <= s->j() && abs(s->m() - n) <= s->j() - 0.5 ){
				Angular_Momentum_Ptr s1 = std::make_shared<Angular_Momentum>(
						Angular_Momentum(s->j() - 0.5, s->m() - n));
				Angular_Momentum_Ptr s2 = std::make_shared<Angular_Momentum>(Angular_Momentum(0.5, n));
				result +=
						clebsch_gordan(s1->j(), s->m() - n, 1, n, s->j(), s->m())
						* epsilon(P, p, s1, lorentz_indices)
						* v(P, p, s2, {});
			}
		}
		return result;
	}

	Matrix u(Particle_Ptr const& P, const Four_Vector& p, Angular_Momentum_Ptr s,
	         const std::vector<Lorentz_Index_Ptr>& lorentz_indices)
	{
		if( p.E().real() < 0 ){
			Four_Vector p_new = -p;
			return u(P, p_new, s, lorentz_indices);
		}
		if( s->j() == 1/2. ){
			auto N = p.E() + P->mass();
			if( s->m() == 1/2. ){
				auto result = std::sqrt(N) * Matrix(4, 1, {1, 0, p.z() / N, ( p.x() + 1.i * p.y()) / N});
				return result;
			}
			if( s->m() == -1/2. ){
				return std::sqrt(N) * Matrix(4, 1, {0, 1, ( p.x() - 1.i * p.y()) / N, -p.z() / N});
			}
			critical_error("Invalid state in polarisation vector.\n");
		}
		Matrix result(4, 1);
		for( double k = -1; k <= 1; k += 2 ){
			double const n = k / 2.;
			if( abs(s->m()) <= s->j() && abs(s->m() - n) <= s->j() - 0.5 ){
				Angular_Momentum_Ptr s1 = std::make_shared<Angular_Momentum>(
						Angular_Momentum(s->j() - 0.5, s->m() - n));
				Angular_Momentum_Ptr s2 = std::make_shared<Angular_Momentum>(Angular_Momentum(0.5, n));
				auto clebsch = clebsch_gordan(s1->j(), s1->m(), s2->j(), s2->m(), s->j(), s->m());
				auto eps = epsilon(P, p, s1, lorentz_indices);
				auto uu = u(P, p, s2, {});;
				result += clebsch * eps * uu;
			}
		}
		return result;
	}

	Matrix vbar(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin)
	{
		return vbar(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->spin(), edge_ptr->lorentz_indices());
	}

	Matrix vbar(Particle_Ptr const& P, const Four_Vector& p, Angular_Momentum_Ptr s,
	            const std::vector<Lorentz_Index_Ptr>& lorentz_indices)
	{
		return v(P, p, s, lorentz_indices).T().apply([](Complex const& z){ return std::conj(z); }) * GA[0];
	}

	Matrix ubar(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin)
	{
		return ubar(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->spin(), edge_ptr->lorentz_indices());
	}

	Matrix ubar(Particle_Ptr const& P, const Four_Vector& p, Angular_Momentum_Ptr s,
	            const std::vector<Lorentz_Index_Ptr>& lorentz_indices)
	{
		return u(P, p, s, lorentz_indices).T().apply([](Complex const& z){ return std::conj(z); }) * GA[0];
	}

	Matrix epsilon(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin)
	{
		return epsilon(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->spin(),
		               edge_ptr->lorentz_indices());
	}

	Matrix epsilon_star(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin)
	{
		return epsilon_star(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->spin(),
		                    edge_ptr->lorentz_indices());
	}

    Complex epsilon(Angular_Momentum_Ptr const& spin, double q, double mass, double cos_theta, double cos_phi, std::vector<Lorentz_Index_Ptr> const& mus)
    {
		if( spin->j() == 1 )
		{
			static std::array<std::function<std::vector<Complex>(double, double, double, double)>,3> f_i =
					{
							[](double, double, double ct, double cp){
								double st = std::sqrt(1-ct*ct);
								double sp = std::sqrt(1-cp*cp);
								return std::vector<Complex>{0, 1./constexpr_sqrt(2.) * (ct * cp + 1.i * sp), 1./constexpr_sqrt(2.) * (ct * sp - 1.i * cp), -1./constexpr_sqrt(2.) * (st)};
							},
							[](double q, double m, double ct, double cp){
								double st = std::sqrt(1-ct*ct);
								double sp = std::sqrt(1-cp*cp);
								double im = 1./m;
								double Em = std::sqrt(q*q + m*m) * im;
								return std::vector<Complex>{q*im, Em * st * cp, Em * st * sp, Em * ct};
							},
							[](double, double, double ct, double cp){
								double st = std::sqrt(1-ct*ct);
								double sp = std::sqrt(1-cp*cp);
								return std::vector<Complex>{0, 1./constexpr_sqrt(2.) * (-ct * cp + 1.i * sp), 1./constexpr_sqrt(2.)*(-ct * sp - 1.i * cp), 1./constexpr_sqrt(2.) * st};
							}
					};
			auto vector = (f_i[spin->m() + 1](q, mass, cos_theta, cos_phi));
			return vector[*(mus[0])];
		}
		Complex result{0., 0.};
	    for( int n = -1; n <= 1; n++ )
	    {
		    if( abs(spin->m()) <= spin->j() && abs(spin->m() - n) <= spin->j() - 1 )
		    {
			    Angular_Momentum_Ptr s1 = std::make_shared<Angular_Momentum>(Angular_Momentum(spin->j() - 1, spin->m() - n));
			    Angular_Momentum_Ptr s2 = std::make_shared<Angular_Momentum>(Angular_Momentum(1, n));
			    auto indices1 = std::vector<Lorentz_Index_Ptr>(mus.begin(), mus.end()-1);
			    auto indices2 = std::vector<Lorentz_Index_Ptr>(mus.end()-1, mus.end());
			    result +=
					    clebsch_gordan(s1->j(), s1->m(), s2->j(), s2->m(), spin->j(), spin->m())
					    * epsilon(s1, q, mass, cos_theta, cos_phi, indices1)
					    * epsilon(s2, q, mass, cos_theta, cos_phi, indices2);
		    }
	    }
	    return result;
    }

    Matrix O32(Four_Vector const& p, Lorentz_Index_Ptr const& mu, Lorentz_Index_Ptr const& lambda)
	{
		return p.contra(mu) * GA[*(lambda)] - GS(p) * MT[*mu][*lambda];
	}

	Matrix O32c(Four_Vector const& p, Lorentz_Index_Ptr const& mu, Lorentz_Index_Ptr const& lambda)
	{
		return p.co(mu) * GAC[*(lambda)] - GS(p) * MT[*mu][*lambda];
	}

	Matrix O(Angular_Momentum_Ptr const& s, Four_Vector const& p, std::vector<Lorentz_Index_Ptr> const& mu, std::vector<Lorentz_Index_Ptr> lambda)
	{
		if( mu.size() != lambda.size() )
		{
			critical_error("Projector 0 requires indices mu and lambda to have same length.");
		}
		if( s->j()-0.5 != mu.size() )
		{
			critical_error(FORMAT("Projector O with spin {} requires {} indices but {} were given.", s->j(), s->j()-0.5, mu.size()));
		}

		Matrix result(4, 4);
		// since the whole projector is symmetric, we sort the indices so we can use std::next_permutation
		std::vector<std::vector<Lorentz_Index_Ptr>> permutations;
		std::size_t const n = constexpr_factorial(lambda.size());
		permutations.reserve(n);
		do{
			permutations.push_back(lambda);
		} while( std::next_permutation(lambda.begin(), lambda.end()) );

		for( auto const& permutation : permutations )
		{
			Matrix temp(4, 4, 1);
			for( std::size_t i = 0; i < lambda.size(); ++i )
			{
				temp *= O32(p, mu[i], permutation[i]);
			}
			result += temp;
		}
		return result;
	}

	Matrix Oc(Angular_Momentum_Ptr const& s, Four_Vector const& p, std::vector<Lorentz_Index_Ptr> const& mu, std::vector<Lorentz_Index_Ptr> lambda)
	{
		if( mu.size() != lambda.size() )
		{
			critical_error("Projector 0 requires indices mu and lambda to have same length.");
		}
		if( s->j()-0.5 != mu.size() )
		{
			critical_error(FORMAT("Projector O with spin {} requires {} indices but {} were given.", s->j(), s->j()-0.5, mu.size()));
		}

		Matrix result(4, 4);
		// since the whole projector is symmetric, we sort the indices so we can use std::next_permutation
		std::vector<std::vector<Lorentz_Index_Ptr>> permutations;
		std::size_t const n = constexpr_factorial(lambda.size());
		permutations.reserve(n);
		do{
			permutations.push_back(lambda);
		} while( std::next_permutation(lambda.begin(), lambda.end()) );

		for( auto const& permutation : permutations )
		{
			Matrix temp(4, 4, 1);
			for( std::size_t i = 0; i < lambda.size(); ++i )
			{
				temp *= O32c(p, mu[i], permutation[i]);
			}
			result += temp;
		}
		return result;
	}


	Matrix epsilon(Particle_Ptr const& P, const Four_Vector &pp, Angular_Momentum_Ptr s, const std::vector<Lorentz_Index_Ptr> &lorentz_indices)
    {
		Four_Vector const p = pp.E().real() > 0 ? pp : -pp;
        if( s->j() == 0 )
        {
            return Matrix(1,1,1);
        }
        if( P->mass() == 0 && std::abs(s->m()) != s->j())
        {
        	critical_error("Massless particle must have |m| = j.\n");
        }
        auto const qq = p.spatial().squared();
	    auto const q = std::sqrt(qq);
        return epsilon(s, q, P->mass(), p.cos_theta(), p.cos_phi(), lorentz_indices);
    }


    Matrix epsilon_star(Particle_Ptr const& P, const Four_Vector &p, Angular_Momentum_Ptr s, const std::vector<Lorentz_Index_Ptr> &lorentz_indices)
    {
        return epsilon(P, p, s, lorentz_indices).apply([](Complex const& z){return std::conj(z);});
    }



	Matrix Projector(std::shared_ptr<Graph_Edge> edge_ptr, const Kinematics& kin, bool ignore_momentum){
		return Projector(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->lorentz_indices(), ignore_momentum);
	}

	Matrix Propagator(std::shared_ptr<Graph_Edge> edge_ptr, const Kinematics& kin, bool ignore_momentum){
		return Propagator(edge_ptr->particle(), edge_ptr->four_momentum(kin), edge_ptr->lorentz_indices(), ignore_momentum);
	}

	Matrix Projector(Particle_Ptr const& P, const Four_Vector &p, const std::vector<Lorentz_Index_Ptr> &lorentz_indices, bool ignore_momentum)
	{
		return Projector(P, P->spin(), p, lorentz_indices, ignore_momentum);
	}

	Matrix Projector(Particle_Ptr const& P, Angular_Momentum const& spin, const Four_Vector &p, const std::vector<Lorentz_Index_Ptr> &lorentz_indices, bool ignore_momentum)
    {
	    int const n = static_cast<int>(spin.j());
	    if( spin.is_half_odd_integer() )
	    {
	    	Angular_Momentum new_spin(spin.j() + 0.5, spin.j() + 0.5);
	    	auto copy = lorentz_indices;
	    	auto mu = std::make_shared<Lorentz_Index>();
	    	auto nu = std::make_shared<Lorentz_Index>();
	    	copy.insert(copy.begin() + copy.size()/2, nu);
	    	copy.insert(copy.begin(), mu);
	    	Matrix contraction(4, 4, 0);
	    	for( std::size_t i = 0; i < 4; ++i )
		    {
	    	    for( std::size_t j = 0; j < 4; ++j )
		        {
					contraction += -(GAC[*mu] * GAC[*nu] * Projector(P, new_spin, p, copy, ignore_momentum));
					++(*nu);
		        }
	    	    ++(*mu);
		    }
		    return ( P->is_true_fermion() ? ( GS(p) + P->mass()) : ( GS(p) - P->mass()) )
		                      * ( n + 1. ) / ( 2. * n + 3. ) * contraction;
	    }

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


        auto T = [&](Lorentz_Index_Ptr const& mu, Lorentz_Index_Ptr const& nu){
        	if( mu == nullptr || nu == nullptr )
	        {
        		critical_error("Lorentz_Index_Ptr is nullptr.\n");
	        }
	        auto const mt = metric_tensor.at(*mu, *nu);
        	if( ignore_momentum )
	        {
				return -mt;
	        }
        	auto pmu = p.contra(mu);
        	auto pnu = p.contra(nu);
        	auto denominator = p.squared();
        	auto numerator = pmu * pnu;
			auto momentum_contribution = (numerator.real() == 0 && numerator.imag() == 0 && denominator == 0) ? 1 : numerator/denominator;
        	auto result = -mt + momentum_contribution;
            return result;
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
		Complex breit_wigner;
		double p2 = p.squared();
		if( p2 > 0 )
		{
			breit_wigner = -1.i/(p.squared() - P->mass() * P->mass() + 1.i * std::sqrt(p2) * P->width(p2));
		}
		else
		{
			breit_wigner = -1.i/(p.squared() - P->mass() * P->mass());
		}
        return  breit_wigner *  Projector(P, p, lorentz_indices, ignore_momentum);
    }
}