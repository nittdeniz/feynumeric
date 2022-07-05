#ifndef FEYNUMERIC_AMPLITUDE_HPP
#define FEYNUMERIC_AMPLITUDE_HPP

#include "functions.hpp"
#include "particle.hpp"
#include "polynomial.hpp"
#include "phase_space.hpp"

#include <vector>

namespace Feynumeric
{
	template<std::size_t N>
	class Amplitude
	{
	private:
		std::vector<FPolynomial<N>> _fpolynomials;
		std::vector<Particle_Ptr> _incoming;
		std::vector<Particle_Ptr> _virtual;
		std::vector<Particle_Ptr> _outgoing;

	public:
		Amplitude(std::vector<std::vector<Polynomial>> const& polynomials, std::vector<Particle_Ptr> const& incoming, std::vector<Particle_Ptr> const& virtuals, std::vector<Particle_Ptr> const& outgoing);

		void scale(func_t<N> const& func){
			for( auto& fp : _fpolynomials ){
				for( auto& f : fp._coefficients ){
					f = f * func;
				}
			}
		}

		double operator()(double x, auto&& args...){
			double result{0.};
			for( auto& fp : _fpolynomials ){
				Complex temp = fp(x, args);
				std::cout << "temp: " << temp << "\n";
				result += (temp * std::conj(temp)).real();
			}
			return result;
		}

		std::map<double, double> width(std::vector<double> const& s_values){
			using namespace std::placeholders;
			int N_polarisations = 1;
			for( auto const& p : _incoming ){
				N_polarisations *= p->spin().n_states();
			}

			std::map<double, double> result;

			if( N == 0 ){
				auto q = std::bind(momentum, _1, _outgoing[0]->mass(), _outgoing[1]->mass());
				auto phase_space = [&](double s){
					return 1. / N_polarisations * q(s) / ( 8 * M_PI * s * s );
				};
				for( auto const& s : s_values ){
					double sum{0.};
					for( auto const& f : _fpolynomials ){
						auto temp = f(s);
						sum += (temp * std::conj(temp)).real();
					}
					result[s] = phase_space(s) * sum;
				}
			}

			return result;
		}

		std::map<double, double> scattering(std::vector<double> const& s_values){
			using namespace std::placeholders;
			int N_polarisations = 1;
			for( auto const& p : _incoming ){
				N_polarisations *= p->spin().n_states();
			}

			std::map<double, double> result;

			std::function<Complex(Complex)> conjugate = [](Complex c){ return std::conj(c);};

			if( N == 1 ){
				auto qout = std::bind(momentum, _1, _outgoing[0]->mass(), _outgoing[1]->mass());
				auto qin  = std::bind(momentum, _1, _incoming[0]->mass(), _incoming[1]->mass());
				auto phase_space = [&](double s){
					return phase_space2(N_polarisations,s, qout(s), qin(s));
				};

				for( auto const& s : s_values ){
					double sum{0.};
					for( auto const& f : _fpolynomials ){
						auto fc = f >> conjugate;
						auto full = fc * f;
						auto integrated = full.integrate(s, -1, 1);
						auto temp = integrated.real();
						sum += integrated.real();
					}
//					std::cout << FORMAT("s: {}\t phi: {}\n", s, phase_space(s));
					result[s] = phase_space(s) * sum;
//					result[s] = sum;
				}
			}

			return result;
		}
	};

	template<>
	Amplitude<0>::Amplitude(std::vector<std::vector<Polynomial>> const& polynomials, std::vector<Particle_Ptr> const& incoming, std::vector<Particle_Ptr> const& virtuals, std::vector<Particle_Ptr> const& outgoing)
	: _incoming(incoming)
	, _virtual(virtuals)
	, _outgoing(outgoing)
	{
		_fpolynomials.reserve(polynomials[0].size());
		for( std::size_t i = 0; i < polynomials[0].size(); ++i ){
			func_t<0> f = [&](){return Complex(1.);};
			_fpolynomials.emplace_back(f, polynomials[0][i]);
		}
	}

	template<>
	Amplitude<1>::Amplitude(std::vector<std::vector<Polynomial>> const& polynomials, std::vector<Particle_Ptr> const& incoming, std::vector<Particle_Ptr> const& virtuals, std::vector<Particle_Ptr> const& outgoing)
	: _incoming(incoming)
    , _virtual(virtuals)
    , _outgoing(outgoing)
	{
		_fpolynomials.reserve(polynomials[0].size());
		for( std::size_t i = 0; i < polynomials[0].size(); ++i ){
			func_t<1> f = polynomials[0][i];
			FPolynomial<1> fp(f, polynomials[1][i]);
			_fpolynomials.push_back(fp);
//			_fpolynomials.emplace_back(f, polynomials[1][i]);
		}
	}
}


#endif