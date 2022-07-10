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
		std::vector<Particle_Ptr> _outgoing;
		int n_pol = 1;

	public:
		Amplitude() = default;
		Amplitude(Amplitude<N> const& other )
		: _fpolynomials(other._fpolynomials)
		, _incoming(other._incoming)
		, _outgoing(other._outgoing)
		{
			for( auto const& p : _incoming ){
				n_pol *= p->spin().n_states();
			}
		}


		Amplitude<N>& operator=(Amplitude<N> const& other){
			_fpolynomials = other._fpolynomials;
			_incoming = other._incoming;
			_outgoing = other._outgoing;
			n_pol = other.n_pol;
			return *this;
		}

		Amplitude(std::vector<std::vector<Polynomial>> const& polynomials, std::vector<Particle_Ptr> const& incoming, std::vector<Particle_Ptr> const& outgoing);

		void scale(func_t<N> const& func){
			for( auto& fp : _fpolynomials ){
				for( auto& f : fp._coefficients ){
					f = f * func;
				}
			}
		}

		Amplitude<N> conjugate() const{

		}

		double width(double s){
			using namespace std::placeholders;
			if( N == 0 ){
				auto const q = momentum(s, _outgoing[0]->mass(), _outgoing[1]->mass());
				auto phase_space = q/n_pol/(8. * M_PI * s * s);
				double sum{0.};
				for( auto const& f : _fpolynomials ){
					auto temp = f(s);
					sum += (temp * std::conj(temp)).real();
				}
				return phase_space * sum;
			}
		}

		func_t<1> scattering(){
			using namespace std::placeholders;
			std::function<Complex(Complex)> conjugate = [](Complex c){ return std::conj(c); };

			if( N == 1 ){
				auto qout = std::bind(momentum, _1, _outgoing[0]->mass(), _outgoing[1]->mass());
				auto qin = std::bind(momentum, _1, _incoming[0]->mass(), _incoming[1]->mass());
				func_t<1> phase_space = [=, this](double s) mutable{
					return phase_space2(n_pol, s, qout(s), qin(s));
				};
				func_t<1> sum = [](double){return 0.;};
				for( auto const& f : _fpolynomials ){
					auto fc = f >> conjugate;
					auto full = fc * f;
					auto integrated = full.integrate(-1, 1);
					sum = sum + integrated;
				}
				return phase_space * sum;
			}
			critical_error("scattering only supports 2 final state particles (i.e. N=1)");
		}

		std::map<double, double> scattering(std::vector<double> const& s_values){
			using namespace std::placeholders;
			int N_polarisations = 1;
			for( auto const& p : _incoming ){
				N_polarisations *= p->spin().n_states();
			}

			std::map<double, double> result;

			std::function<Complex(Complex)> conjugate = [](Complex c){ return std::conj(c); };

			if( N == 1 ){
				auto qout = std::bind(momentum, _1, _outgoing[0]->mass(), _outgoing[1]->mass());
				auto qin = std::bind(momentum, _1, _incoming[0]->mass(), _incoming[1]->mass());
				auto phase_space = [=](double s) mutable{
					return phase_space2(N_polarisations, s, qout(s), qin(s));
				};

				#pragma omp parallel for
				for( std::size_t i = 0; i < s_values.size(); ++i ){
					auto const s = s_values[i];
					double sum{0.};
					for( auto const& f : _fpolynomials ){
						auto fc = f >> conjugate;
//						std::cout << "f: " << f(1.2, 0.5) << "\n";
//						std::cout << "fc: " << fc(1.2, 0.5) << "\n";
						auto full = fc * f;
//						std::cout << "full: " << full(1.2, 0.5) << "\n";
						auto integrated = full.integrate(s, -1, 1);
//						auto i12 = full.integrate(1.2, -1, 1);
//						std::cout << "int: " << integrated << "\n";
						auto temp = integrated.real();
						sum += integrated.real();
					}
//					std::cout << FORMAT("s: {}\t phi: {}\n", s, phase_space(s));
					#pragma omp critical
					{
						result[s] = phase_space(s) * sum;
					}
//					result[s] = sum;
				}
			}
			return result;
		}

		template <std::size_t U>
		friend Amplitude<U> operator+(Amplitude<U> const& lhs, Amplitude<U> const& rhs);
	};

	template<>
	Amplitude<0>::Amplitude(std::vector<std::vector<Polynomial>> const& polynomials, std::vector<Particle_Ptr> const& incoming, std::vector<Particle_Ptr> const& outgoing)
	: _incoming(incoming)
	, _outgoing(outgoing)
	{
		for( auto& p : incoming ){
			n_pol *= p->spin().n_states();
		}
		_fpolynomials.reserve(polynomials[0].size());
		for( std::size_t i = 0; i < polynomials[0].size(); ++i ){
			func_t<0> f = [&](){return Complex(1.);};
			_fpolynomials.emplace_back(f, polynomials[0][i]);
		}
	}

	template<>
	Amplitude<1>::Amplitude(std::vector<std::vector<Polynomial>> const& polynomials, std::vector<Particle_Ptr> const& incoming, std::vector<Particle_Ptr> const& outgoing)
	: _incoming(incoming)
    , _outgoing(outgoing)
	{
		for( auto& p : incoming ){
			n_pol *= p->spin().n_states();
		}
		_fpolynomials.reserve(polynomials[0].size());
		for( std::size_t i = 0; i < polynomials[0].size(); ++i ){
			func_t<1> f = polynomials[0][i];
			FPolynomial<1> fp(f, polynomials[1][i]);
			_fpolynomials.emplace_back(f, polynomials[1][i]);
		}
	}

	template <std::size_t U>
	inline Amplitude<U> operator+(Amplitude<U> const& lhs, Amplitude<U> const& rhs){
		if( lhs._incoming.size() != rhs._incoming.size() ){
			critical_error("Can not add amplitudes with unequal number of incoming particles.");
		}
		if( lhs._outgoing.size() != rhs._outgoing.size() ){
			critical_error("Can not add amplitudes with unequal number of outgoing particles.");
		}
		for( std::size_t i = 0; i < lhs._incoming.size(); ++i ){
			if( lhs._incoming[i] != rhs._incoming[i] ){
				critical_error("Can not add amplitudes with different incoming particles.");
			}
		}
		for( std::size_t i = 0; i < lhs._outgoing.size(); ++i ){
			if( lhs._outgoing[i] != rhs._outgoing[i] ){
				critical_error("Can not add amplitudes with different outgoing particles.");
			}
		}

		Amplitude<U> result;
		result._incoming = lhs._incoming;
		result._outgoing = rhs._outgoing;
		result.n_pol = lhs.n_pol;
		result._fpolynomials.resize(lhs._fpolynomials.size());
		for( std::size_t i = 0; i < result._fpolynomials.size(); ++i ){
			auto m = lhs._fpolynomials[i]._coefficients.size();
			auto n = rhs._fpolynomials[i]._coefficients.size();
			auto max = std::max(m, n);
			result._fpolynomials[i]._coefficients.resize(max);
			result._fpolynomials[i].n = max;
			result._fpolynomials[i].order = max-1;
			for( std::size_t j = 0; j < max; ++j ){
				if( j < m && j < n ){
					result._fpolynomials[i]._coefficients[j] = lhs._fpolynomials[i]._coefficients[j] + rhs._fpolynomials[i]._coefficients[j];
				}else if( j < m ){
					result._fpolynomials[i]._coefficients[j] = lhs._fpolynomials[i]._coefficients[j];
				}else if( j < n ){
					result._fpolynomials[i]._coefficients[j] = rhs._fpolynomials[i]._coefficients[j];
				}else{
					critical_error("This should not happen.");
				}
			}
		}
		return result;
	}
}


#endif