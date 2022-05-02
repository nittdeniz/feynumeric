#include "form_factors.hpp"

#include <cmath>

#include <feynumeric/momentum.hpp>

FORM_FACTOR_FUNCTION identity = [](Feynumeric::Particle_Ptr const&, Feynumeric::Particle_Ptr const&, Feynumeric::Particle_Ptr const&, double){
	return 1.;
};

FORM_FACTOR_FUNCTION moniz = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E){
	return 1.;
};

FORM_FACTOR_FUNCTION manley = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E){
	return 1.;
};

FORM_FACTOR_FUNCTION cassing = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E){
	return 1.;
};

FORM_FACTOR_FUNCTION cutkosky = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E){
	return 1.;
};

FORM_FACTOR_FUNCTION breit_wigner = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E){
	double const a = R->width() * R->width();
	double const b = E - R->mass();
	double const c = R->spin().j();
	return std::pow(a/(a + b*b), c);
};

FORM_FACTOR_FUNCTION dyson_factor_32p =  [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E){
	using namespace Feynumeric;
	double const& m_R  = R->mass();
	double const& m_N  = N->mass();
	double const& m_pi = pi->mass();
	double const q0 = momentum(m_R, m_N, m_pi);
	double const qE = momentum(E, m_N, m_pi);
	auto f = [&](double m){ return (m + m_N - m_pi) * (m + m_N + m_pi);}; // obtained from mathematica
	return std::pow(qE/q0, 3) * f(E) / f(m_R);
};