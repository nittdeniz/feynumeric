#include "form_factors.hpp"

#include <cmath>

#include <feynumeric/momentum.hpp>
#include <feynumeric/units.hpp>

FORM_FACTOR_FUNCTION identity = [](Feynumeric::Particle_Ptr const&, Feynumeric::Particle_Ptr const&, Feynumeric::Particle_Ptr const&, double){
	return 1.;
};

FORM_FACTOR_FUNCTION moniz = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E){
	using namespace Feynumeric;
	using namespace Feynumeric::Units;
	double const beta = 300._MeV;
	double const& m_R  = R->mass();
	double const& m_N  = N->mass();
	double const& m_pi = pi->mass();
	double const q0 = momentum(m_R, m_N, m_pi);
	double const qE = momentum(E, m_N, m_pi);
	double const b2 = beta*beta;
	return std::pow( (b2 + q0 * q0) / (b2 + qE * qE), R->user_data<double>("l") + 1);
};

FORM_FACTOR_FUNCTION manley = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E){
	using namespace Feynumeric;
	using namespace Feynumeric::Units;
	double const beta = 400._MeV;
	double const& m_R  = R->mass();
	double const& m_N  = N->mass();
	double const& m_pi = pi->mass();
	double const q0 = momentum(m_R, m_N, m_pi);
	double const qE = momentum(E, m_N, m_pi);
	double const b2 = beta*beta;
	return std::pow( (b2 + q0 * q0) / (b2 + qE * qE), R->user_data<double>("l"));
};

FORM_FACTOR_FUNCTION cassing = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E){
	using namespace Feynumeric;
	using namespace Feynumeric::Units;
	double const beta = 164._MeV;
	double const& m_R  = R->mass();
	double const& m_N  = N->mass();
	double const& m_pi = pi->mass();
	double const q0 = momentum(m_R, m_N, m_pi);
	double const qE = momentum(E, m_N, m_pi);
	double const b2 = beta*beta;
	return std::pow( (b2 + q0 * q0) / (b2 + qE * qE), (R->user_data<double>("l") + 1)/2.);
};

FORM_FACTOR_FUNCTION cutkosky = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E){
	using namespace Feynumeric;
	using namespace Feynumeric::Units;
	double const& m_R  = R->mass();
	double const& m_N  = N->mass();
	double const& m_pi = pi->mass();
	double const q0 = momentum(m_R, m_N, m_pi);
	double const qE = momentum(E, m_N, m_pi);
	double const m2 = m_pi * m_pi;
	return std::pow( (m_pi + std::sqrt(m2 + q0 * q0) ) / (m_pi + std::sqrt(m2 + qE * qE) ), 2 * R->user_data<double>("l") * 2);
};

FORM_FACTOR_FUNCTION breit_wigner = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double E){
	double const a = 0.3*0.3;
	double const b = E - R->mass();
	double const c = 1.;
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
	return std::pow(qE/q0, 3);// * f(E) / f(m_R);
};