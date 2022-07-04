#include "form_factors.hpp"

#include <cmath>

#include <feynumeric/momentum.hpp>
#include <feynumeric/units.hpp>

namespace Form_Factor
{

	std::string const CMD_FORM_FACTOR_NONE = "none";
	std::string const CMD_FORM_FACTOR_CASSING = "cassing";
	std::string const CMD_FORM_FACTOR_CUTKOSKY = "cutkosky";
	std::string const CMD_FORM_FACTOR_MANLEY = "manley";
	std::string const CMD_FORM_FACTOR_MONIZ = "moniz";
	std::string const CMD_FORM_FACTOR_BREIT_WIGNER = "breit_wigner";
	std::string const CMD_FORM_FACTOR_GAUSSIAN = "gaussian";
	std::string const CMD_FORM_FACTOR_MULTIPOL_GAUSS = "multipolgauss";
	std::string const CMD_FORM_FACTOR_RAYLEIGH = "rayleigh";
	std::string const CMD_FORM_FACTOR_INVERSE_GAUSS = "inverse_gaussian";

	FORM_FACTOR_FUNCTION identity = [](Feynumeric::Particle_Ptr const&, Feynumeric::Particle_Ptr const&,
	                                   Feynumeric::Particle_Ptr const&, double){
		return 1.;
	};

	FORM_FACTOR_FUNCTION moniz = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N,
	                                Feynumeric::Particle_Ptr const& pi, double E){
		using namespace Feynumeric;
		using namespace Feynumeric::Units;
		double const beta = 300._MeV;
		double const& m_R = R->mass();
		double const& m_N = N->mass();
		double const& m_pi = pi->mass();
		double const q0 = momentum(m_R, m_N, m_pi);
		double const qE = momentum(E, m_N, m_pi);
		double const b2 = beta * beta;
		return std::pow(( b2 + q0 * q0 ) / ( b2 + qE * qE ), R->user_data<double>("l") + 1);
	};

	FORM_FACTOR_FUNCTION manley = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N,
	                                 Feynumeric::Particle_Ptr const& pi, double E){
		using namespace Feynumeric;
		using namespace Feynumeric::Units;
		double const beta = 400._MeV;
		double const& m_R = R->mass();
		double const& m_N = N->mass();
		double const& m_pi = pi->mass();
		double const q0 = momentum(m_R, m_N, m_pi);
		double const qE = momentum(E, m_N, m_pi);
		double const b2 = beta * beta;
		return std::pow(( b2 + q0 * q0 ) / ( b2 + qE * qE ), R->user_data<double>("l"));
	};

	FORM_FACTOR_FUNCTION cassing = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N,
	                                  Feynumeric::Particle_Ptr const& pi, double E){
		using namespace Feynumeric;
		using namespace Feynumeric::Units;
		double const beta = 164._MeV;
		double const& m_R = R->mass();
		double const& m_N = N->mass();
		double const& m_pi = pi->mass();
		double const q0 = momentum(m_R, m_N, m_pi);
		double const qE = momentum(E, m_N, m_pi);
		double const b2 = beta * beta;
		return std::pow(( b2 + q0 * q0 ) / ( b2 + qE * qE ), ( R->user_data<double>("l") + 1 ) / 2.);
	};

	FORM_FACTOR_FUNCTION cutkosky = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N,
	                                   Feynumeric::Particle_Ptr const& pi, double E){
		using namespace Feynumeric;
		using namespace Feynumeric::Units;
		double const& m_R = R->mass();
		double const& m_N = N->mass();
		double const& m_pi = pi->mass();
		double const q0 = momentum(m_R, m_N, m_pi);
		double const qE = momentum(E, m_N, m_pi);
		double const m2 = m_pi * m_pi;
		return std::pow(( m_pi + std::sqrt(m2 + q0 * q0)) / ( m_pi + std::sqrt(m2 + qE * qE)),
		                2 * R->user_data<double>("l") * 2);
	};

	FORM_FACTOR_FUNCTION breit_wigner = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N,
	                                       Feynumeric::Particle_Ptr const& pi, double E){
		double const l = 0.8;//4*R->spin().j()*R->width();
		double const l4 = std::pow(l, 4);
		double const b = E - R->mass();
		double const c = 1.;
		return l4 / ( std::pow(E * E - R->mass() * R->mass(), 2) + l4 );
	};

	FORM_FACTOR_FUNCTION dyson_factor_32 = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N,
	                                          Feynumeric::Particle_Ptr const& pi, double E){
		using namespace Feynumeric;
		double const& m_R = R->mass();
		double const& m_N = N->mass();
		double const& m_pi = pi->mass();
		double const q0 = momentum(m_R, m_N, m_pi);
		double const qE = momentum(E, m_N, m_pi);
		auto f = [](double a, double b){ return std::sqrt(a * a + b * b); };
		return std::pow(qE / q0, 3) * E / m_R * ( m_N + f(m_N, qE)) / ( m_N + f(m_N, q0));
	};

	FORM_FACTOR_FUNCTION dyson_factor_12 = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N,
	                                          Feynumeric::Particle_Ptr const& pi, double E){
		using namespace Feynumeric;
		double const& m_R = R->mass();
		double const& m_N = N->mass();
		double const& m_pi = pi->mass();
		double const q0 = momentum(m_R, m_N, m_pi);
		double const qE = momentum(E, m_N, m_pi);
		auto f = [](double a, double b){ return std::sqrt(a * a + b * b); };
		auto g = [&](double q){
			double mp2 = m_pi * m_pi;
			double m2 = m_N * m_N;
			return q * ( m_N * mp2 + mp2 * f(m_N, q) + 2 * q * q * ( f(m_N, q) + f(m_pi, q)));
		};
		return m_R / E * g(qE) / g(q0);
	};

	FORM_FACTOR_FUNCTION gaussian = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N,
	                                   Feynumeric::Particle_Ptr const& pi, double E){
		return std::exp(-std::pow(E * E - R->mass() * R->mass(), 2) / ( std::pow(0.7, 4)));
	};

	FORM_FACTOR_FUNCTION rayleigh = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N,
	                                   Feynumeric::Particle_Ptr const& pi, double E){
		auto l = 0.61;
		auto l4 = std::pow(l, 4);
		auto kin_limit = 0.95 * ( N->mass() + pi->mass());
		return ( R->mass() - kin_limit ) / ( E - kin_limit ) *
		       std::exp(-std::pow(E * E - R->mass() * R->mass(), 2) / l4);
	};

	FORM_FACTOR_FUNCTION inverse_gaussian = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N,
	                                           Feynumeric::Particle_Ptr const& pi, double E){
		auto l = 10.;
		return std::pow(E / R->mass(), -1.5) * std::exp(-l * std::pow(E - R->mass(), 2) / ( R->mass() * E ));
	};

	FORM_FACTOR_FUNCTION multipol_gauss = [](Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N,
	                                         Feynumeric::Particle_Ptr const& pi, double E){
		auto Gamma_tilde = R->width() / std::sqrt(std::pow(2, 1. / ( 2 * R->spin().j() * 2 ) - 1));
		auto M2 = R->mass() * R->mass();
		auto gamma_M = Gamma_tilde * Gamma_tilde * M2;
		auto s = E * E;
		auto Lambda4 = 1;
		auto multipol = std::pow(gamma_M / ( std::pow(s - M2, 2) + gamma_M ), R->spin().j() - 0.5);
		auto gauss = std::exp(-std::pow(s - M2, 2) / Lambda4);
		return std::sqrt(std::sqrt(2)) * multipol * gauss;
	};

	std::map<std::string, FORM_FACTOR_FUNCTION> ff_dict = {
			{CMD_FORM_FACTOR_NONE,           identity},
			{CMD_FORM_FACTOR_MULTIPOL_GAUSS, multipol_gauss},
			{CMD_FORM_FACTOR_MONIZ,          moniz},
			{CMD_FORM_FACTOR_CUTKOSKY,       cutkosky},
			{CMD_FORM_FACTOR_CASSING,        cassing},
			{CMD_FORM_FACTOR_BREIT_WIGNER,   breit_wigner},
			{CMD_FORM_FACTOR_GAUSSIAN,       gaussian},
			{CMD_FORM_FACTOR_MANLEY,         manley},
			{CMD_FORM_FACTOR_RAYLEIGH,       rayleigh},
			{CMD_FORM_FACTOR_INVERSE_GAUSS,  inverse_gaussian}
	};
}