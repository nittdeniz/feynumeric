#include <feynumeric/contract.hpp>
#include <feynumeric/dirac.hpp>
#include <feynumeric/feynman_graph.hpp>
#include <feynumeric/particle_manager.hpp>
#include <feynumeric/qed.hpp>
#include <feynumeric/constexpr_math.hpp>
#include <feynumeric/units.hpp>


#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"
#include <iostream>

using namespace Feynumeric::Units;
using Feynumeric::Particle;
using Feynumeric::Matrix;

Feynumeric::Vertex_Manager_Ptr VMP = std::make_shared<Feynumeric::Vertex_Manager>();

double isospin2_2(std::shared_ptr<Feynumeric::Graph_Edge> a, std::shared_ptr<Feynumeric::Graph_Edge> b)
{
	static Matrix tau(2, 2, {1, std::sqrt(2.), std::sqrt(2.), -1});
	if( a->particle()->isospin().j() != 0.5 || b->particle()->isospin().j() != 0.5 )
	{
		Feynumeric::critical_error(FORMAT("Both {} and {} must have isospin_j = 0.5\n", a->particle()->name(), b->particle()->name()));
	}
	std::size_t in, out;
	if( a->front() == b->back() )
	{
		in = 1-(a->particle()->isospin().m() + 0.5);
		out = 1-(b->particle()->isospin().m() + 0.5);
	}
	else if( a->back() == b->front() )
	{
		in = 1-(b->particle()->isospin().m() + 0.5);
		out = 1-(a->particle()->isospin().m() + 0.5);
	}
	else
	{
		Feynumeric::critical_error("Unsupported isospin for meeting fermions.");
	}
	return tau.at(out, in).real();
}

void init_vertices(Feynumeric::Particle_Manager const& P)
{
	using Feynumeric::Edge_Direction;
	using namespace Feynumeric::Units;

	Feynumeric::QED::init_vertices();
	VMP->import(*Feynumeric::QED::VMP);

	VMP->add(Feynumeric::Vertex(
			{
					{P["N12p"]},
					{P["N"]},
					{P["Pion"]}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& pi = edges[2];
				auto const g = R->particle()->user_data<double>("gRNpi");
				auto const m_pi = pi->particle()->mass();
				auto isospin = isospin2_2(R, N);
				return -g/m_pi * isospin * GA5 * GS(pi->four_momentum(kin));
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{P["N12p"]},
					{P["N"]},
					{Feynumeric::QED::Photon}
			},
			[&](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& gamma = edges[2];
				auto const g = R->particle()->user_data<double>("gRNphoton");
				auto const m_rho = P["rho0"]->mass();
				return -g/m_rho * dirac_sigmac(gamma->four_momentum(kin), gamma->lorentz_indices()[0]);
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{P["N32m"]},
					{P["N"]},
					{P["Pion"]}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& pi = edges[2];
				auto mu = R->lorentz_indices()[0];
				auto const& pR = R->four_momentum(kin);
				auto const& pPi = pi->four_momentum(kin);
				auto const g = R->particle()->user_data<double>("gRNpi");
				auto const m_pi = pi->particle()->mass();
				auto const isospin = isospin2_2(R, N);
				Lorentz_Index_Ptr nu = std::make_shared<Lorentz_Index>();
				auto const check0 = O32c(pR, nu, mu) * pPi.contra(nu);
				++(*nu);
				auto const check1 = O32c(pR, nu, mu) * pPi.contra(nu);
				++(*nu);
				auto const check2 = O32c(pR, nu, mu) * pPi.contra(nu);
				++(*nu);
				auto const check3 = O32c(pR, nu, mu) * pPi.contra(nu);
				auto result = CONTRACT_MATRIX(O32c(pR, lambda, mu) * pPi.contra(lambda), lambda) * GA5;
				double beta = 300._MeV;
				double q0 = 1.;
                double q = 2.;
				int l = R->particle()->user_data<double>("l");
				auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(beta, q0, q, l);
//				std::cerr << FORMAT("pion: {} {} {} {}\n", isospin * g/(m_pi*m_pi),isospin, g, m_pi);
				return 1.i * isospin * g/(m_pi*m_pi) *
					CONTRACT_MATRIX(O32c(pR, lambda, mu) * pPi.contra(lambda), lambda) * GA5;

			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{P["N32m"], Edge_Direction::IN},
					{P["N"], Edge_Direction::OUT},
					{Feynumeric::QED::Photon}
			},
			[&](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& photon = edges[2];
				auto kappa = photon->lorentz_indices()[0];
				auto lambda = R->lorentz_indices()[0];
				auto const& pR = R->four_momentum(kin);
				auto const& pPhoton = photon->four_momentum(kin);
				auto const g = R->particle()->user_data<double>("gRNphoton");
				auto const m_rho = P["rho0"]->mass();
//				std::cerr << FORMAT("photon: {} {} {}\n", g/(4*m_rho*m_rho), g, m_rho);
				auto result = g/(4*m_rho*m_rho) * CONTRACT_MATRIX(O32c(pPhoton, mu, kappa) * O32c(pR, mu, lambda) * MT[*mu][*mu], mu);
				return result;
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{P["N32m"], Edge_Direction::OUT},
					{P["N"], Edge_Direction::IN},
					{Feynumeric::QED::Photon}
			},
			[&](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& photon = edges[2];
				auto kappa = photon->lorentz_indices()[0];
				auto lambda = R->lorentz_indices()[0];
				auto const& pR = R->four_momentum(kin);
				auto const& pPhoton = photon->four_momentum(kin);
				auto const g = R->particle()->user_data<double>("gRNphoton");
//				auto const g = std::any_cast<double>(R->particle()->user_data("gRNpi"));
				auto const m_rho = P["rho0"]->mass();
//				std::cerr << FORMAT("photon: {} {} {}\n", g/(4*m_rho*m_rho), g, m_rho);
				auto result = g/(4*m_rho*m_rho) *  CONTRACT_MATRIX(O32c(pR, mu, lambda) * O32c(pPhoton, mu, kappa) * MT[*mu][*mu], mu);
				return result;

			}
	));
}