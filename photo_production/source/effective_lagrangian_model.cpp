#include <feynumeric/contract.hpp>
#include <feynumeric/dirac.hpp>
#include <feynumeric/feynman_graph.hpp>
#include <feynumeric/isospin.hpp>
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

void init_vertices(Feynumeric::Particle_Manager const& P)
{
	using Feynumeric::Edge_Direction;
	using namespace Feynumeric::Units;

//	Feynumeric::QED::init_vertices();
//	VMP->import(*Feynumeric::QED::VMP);

	VMP->add(Feynumeric::Vertex(
			{
					{P["R12"]},
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
				auto isospin = R->particle()->isospin().j() == 1.5 ? isospin_T(R, N) : isospin_tau(R, N);
				return -g/m_pi * isospin * GA5 * GS(pi->four_momentum(kin));
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{P["R12"]},
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
					{P["R32"]},
					{P["N"]},
					{P["Pion"]}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				static auto const Identity = Matrix(4, 4, 1);
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& pi = edges[2];
				auto mu = R->lorentz_indices()[0];
				auto const& pR = R->four_momentum(kin);
				auto const& pPi = pi->four_momentum(kin);
				auto const g = R->particle()->user_data<double>("gRNpi");
				auto const m_pi = pi->particle()->mass();
				auto const isospin = R->particle()->isospin().j() == 1.5 ? isospin_T(R, N) : isospin_tau(R, N);
				auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), pi->particle(), pR.E().real());
				auto const ga5 = (R->particle()->parity() == -1 ? GA5 : Identity);
				return -1.i * form_factor * isospin * g/(m_pi*m_pi) * O32c(pR, pPi, mu) * ga5;
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{P["R32"], Edge_Direction::IN},
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
				auto result = 1.i * g/(4*m_rho*m_rho) * CONTRACT_MATRIX(O32c(pPhoton, mu, kappa) * O32c(pR, mu, lambda) * MT[*mu][*mu], mu);
				return result;
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{P["R32"], Edge_Direction::OUT},
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
				auto const m_rho = P["rho0"]->mass();
				auto result = 1.i * g/(4*m_rho*m_rho) *  CONTRACT_MATRIX(O32c(pR, mu, lambda) * O32c(pPhoton, mu, kappa) * MT[*mu][*mu], mu);
				return result;
			}
	));
}