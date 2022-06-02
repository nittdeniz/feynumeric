#include <feynumeric/contract.hpp>
#include <feynumeric/dirac.hpp>
#include <feynumeric/feynman_graph.hpp>
#include <feynumeric/graph_edge.hpp>
#include <feynumeric/graph_vertex.hpp>
#include <feynumeric/isospin.hpp>
#include <feynumeric/messages.hpp>
#include <feynumeric/particle_manager.hpp>
#include <feynumeric/qed.hpp>
#include <feynumeric/constexpr_math.hpp>

#include <fstream>
#include <iostream>
#include <sstream>


#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"
#include "coupling_constants.hpp"

using Feynumeric::Particle;
using Feynumeric::Matrix;

Feynumeric::Vertex_Manager_Ptr VMP = std::make_shared<Feynumeric::Vertex_Manager>();

Feynumeric::Matrix gamma5(int lorentz_parity, int particle_parity){
	return lorentz_parity == particle_parity ? Matrix(1,1,1) : Feynumeric::GA5;
}

void init_vertices(Feynumeric::Particle_Manager const& P)
{
	using Feynumeric::Edge_Direction;

	Feynumeric::QED::init_vertices();
	VMP->import(*Feynumeric::QED::VMP);

	std::map<std::string, double> coupling_constants;
	std::string const file_name = "data/coupling_constants_isospin_symmetry.txt";
	std::ifstream ifs(file_name);
	if( !ifs ){
	    Feynumeric::critical_error(FORMAT("Could not open {}.", file_name));
	}

	std::string buffer;
	while( std::getline(ifs, buffer) ){
        std::string key;
        double value;
        std::stringstream sstream(buffer);
        sstream >> key >> value;
        coupling_constants[key] = value;
        std::cout << key << "/" << value << "\n";
	}

	std::cout << coupling_constants.size() << "\n";

	VMP->add(Feynumeric::Vertex(
			{
					{P["N"]},
					{P["N"]},
					{P["Pion"]}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& pi = edges[2];
				auto const g = 0.97;
				auto const m_pi = pi->particle()->mass();
				auto const iso = (R->back() == N->front())? isospin(R, N, pi) : isospin(N, R, pi);
				return g/m_pi * iso * GA5 * GS(pi->four_momentum(kin));
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{P["N"]},
					{P["N"]},
					{P["Rho"]}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				auto const& N1 = edges[0];
				auto const& N2 = edges[1];
				auto const& rho = edges[2];
				auto mu = rho->lorentz_indices()[0];
				auto const g = 5.96;
				auto const iso = (N1->back() == N2->front())? isospin(N1, N2, rho) : isospin(N2, N1, rho);
				return g/2. * iso * GAC[*mu];
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{P["Pion"]},
					{P["Pion"]},
					{P["Rho"]}
			},
			[coupling_constants](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				auto const& piP = edges[0];
				auto const& piM = edges[1];
				auto const& rho = edges[2];
				auto const& p = piP->particle()->charge() >= 0? piP->four_momentum(kin) : piM->four_momentum(kin);
				auto const& q = piP->particle()->charge() >= 0? piM->four_momentum(kin) : piP->four_momentum(kin);
				auto mu = rho->lorentz_indices()[0];
				//auto const g = rho->particle()->user_data<double>("g_pipi");
				auto const g = coupling_constants.at("grhopipi");
				//int l = static_cast<int>(R->particle()->user_data<double>("l"));
				//auto isospin = R->particle()->isospin().j() == 1.5 ? isospin_T(R, N) : isospin_tau(R, N);
				return g * (p-q).co(mu);
			}
	));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["Pion"], Edge_Direction::ANY},
                    {P["Pion"], Edge_Direction::ANY},
                    {P["f0_500"], Edge_Direction::ANY}
            },
            [](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                auto const& piP = edges[0];
                auto const& piM = edges[1];
                auto const& f0  = edges[2];
                auto const& p = f0->four_momentum(kin);
                auto const m = piP->particle()->mass();
                auto const g = f0->particle()->user_data<double>("g_pipi");
                return Matrix(1, 1, g/m * p.squared());
            }
    ));

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
				auto const iso = (R->back() == N->front())? isospin(R, N, pi) : isospin(N, R, pi);
				int const lorentz_parity = 1;
				int const particle_parity = R->particle()->parity() * N->particle()->parity() * pi->particle()->parity();
				//auto isospin = R->particle()->isospin().j() == 1.5 ? isospin_T(R, N) : isospin_tau(R, N);
				return -g / m_pi * iso * gamma5(lorentz_parity, particle_parity) * GS(pi->four_momentum(kin));
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{P.get("R12")},
					{P.get("N")},
					{P.get("f0_500")}
			},
			[&](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& f0 = edges[2];
				auto const g = R->particle()->user_data<double>("gRNf0_500");
				auto const m_pi = P.get("pi0")->mass();
				int const lorentz_parity = 1;
				int const particle_parity = R->particle()->parity() * N->particle()->parity() * f0->particle()->parity();
				return g/m_pi * gamma5(lorentz_parity, particle_parity) * GS(f0->four_momentum(kin));
			}
	));

    VMP->add(Feynumeric::Vertex(
            {
                    {P.get("R12")},
                    {P.get("N")},
                    {P.get("eta")}
            },
            [&](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                auto const& R = edges[0];
                auto const& N = edges[1];
                auto const& eta = edges[2];
                auto const g = R->particle()->user_data<double>("gRNeta");
                auto const m_pi = P.get("pi0")->mass();
                int const lorentz_parity = 1;
                int const particle_parity = R->particle()->parity() * N->particle()->parity() * eta->particle()->parity();
                return g/m_pi * gamma5(lorentz_parity, particle_parity) * GS(eta->four_momentum(kin));
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
				//auto const isospin = R->particle()->isospin().j() == 1.5 ? isospin_T(R, N) : isospin_tau(R, N);
				auto const iso = (R->back() == N->front())? isospin(R, N, pi) : isospin(N, R, pi);
				auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), pi->particle(), pR.E().real());
				int l = static_cast<int>(R->particle()->user_data<double>("l"));
				Complex phase = std::exp(2.i * M_PI/360. * R->particle()->user_data<double>("phase"));
				return -1.i * phase * form_factor * iso * g/(m_pi*m_pi) * O32c(pR, pPi, mu) * gamma5(l, R->particle()->parity());
			}
	));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["R32"]},
                    {P["N"]},
                    {P["eta"]}
            },
            [](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                static auto const Identity = Matrix(4, 4, 1);
                auto const& R = edges[0];
                auto const& N = edges[1];
                auto const& eta = edges[2];
                auto mu = R->lorentz_indices()[0];
                auto const& pR = R->four_momentum(kin);
                auto const& pPi = eta->four_momentum(kin);
                auto const g = R->particle()->user_data<double>("gRNeta");
                auto const m_pi = eta->particle()->mass();
                //auto const isospin = R->particle()->isospin().j() == 1.5 ? isospin_T(R, N) : isospin_tau(R, N);
                //auto const iso = (R->back() == N->front()) ? isospin(R, N, eta) : isospin(N, R, eta);
                //auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), eta->particle(), pR.E().real());
                int l = static_cast<int>(R->particle()->user_data<double>("l"));
                Complex phase = std::exp(2.i * M_PI/360. * R->particle()->user_data<double>("phase"));
                return -1.i * phase * g/(m_pi*m_pi) * O32c(pR, pPi, mu) * gamma5(l, R->particle()->parity());
            }
    ));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["R32"]},
                    {P["R12"]},
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
                auto const& pN = N->four_momentum(kin);
                auto const& pPi = pi->four_momentum(kin);
                std::array<std::string, 2> names{R->particle()->name().substr(0, 5), N->particle()->name().substr(0, 5)};
                std::sort(names.begin(), names.end());
                auto const g = R->particle()->user_data<double>(FORMAT("g{}{}pi", names[0], names[1]));
                auto const m_pi = pi->particle()->mass();
                //auto const isospin = R->particle()->isospin().j() == 1.5 ? isospin_T(R, N) : isospin_tau(R, N);
                auto const iso = (R->back() == N->front())? isospin(R, N, pi) : isospin(N, R, pi);
                auto form_factor1 = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), pi->particle(), pR.E().real());
                auto form_factor2 = N->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(N->particle(), R->particle(), pi->particle(), pN.E().real());
                int l = 1;
                Complex phase = std::exp(2.i * M_PI/360. * R->particle()->user_data<double>("phase"));
                return -1.i * phase * form_factor1 * form_factor2 * iso * g/(m_pi*m_pi) * O32c(pR, pPi, mu) * gamma5(l, R->particle()->parity());
            }
    ));

	VMP->add(Feynumeric::Vertex(
			{
					{P.get("R32"), Edge_Direction::OUT},
					{P.get("R32"), Edge_Direction::IN},
					{P.get("Pion")}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				static auto const Identity = Matrix(4, 4, 1);
				auto const& Rout = edges[0];
				auto const& Rin = edges[1];
				auto const& pi = edges[2];
				auto mu = Rout->lorentz_indices()[0];
				auto nu = Rin->lorentz_indices()[0];
				auto const& pRout = Rout->four_momentum(kin);
				auto const& pRin = Rin->four_momentum(kin);
				auto const& pPi = pi->four_momentum(kin);

				auto const g = Rout->particle()->user_data<double>("gD1232N1520pi");
				auto const m_pi = pi->particle()->mass();
				//auto const isospin = R->particle()->isospin().j() == 1.5 ? isospin_T(R, N) : isospin_tau(R, N);
				auto const iso = isospin(Rout, Rin, pi);
				//auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), pi->particle(), pR.E().real());
				int const lorentz_parity = 1;
				//Complex phase = std::exp(2.i * M_PI/360. * R->particle()->user_data<double>("phase"));
				int const particle_parity = Rout->particle()->parity() * Rin->particle()->parity() * pi->particle()->parity();
				return -1.i * iso * g/(m_pi*m_pi*m_pi*m_pi ) * O32c(pRout, pPi, mu) * gamma5(lorentz_parity, particle_parity) * O32c(pRin, pPi, nu);
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


	VMP->add(Feynumeric::Vertex(
			{
					{P.get("R32"), Edge_Direction::IN},
					{P.get("N"), Edge_Direction::OUT},
					{P.get("Rho")}
			},
			[&](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& rho = edges[2];
				auto kappa = rho->lorentz_indices()[0];
				auto lambda = R->lorentz_indices()[0];
				auto const& pR = R->four_momentum(kin);
				auto const& pRho = rho->four_momentum(kin);
				auto const g = R->particle()->user_data<double>("gRNphoton");
				auto const m_rho = rho->particle()->mass();
				auto result = 1.i * g/(4*m_rho*m_rho) * CONTRACT_MATRIX(O32c(pRho, mu, kappa) * O32c(pR, mu, lambda) * MT[*mu][*mu], mu);
				return result;
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{P["R32"], Edge_Direction::OUT},
					{P["N"], Edge_Direction::IN},
					{P.get("Rho")}
			},
			[&](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& Rho = edges[2];
				auto kappa = Rho->lorentz_indices()[0];
				auto lambda = R->lorentz_indices()[0];
				auto const& pR = R->four_momentum(kin);
				auto const& pRho = Rho->four_momentum(kin);
				auto const g = R->particle()->user_data<double>("gRNphoton");
				auto const m_rho = Rho->particle()->mass();
				auto result = 1.i * g/(4*m_rho*m_rho) *  CONTRACT_MATRIX(O32c(pR, mu, lambda) * O32c(pRho, mu, kappa) * MT[*mu][*mu], mu);
				return result;
			}
	));

}