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
#include <feynumeric/units.hpp>

#include <fstream>
#include <iostream>
#include <sstream>

#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"

using Feynumeric::Particle;
using Feynumeric::Matrix;

Feynumeric::Vertex_Manager_Ptr VMP = std::make_shared<Feynumeric::Vertex_Manager>();

Couplings couplings;

Feynumeric::Matrix gamma5(int lorentz_parity, int particle_parity){
	return lorentz_parity == particle_parity ? Matrix(1,1,1) : Feynumeric::GA5;
}

void init_vertices(Feynumeric::Particle_Manager const& P, std::string const& coupling_constant_file)
{
	using Feynumeric::Edge_Direction;
    couplings.load(coupling_constant_file);
	Feynumeric::QED::init_vertices();
	VMP->import(*Feynumeric::QED::VMP);

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
				return -1. * g/m_pi * iso * GA5 * GS(pi->four_momentum(kin));
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{P["N"], Edge_Direction::OUT},
					{P["N"], Edge_Direction::IN},
					{P["Rho"]}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				auto const& N1 = edges[0];
				auto const& N2 = edges[1];
				auto const& rho = edges[2];
				auto mu = rho->lorentz_indices()[0];
				auto const g = 5.96;
                auto const kappa = 0;// 3.71;
				auto const iso = isospin(N1, N2, rho);//(N1->back() == N2->front())? isospin(N1, N2, rho) : isospin(N2, N1, rho);
				return 1.i * g/2. * iso * (GAC[*mu] + kappa/(.1*rho->particle()->mass()) * dirac_sigma(mu, N1->four_momentum(kin)));
			}
	));

    VMP->add(Feynumeric::Vertex(
            {
                    {Feynumeric::QED::Photon},
                    {P["Rho"]}
            },
            [](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                using namespace Feynumeric::Units;
                auto const& photon = edges[0];
                auto const& rho = edges[1];
                auto mu = photon->lorentz_indices()[0];
                auto nu = rho->lorentz_indices()[0];
                auto q2 = photon->four_momentum(kin).squared();
                auto const coupl_str = coupling_string("Rho", "Pion", "Pion");
                auto const g = couplings.get(coupl_str);
                double e = 1._e;
                return Matrix(1, 1, -1.i * e * q2/g * MT[*mu][*nu]);
            }
    ));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["N"]},
                    {P["N"]},
                    {Feynumeric::QED::Photon}
            },
            [](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                using namespace Feynumeric::Units;
                double const g = 1._e;
                auto const& photon = edges[2];
                auto const& mu = photon->lorentz_indices()[0];
                return g * GAC[*mu];
            }
    ));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["Pion"]},
                    {P["Pion"]},
                    {Feynumeric::QED::Photon}
            },
            [](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                using namespace Feynumeric::Units;
                auto const& piP = edges[0];
                auto const& piM = edges[1];
                auto const& photon = edges[2];
                auto const& p = piP->particle()->charge() >= 0? piP->four_momentum(kin) : piM->four_momentum(kin);
                auto const& q = piP->particle()->charge() >= 0? piM->four_momentum(kin) : piP->four_momentum(kin);
                auto mu = photon->lorentz_indices()[0];
                double const g = 1._e;
                return 1. * g * (p-q).co(mu);
            }
    ));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["N"]},
                    {P["N"]},
                    {P["f0_500"]}
            },
            [](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                using namespace Feynumeric::Units;
                auto const coupl_str = coupling_string("f0_500", "N", "N");
                auto const g = couplings.get(coupl_str);
//                auto const& N1 = edges[0];
//                auto const& N2 = edges[1];
//                auto const f0 = edges[2]->particle();
//                auto const m2 = f0->mass() * f0->mass();
//                auto const p1 = N1->four_momentum(kin);
//                auto const p2 = N2->four_momentum(kin);
//                auto const q2 = (p1-p2).squared();
//                auto const Lambda = 2._GeV;
//                auto const Lambda2 = Lambda*Lambda;
//                auto const f = static_cast<double>((Lambda2 - m2)/(Lambda2 - q2));
                return Matrix(1, 1, -1.i * g);
            }
    ));

	VMP->add(Feynumeric::Vertex(
			{
					{P["Pion"]},
					{P["Pion"]},
					{P["Rho"]}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
				using namespace Feynumeric;
				auto const& piP = edges[0];
				auto const& piM = edges[1];
				auto const& rho = edges[2];
				auto const& p = piP->particle()->charge() >= 0? piP->four_momentum(kin) : piM->four_momentum(kin);
				auto const& q = piP->particle()->charge() >= 0? piM->four_momentum(kin) : piP->four_momentum(kin);
				auto mu = rho->lorentz_indices()[0];
				//auto const g = rho->particle()->user_data<double>("g_pipi");
                auto const coupl_str = coupling_string("Rho", "Pion", "Pion");
                auto const g = couplings.get(coupl_str);
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
                //auto const g = f0->particle()->user_data<double>("g_pipi");
                auto const coupl_str = coupling_string("f0_500", "Pion", "Pion");
                auto const g = couplings.get(coupl_str);
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
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "Pion");
				auto const g = couplings.get(coupl_str);
				auto const m_pi = pi->particle()->mass();
				auto const iso = (R->back() == N->front())? isospin(R, N, pi) : isospin(N, R, pi);
				int const lorentz_parity = 1;
				int const particle_parity = R->particle()->parity() * N->particle()->parity() * pi->particle()->parity();
				auto const pR = R->four_momentum(kin);
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), pi->particle(), pR.E().real());
//                std::cout << pi->four_momentum(kin) << "\n";
				//auto isospin = R->particle()->isospin().j() == 1.5 ? isospin_T(R, N) : isospin_tau(R, N);
				return 1.i * g / m_pi * iso * form_factor * gamma5(lorentz_parity, particle_parity) * GS(pi->four_momentum(kin));
			}
	));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["R12"]},
                    {P["R12"]},
                    {P["Pion"]}
            },
            [](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                auto const& R = edges[0];
                auto const& N = edges[1];
                auto const& pi = edges[2];
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "Pion");
                auto const g = couplings.get(coupl_str);
                auto const m_pi = pi->particle()->mass();
                auto const iso = (R->back() == N->front())? isospin(R, N, pi) : isospin(N, R, pi);
                int const lorentz_parity = 1;
                auto const& pR = R->four_momentum(kin);
                auto const& pN = N->four_momentum(kin);
                auto form_factor1 = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), pi->particle(), pR.E().real());
                auto form_factor2 = N->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(N->particle(), R->particle(), pi->particle(), pN.E().real());
                int const particle_parity = R->particle()->parity() * N->particle()->parity() * pi->particle()->parity();
                //auto isospin = R->particle()->isospin().j() == 1.5 ? isospin_T(R, N) : isospin_tau(R, N);
                return 1.i * g / m_pi * form_factor1 * form_factor2 * iso * gamma5(lorentz_parity, particle_parity) * GS(pi->four_momentum(kin));
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
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "f0_500");
                auto const g = couplings.get(coupl_str);
				auto const m_pi = P.get("pi0")->mass();
				int const lorentz_parity = 1;
				auto const pR = R->four_momentum(kin);
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), f0->particle(), pR.E().real());
				int const particle_parity = R->particle()->parity() * N->particle()->parity() * f0->particle()->parity();
				return g/m_pi * form_factor * gamma5(lorentz_parity, particle_parity) * GS(f0->four_momentum(kin));
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
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "eta");
                auto const g = couplings.get(coupl_str);
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
                auto const coupl_str = coupling_string(R->particle()->name(), N->particle()->name(), gamma->particle()->name());
                auto const g = couplings.get(coupl_str);
				auto const m_rho = P["rho0"]->mass();
				auto const pR = R->four_momentum(kin);
				int const lorentz_parity = -1;
				int const particle_parity = R->particle()->parity() * N->particle()->parity() * gamma->particle()->parity();
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), gamma->particle(), pR.E().real());
				return -g/m_rho * form_factor * dirac_sigmac(gamma->four_momentum(kin), gamma->lorentz_indices()[0]) * gamma5(lorentz_parity, particle_parity);
			}
	));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["R12"]},
                    {P["N"]},
                    {P.get("Rho")}
            },
            [&](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                auto const& R = edges[0];
                auto const& N = edges[1];
                auto const& Rho = edges[2];
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "Rho");
                auto const g = couplings.get(coupl_str);
                auto const m_rho = Rho->particle()->mass();
                auto const pR = R->four_momentum(kin);
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), Rho->particle(), pR.E().real());
	            int const lorentz_parity = -1;
	            int const particle_parity = R->particle()->parity() * N->particle()->parity() * Rho->particle()->parity();
                return -g/m_rho * form_factor * dirac_sigmac(Rho->four_momentum(kin), Rho->lorentz_indices()[0]) * gamma5(lorentz_parity, particle_parity);
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
				auto const mu = R->lorentz_indices()[0];
				auto const& pR = R->four_momentum(kin);
				auto const& pPi = pi->four_momentum(kin);
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "Pion");
                auto const g = couplings.get(coupl_str);
				auto const m_pi = pi->particle()->mass();
				//auto const isospin = R->particle()->isospin().j() == 1.5 ? isospin_T(R, N) : isospin_tau(R, N);
				auto const iso = (R->back() == N->front())? isospin(R, N, pi) : isospin(N, R, pi);
				auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), pi->particle(), pR.E().real());
				int const lorentz_parity = -1;
				int const particle_parity = R->particle()->parity() * N->particle()->parity() * pi->particle()->parity();
				Complex phase = std::exp(2.i * M_PI/360. * R->particle()->user_data<double>("phase"));
				auto const g5 = gamma5(lorentz_parity, particle_parity);
				return 1.i * phase * form_factor * iso * g/(m_pi*m_pi) * O32c(pR, pPi, mu) * g5;
			}
	));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["Pion"], Edge_Direction::OUT},
                    {P["N"], Edge_Direction::IN},
                    {P["N"], Edge_Direction::OUT},
                    {Feynumeric::QED::Photon, Edge_Direction::IN}
            },
            [&](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                using namespace Feynumeric::Units;

                auto const pion = edges[0];
                auto const Ni= edges[1];
                auto const Nf = edges[2];
                auto const photon = edges[3];
                auto const mu = photon->lorentz_indices()[0];

                auto const mpi = pion->particle()->mass();
                auto const mpi2 = mpi*mpi;
                auto const mN = Ni->particle()->mass();
                auto const mN2 = mN*mN;

                auto const q = pion->four_momentum(kin);
                auto const k = photon->four_momentum(kin);
                auto const pi = Ni->four_momentum(kin);
                auto const pf = Nf->four_momentum(kin);

                auto const s = dot(pi+k, pi+k);
                auto const t = dot(pi-pf, pi-pf);
                auto const u = dot(pi-q, pi-q);


                auto const Lambda = 0.9;
                auto const Lambda4 = Lambda * Lambda * Lambda * Lambda;

                auto const F1 = 1./(1.+(s-mN2)*(s-mN2)/Lambda4);
                auto const F2 = 1./(1.+(u-mN2)*(u-mN2)/Lambda4);;
                auto const F3 = 1./(1.+(t-mpi2)*(t-mpi2)/Lambda4);;
                auto const Fhat = F1 + F2 + F3 - F1*F2 - F1*F3 - F2*F3 + F1 * F2 * F3;

                auto const fNNpi = 0.97;
                auto const kappa_p = 1.;
                auto const kappa_n = 1.;

                auto const A1 = -(
                        1/(2.*mN) * (F1 * kappa_n + F2 * kappa_p) +
                        2.*mN*F2/(u-mN2) * (1.+kappa_p) +
                        2.*mN*F1/(s+mN2) * kappa_n
                        );
                auto const A2 = 4*mN*Fhat/((t-mpi2)*(u-mN2));
                auto const A3 = 2*kappa_n * F1/(s-mN2);
                auto const A4 = 2*kappa_p * F2/(u-mN2);
                auto const M1 = GA5 * (GAC[*mu] * GS(k) - k.co(mu));
                auto const M2 = GA5/2. * (
                        (2.*pi.co(mu) - k.co(mu))*(2.*dot(q, k) - dot(k, k)) -
                        (2.*q.co(mu) - k.co(mu))*(2.*dot(pi, k) - dot(k, k))
                        );
                auto const M3 = GA5/2. * (
                        GAC[*mu] * (2.*dot(pf, k) + dot(k, k)) - (2.*pf.co(mu) + k.co(mu))*GS(k)
                        );
                auto const M4 = GA5/2. * (
                        GAC[*mu] * (2.*dot(pi, k) + dot(k, k)) - (2.*pi.co(mu) + k.co(mu))*GS(k)
                );

                return std::sqrt(2.) * 1._e * fNNpi/mpi * (A1 * M1 + A2 * M2 + A3 * M3 + A4 * M4);
            }
    ));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["Pion"], Edge_Direction::OUT},
                    {P["N"], Edge_Direction::IN},
                    {P["N"], Edge_Direction::OUT},
                    {P["Rho"], Edge_Direction::IN}
            },
            [&](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                using namespace Feynumeric::Units;

                auto const pion = edges[0];
                auto const Ni= edges[1];
                auto const Nf = edges[2];
                auto const photon = edges[3];
                auto const mu = photon->lorentz_indices()[0];

                auto const mpi = pion->particle()->mass();
                auto const mpi2 = mpi*mpi;
                auto const mN = Ni->particle()->mass();
                auto const mN2 = mN*mN;

                auto const q = pion->four_momentum(kin);
                auto const k = photon->four_momentum(kin);
                auto const pi = Ni->four_momentum(kin);
                auto const pf = Nf->four_momentum(kin);

                auto const s = dot(pi+k, pi+k);
                auto const t = dot(pi-pf, pi-pf);
                auto const u = dot(pi-q, pi-q);


                auto const Lambda = 0.9;
                auto const Lambda4 = Lambda * Lambda * Lambda * Lambda;

                auto const F1 = 1./(1.+(s-mN2)*(s-mN2)/Lambda4);
                auto const F2 = 1./(1.+(u-mN2)*(u-mN2)/Lambda4);;
                auto const F3 = 1./(1.+(t-mpi2)*(t-mpi2)/Lambda4);;
                auto const Fhat = F1 + F2 + F3 - F1*F2 - F1*F3 - F2*F3 + F1 * F2 * F3;

                auto const fNNpi = 0.97;
                auto const kappa_p = 1.;
                auto const kappa_n = 1.;
                auto const kappa_rho = 1.;
                auto const gtildeRho = 1.;

                auto const A1 = (
                        kappa_rho/(2.*mN) * (F2-F1) +
                        2.*mN * (1.+kappa_rho) * ( F2/(u-mN2) - F1/(s-mN2))
                );
                auto const A2 = -4*mN*Fhat/((t-mpi2)*(u-mN2));
                auto const A3 = 2*kappa_rho * F1/(s-mN2);
                auto const A4 = -2.*kappa_rho * F2/(u-mN2);
                auto const A5 = - 4.*mN * Fhat / ( (t-mpi2)*(s-mN2));

                auto const M1 = GA5 * (GAC[*mu] * GS(k) - k.co(mu));
                auto const M2 = GA5/2. * (
                        (2.*pi.co(mu) - k.co(mu))*(2.*dot(q, k) - dot(k, k)) -
                        (2.*q.co(mu) - k.co(mu))*(2.*dot(pi, k) - dot(k, k))
                );
                auto const M3 = GA5/2. * (
                        GAC[*mu] * (2.*dot(pf, k) + dot(k, k)) - (2.*pf.co(mu) + k.co(mu))*GS(k)
                );
                auto const M4 = GA5/2. * (
                        GAC[*mu] * (2.*dot(pi, k) + dot(k, k)) - (2.*pi.co(mu) + k.co(mu))*GS(k)
                );
                auto const M5 = GA5/2. * (
                    (2.*pf.co(mu) + k.co(mu))*(2.*dot(q, k) - dot(k, k)) -
                    (2.*q.co(mu) - k.co(mu))*(2.*dot(pf, k) - dot(k, k))
            );

                return  gtildeRho*fNNpi/(std::sqrt(2.) * mpi) * (A1 * M1 + A2 * M2 + A3 * M3 + A4 * M4);
            }
    ));

//
//    VMP->add(Feynumeric::Vertex(
//            {
//                    {P["R32"], Feynumeric::Edge_Direction::OUT},
//                    {P["N"], Feynumeric::Edge_Direction::IN},
//                    {P["Pion"]}
//            },
//            [](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
//                using namespace Feynumeric;
//                static auto const Identity = Matrix(4, 4, 1);
//                auto const& R = edges[0];
//                auto const& N = edges[1];
//                auto const& pi = edges[2];
//                auto mu = R->lorentz_indices()[0];
//                auto const& pR = R->four_momentum(kin);
//                auto const& pPi = pi->four_momentum(kin);
//                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "Pion");
//                auto const g = couplings.get(coupl_str);
//                auto const m_pi = pi->particle()->mass();
//                //auto const isospin = R->particle()->isospin().j() == 1.5 ? isospin_T(R, N) : isospin_tau(R, N);
//                auto const iso = (R->back() == N->front())? isospin(R, N, pi) : isospin(N, R, pi);
//                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), pi->particle(), pR.E().real());
//                int const lorentz_parity = -1;
//                int const particle_parity = R->particle()->parity() * N->particle()->parity() * pi->particle()->parity();
//                Complex phase = std::exp(2.i * M_PI/360. * R->particle()->user_data<double>("phase"));
//                auto const g5 = gamma5(lorentz_parity, particle_parity);
//                return 1.i * phase * form_factor * iso * g/(m_pi*m_pi) * O32c(pR, pPi, mu) * g5;
//            }
//    ));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["R52"]},
                    {P["N"]},
                    {P["Pion"]}
            },
            [](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                static auto const Identity = Matrix(4, 4, 1);
                auto const& R = edges[0];
                auto const& N = edges[1];
                auto const& pi = edges[2];
//                auto mu = R->lorentz_indices()[0];
//                auto nu = R->lorentz_indices()[1];
                auto const& pR = R->four_momentum(kin);
                auto const& pPi = pi->four_momentum(kin);
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "Pion");
                auto const g = couplings.get(coupl_str);
                auto const m_pi = pi->particle()->mass();
                auto const iso = (R->back() == N->front())? isospin(R, N, pi) : isospin(N, R, pi);
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), pi->particle(), pR.E().real());
                int const lorentz_parity = 1;
                int const particle_parity = R->particle()->parity() * N->particle()->parity() * pi->particle()->parity();
                Complex phase = std::exp(2.i * M_PI/360. * R->particle()->user_data<double>("phase"));
                auto const g5 = gamma5(lorentz_parity, particle_parity);
                return -1.i * phase * form_factor * iso * g/(m_pi*m_pi) * Oc(pR, {pPi, pPi}, R->lorentz_indices()) * g5;
            }
    ));


    VMP->add(Feynumeric::Vertex(
            {
                    {P["R72"]},
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
                auto nu = R->lorentz_indices()[1];
                auto la = R->lorentz_indices()[2];
                auto const& pR = R->four_momentum(kin);
                auto const& pPi = pi->four_momentum(kin);
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "Pion");
                auto const g = couplings.get(coupl_str);
                auto const m_pi = pi->particle()->mass();
                auto const iso = (R->back() == N->front())? isospin(R, N, pi) : isospin(N, R, pi);
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), pi->particle(), pR.E().real());
                int const lorentz_parity = -1;
                int const particle_parity = R->particle()->parity() * N->particle()->parity() * pi->particle()->parity();
                Complex phase = std::exp(2.i * M_PI/360. * R->particle()->user_data<double>("phase"));
                auto const g5 = gamma5(lorentz_parity, particle_parity);
                return -1.i * phase * form_factor * iso * g/(m_pi*m_pi) * Oc(pR, {pPi, pPi, pPi}, {mu, nu, la}) * g5;
            }
    ));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["R32"]},
                    {P["N"]},
                    {P["eta"]}
            },
            [&](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                static auto const Identity = Matrix(4, 4, 1);
                auto const& R = edges[0];
                auto const& N = edges[1];
                auto const& eta = edges[2];
                auto mu = R->lorentz_indices()[0];
                auto const& pR = R->four_momentum(kin);
                auto const& pPi = eta->four_momentum(kin);
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "eta");
                auto const g = couplings.get(coupl_str);
                auto const m_pi = P.get("pi0")->mass();
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), eta->particle(), pR.E().real());
	            int const lorentz_parity = -1;
	            int const particle_parity = R->particle()->parity() * N->particle()->parity() * eta->particle()->parity();
                Complex phase = std::exp(2.i * M_PI/360. * R->particle()->user_data<double>("phase"));
                return -1.i * phase * form_factor * g / (m_pi * m_pi) * O32c(pR, pPi, mu) * gamma5(lorentz_parity, particle_parity);
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
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), remove_charge_ending(N->particle()->name()), "Pion");
                auto const g = couplings.get(coupl_str);
                auto const m_pi = pi->particle()->mass();
	            int const lorentz_parity = -1;
	            int const particle_parity = R->particle()->parity() * N->particle()->parity() * pi->particle()->parity();
                auto const iso = (R->back() == N->front())? isospin(R, N, pi) : isospin(N, R, pi);
                auto form_factor1 = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), pi->particle(), pR.E().real());
                auto form_factor2 = N->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(N->particle(), R->particle(), pi->particle(), pN.E().real());
                Complex phase = std::exp(2.i * M_PI/360. * R->particle()->user_data<double>("phase"));
                return -1.i * phase * form_factor1 * form_factor2 * iso * g/(m_pi*m_pi) * O32c(pR, pPi, mu) * gamma5(lorentz_parity, particle_parity);
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

                auto const coupl_str = coupling_string(remove_charge_ending(Rout->particle()->name()), remove_charge_ending(Rin->particle()->name()), "Pion");
                auto const g = couplings.get(coupl_str);
				auto const m_pi = pi->particle()->mass();
				//auto const isospin = R->particle()->isospin().j() == 1.5 ? isospin_T(R, N) : isospin_tau(R, N);
				auto const iso = isospin(Rout, Rin, pi);
                auto form_factor1 = Rout->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(Rout->particle(), Rin->particle(), pi->particle(), pRout.E().real());
                auto form_factor2 = Rin->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(Rin->particle(), Rout->particle(), pi->particle(), pRin.E().real());
				//Complex phase = std::exp(2.i * M_PI/360. * R->particle()->user_data<double>("phase"));
				int const lorentz_parity = 1;
				int const particle_parity = Rout->particle()->parity() * Rin->particle()->parity() * pi->particle()->parity();
				return -1.i * iso * form_factor1 * form_factor2 * g/(m_pi*m_pi*m_pi*m_pi ) * O32c(pRout, pPi, mu) * gamma5(lorentz_parity, particle_parity) * O32c(pRin, pPi, nu);
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
                auto const coupl_str = coupling_string(R->particle()->name(), N->particle()->name(), photon->particle()->name());
                auto const g = couplings.get(coupl_str);
				auto const m_rho = P["rho0"]->mass();
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), photon->particle(), pR.E().real());
				int const lorentz_parity = 1;
				int const particle_parity = R->particle()->parity() * N->particle()->parity() * photon->particle()->parity();
				auto result = 1.i * g/(4*m_rho*m_rho) * form_factor * CONTRACT_MATRIX(O32c(pPhoton, mu, kappa) * O32c(pR, mu, lambda) * MT[*mu][*mu], mu)  * gamma5(lorentz_parity, particle_parity);
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
                auto const coupl_str = coupling_string(R->particle()->name(), N->particle()->name(), photon->particle()->name());
                auto const g = couplings.get(coupl_str);
				auto const m_rho = P["rho0"]->mass();
				int const lorentz_parity = 1;
				int const particle_parity = R->particle()->parity() * N->particle()->parity() * photon->particle()->parity();
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), photon->particle(), pR.E().real());
				auto result = 1.i * g/(4*m_rho*m_rho) * form_factor * CONTRACT_MATRIX(O32c(pR, mu,lambda) * O32c(pPhoton, mu, kappa) * MT[*mu][*mu], mu) * gamma5(lorentz_parity, particle_parity);
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
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "Rho");
                auto const g = couplings.get(coupl_str);
				auto const m_rho = rho->particle()->mass();
				int const lorentz_parity = 1;
				int const particle_parity = R->particle()->parity() * N->particle()->parity() * rho->particle()->parity();
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), rho->particle(), pR.E().real());
				auto result = 1.i * g/(4*m_rho*m_rho) * form_factor * CONTRACT_MATRIX(O32c(pRho, mu, kappa) * O32c(pR, mu, lambda) * MT[*mu][*mu], mu) * gamma5(lorentz_parity, particle_parity);
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
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "Rho");
                auto const g = couplings.get(coupl_str);
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), Rho->particle(), pR.E().real());
				auto const m_rho = Rho->particle()->mass();
				int const lorentz_parity = 1;
				int const particle_parity = R->particle()->parity() * N->particle()->parity() * Rho->particle()->parity();
				auto result = 1.i * g/(4*m_rho*m_rho) * form_factor * CONTRACT_MATRIX(O32c(pR, mu, lambda) * O32c(pRho, mu, kappa) * MT[*mu][*mu], mu) * gamma5(lorentz_parity, particle_parity);
				return result;
			}
	));


    VMP->add(Feynumeric::Vertex(
            {
                    {P["R52"], Edge_Direction::IN},
                    {P["N"], Edge_Direction::OUT},
                    {Feynumeric::QED::Photon}
            },
            [&](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                auto const& R = edges[0];
                auto const& N = edges[1];
                auto const& photon = edges[2];
                auto kappa = photon->lorentz_indices()[0];
                auto mu = R->lorentz_indices()[0];
                auto nu = R->lorentz_indices()[1];
                auto const& pN = N->four_momentum(kin);
                auto const& pR = R->four_momentum(kin);
                auto const& pPhoton = photon->four_momentum(kin);
                auto const coupl_str = coupling_string(R->particle()->name(), N->particle()->name(), photon->particle()->name());
                auto const g = couplings.get(coupl_str);
                auto const m_rho = P["rho0"]->mass();
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), photon->particle(), pR.E().real());
                int const lorentz_parity = 1;
                int const particle_parity = R->particle()->parity() * N->particle()->parity() * photon->particle()->parity();
//                auto result = 1.i * g/(4*m_rho*m_rho) * form_factor * CONTRACT_MATRIX(O32c(pPhoton, mu, kappa) * O32c(pR, mu, lambda) * MT[*mu][*mu], mu)  * gamma5(lorentz_parity, particle_parity);

                Matrix result(4, 4, 1.);

                auto alpha = std::make_shared<Lorentz_Index>();
                auto beta = std::make_shared<Lorentz_Index>();

                for( int i = 0; i < 4; ++i ){
                    for( int j = 0; j < 4; ++j ){
                        result += O32c(pPhoton, alpha, kappa) * Oc(pR, {alpha, beta}, {mu, nu}) * pN.contra(beta) * MT[*alpha][*alpha];
                        ++(*beta);
                    }
                    ++(*alpha);
                }
                return 1.i * g/std::pow(2*m_rho, 3) * form_factor * result;
            }
    ));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["R52"], Edge_Direction::OUT},
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
                auto mu = R->lorentz_indices()[0];
                auto nu = R->lorentz_indices()[0];
                auto const& pN = R->four_momentum(kin);
                auto const& pR = R->four_momentum(kin);
                auto const& pPhoton = photon->four_momentum(kin);
                auto const coupl_str = coupling_string(R->particle()->name(), N->particle()->name(), photon->particle()->name());
                auto const g = couplings.get(coupl_str);
                auto const m_rho = P["rho0"]->mass();
                int const lorentz_parity = 1;
                int const particle_parity = R->particle()->parity() * N->particle()->parity() * photon->particle()->parity();
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), photon->particle(), pR.E().real());
                Matrix result(4, 4, 1.);

                auto alpha = std::make_shared<Lorentz_Index>();
                auto beta = std::make_shared<Lorentz_Index>();

                for( int i = 0; i < 4; ++i ){
                    for( int j = 0; j < 4; ++j ){
                        result += Oc(pR, {alpha, beta}, {mu, nu}) * O32c(pPhoton, alpha, kappa) * pN.contra(beta) * MT[*alpha][*alpha];
                        ++(*beta);
                    }
                    ++(*alpha);
                }
                return 1.i * g/std::pow(2*m_rho, 3) * form_factor * result;
            }
    ));


    VMP->add(Feynumeric::Vertex(
            {
                    {P.get("R52"), Edge_Direction::IN},
                    {P.get("N"), Edge_Direction::OUT},
                    {P.get("Rho")}
            },
            [&](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                auto const& R = edges[0];
                auto const& N = edges[1];
                auto const& rho = edges[2];
                auto kappa = rho->lorentz_indices()[0];
                auto mu = R->lorentz_indices()[0];
                auto nu = R->lorentz_indices()[0];
                auto const& pN = N->four_momentum(kin);
                auto const& pR = R->four_momentum(kin);
                auto const& pRho = rho->four_momentum(kin);
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "Rho");
                auto const g = couplings.get(coupl_str);
                auto const m_rho = rho->particle()->mass();
                int const lorentz_parity = 1;
                int const particle_parity = R->particle()->parity() * N->particle()->parity() * rho->particle()->parity();
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), rho->particle(), pR.E().real());


                Matrix result(4, 4, 1.);

                auto alpha = std::make_shared<Lorentz_Index>();
                auto beta = std::make_shared<Lorentz_Index>();

                for( int i = 0; i < 4; ++i ){
                    for( int j = 0; j < 4; ++j ){
                        result += O32c(pRho, alpha, kappa) * Oc(pR, {alpha, beta}, {mu, nu}) * pN.contra(beta) * MT[*alpha][*alpha];
                        ++(*beta);
                    }
                    ++(*alpha);
                }
                return 1.i * g/std::pow(2*m_rho, 3) * form_factor * result;
            }
    ));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["R52"], Edge_Direction::OUT},
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
                auto mu = R->lorentz_indices()[0];
                auto nu = R->lorentz_indices()[1];
                auto const& pN = N->four_momentum(kin);
                auto const& pR = R->four_momentum(kin);
                auto const& pRho = Rho->four_momentum(kin);
                auto const coupl_str = coupling_string(remove_charge_ending(R->particle()->name()), "N", "Rho");
                auto const g = couplings.get(coupl_str);
                auto form_factor = R->particle()->user_data<FORM_FACTOR_FUNCTION>("form_factor")(R->particle(), N->particle(), Rho->particle(), pR.E().real());
                auto const m_rho = Rho->particle()->mass();
                int const lorentz_parity = 1;
                int const particle_parity = R->particle()->parity() * N->particle()->parity() * Rho->particle()->parity();

                Matrix result(4, 4, 1.);

                auto alpha = std::make_shared<Lorentz_Index>();
                auto beta = std::make_shared<Lorentz_Index>();

                for( int i = 0; i < 4; ++i ){
                    for( int j = 0; j < 4; ++j ){
                        result += Oc(pR, {alpha, beta}, {mu, nu}) * O32c(pRho, alpha, kappa) * pN.contra(beta) * MT[*alpha][*alpha];
                        ++(*beta);
                    }
                    ++(*alpha);
                }
                return 1.i * g/std::pow(2*m_rho, 3) * form_factor * result;
            }
    ));

    VMP->add(Feynumeric::Vertex(
            {
                    {P["N"], Edge_Direction::IN},
                    {P["N"], Edge_Direction::OUT},
                    {P["Pion"], Edge_Direction::IN},
                    {P["Pion"], Edge_Direction::OUT}
            },
            [](Feynumeric::Kinematics const& kin, std::vector<std::shared_ptr<Feynumeric::Graph_Edge>> const& edges){
                using namespace Feynumeric;
                auto const& N1 = edges[0];
                auto const& N2 = edges[1];
                auto const& pi1 = edges[2];
                auto const& pi2 = edges[3];
                auto const m_pi = pi1->particle()->mass();
                auto const p_N1 = N1->four_momentum(kin);
                auto const p_N2 = N2->four_momentum(kin);
                auto const p_pi1 = pi1->four_momentum(kin);
                auto const p_pi2 = pi2->four_momentum(kin);
                return Matrix(1, 1, dot(p_pi1, p_pi2)/(m_pi*m_pi));
            }
    ));

}