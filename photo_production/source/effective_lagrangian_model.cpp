#include <feynumeric/contract.hpp>
#include <feynumeric/dirac.hpp>
#include <feynumeric/particle.hpp>
#include <feynumeric/particle_manager.hpp>
#include <feynumeric/qed.hpp>
#include <feynumeric/units.hpp>
#include <feynumeric/messages.hpp>
#include <feynumeric/constexpr_math.hpp>


#include "effective_lagrangian_model.hpp"

using namespace Feynumeric::Units;
using Feynumeric::Particle;
using Feynumeric::Matrix;

/**************************************************************************\
|* VALUES TAKEN FROM https://pdg.lbl.gov/2021/tables/contents_tables.html *|
\**************************************************************************/

//Feynumeric::Particle_Ptr Photon   = std::make_shared<Feynumeric::Particle>("Photon", Particle::Type::NeutralBoson, 0, 0, 0, 1);

//Feynumeric::Particle_Ptr Proton   = std::make_shared<Feynumeric::Particle>("Proton", Particle::Type::TrueBaryon, 938.270813_MeV, 0, +1, 0.5);
//Feynumeric::Particle_Ptr Neutron  = std::make_shared<Feynumeric::Particle>("Neutron", Particle::Type::TrueBaryon, 939.5654133_MeV, 0, 0, 0.5);

//Feynumeric::Particle_Ptr Pi_Zero  = std::make_shared<Feynumeric::Particle>("Pi_0", Particle::Type::NeutralMeson, 134.9768_MeV, 0, 0, 0);
//Feynumeric::Particle_Ptr Pi_Plus  = std::make_shared<Feynumeric::Particle>("Pi_+", Particle::Type::AntiMeson, 139.57039_MeV, 0, +1, 0);
//Feynumeric::Particle_Ptr Pi_Minus = std::make_shared<Feynumeric::Particle>("Pi_-", Particle::Type::TrueMeson, 139.57039_MeV, 0, -1, 0);

//Feynumeric::Particle_Ptr Rho_Zero  = std::make_shared<Feynumeric::Particle>("Rho_0", Particle::Type::NeutralMeson, 775.26_MeV, 147.8_MeV, 0, 1);
//Feynumeric::Particle_Ptr Rho_Plus  = std::make_shared<Feynumeric::Particle>("Rho_+", Particle::Type::AntiMeson, 775.11_MeV, 149.1_MeV, +1, 1);
//Feynumeric::Particle_Ptr Rho_Minus = std::make_shared<Feynumeric::Particle>("Rho_-", Particle::Type::TrueMeson, 775.11_MeV, 149.1_MeV, -1, 1);

//Feynumeric::Particle_Ptr N1440p   = std::make_shared<Feynumeric::Particle>("N1440_+", Particle::Type::TrueBaryon, 1.44_GeV, 350._MeV, +1, 0.5);
//Feynumeric::Particle_Ptr N1440n   = std::make_shared<Feynumeric::Particle>("N1440_0", Particle::Type::TrueBaryon, 1.44_GeV, 350._MeV, 0, 0.5);

//Feynumeric::Particle_Ptr N1520p   = std::make_shared<Feynumeric::Particle>("N1520_+", Particle::Type::TrueBaryon, 1.515_GeV, 110._MeV, +1, 1.5);
//Feynumeric::Particle_Ptr N1520n   = std::make_shared<Feynumeric::Particle>("N1520_0", Particle::Type::TrueBaryon, 1.515_GeV, 110._MeV, 0, 1.5);

//Feynumeric::Particle_Ptr N1535p   = std::make_shared<Feynumeric::Particle>("N1535_+", Particle::Type::TrueBaryon, 1.53_GeV, 150._MeV, +1, 0.5);
//Feynumeric::Particle_Ptr N1535n   = std::make_shared<Feynumeric::Particle>("N1535_0", Particle::Type::TrueBaryon, 1.53_GeV, 150._MeV, 0, 0.5);

//Feynumeric::Particle_Ptr N1650p   = std::make_shared<Feynumeric::Particle>("N1535_+", Particle::Type::TrueBaryon, 1.65_GeV, 125._MeV, +1, 0.5);
//Feynumeric::Particle_Ptr N1650n   = std::make_shared<Feynumeric::Particle>("N1535_0", Particle::Type::TrueBaryon, 1.65_GeV, 125._MeV, 0, 0.5);

//double dyson_factor(double p2, double m, int l)
//{
//	double qR = Feynumeric::momentum(std::sqrt(p2), Proton->mass(), Pi_Minus->mass());
//	double q0 = Feynumeric::momentum(m, Proton->mass(), Pi_Minus->mass());
//	return std::pow(qR/q0, 2*l+1);
//}

void init_particles()
{

//    using namespace Feynumeric;
//    QED::init_particles();
//    Proton->isospin(Angular_Momentum(0.5, 0.5));
//	Neutron->isospin(Angular_Momentum(0.5, -0.5));
//    Pi_Zero->isospin(Angular_Momentum(1, 0));
//	Pi_Plus->isospin(Angular_Momentum(1, 1));
//	Pi_Minus->isospin(Angular_Momentum(1, -1));

    /********************************************\
	|*	███    ██  ██ ██   ██ ██   ██  ██████   *|
	|*	████   ██ ███ ██   ██ ██   ██ ██  ████  *|
	|*	██ ██  ██  ██ ███████ ███████ ██ ██ ██  *|
	|*	██  ██ ██  ██      ██      ██ ████  ██  *|
	|*	██   ████  ██      ██      ██  ██████   *|
    \********************************************/

//	N1440p->parity(1);
//    N1440p->user_data("gRNpi", 0.380006);
//	N1440p->user_data("gRNrho", 0.0528755);
//	N1440p->user_data("gRNgamma", 0.0620974);
//    N1440p->isospin(Angular_Momentum(0.5, 0.5));
//    N1440p->user_data("branching_N_pi_upper", 75._percent);
//	N1440p->user_data("branching_N_pi_lower", 55._percent);
//	N1440p->user_data("branching_N_rho_upper", 0._percent);
//	N1440p->user_data("branching_N_rho_lower", 0._percent);
//	N1440p->user_data("branching_proton_photon_upper", 0.048_percent);
//	N1440p->user_data("branching_proton_photon_lower", 0.035_percent);
//	N1440p->user_data("branching_neutron_photon_upper", 0.04_percent);
//	N1440p->user_data("branching_neutron_photon_lower", 0.02_percent);
//    N1440p->width([&](double p2){
//    	if( p2 < 0 )
//	    {
//    		return 0.;
//	    }
//    	int l = 1;
//    	return N1440p->width() * dyson_factor(p2, N1440p->mass(), l);
//    });
//
//	N1440n->parity(1);
//	N1440n->user_data("gRNpi", 0.379618);
//	N1440n->user_data("gRNrho", 1.);
//	N1440n->user_data("gRneutron_photon", 0.0529588);
//	N1440n->isospin(Angular_Momentum(0.5, -0.5));
//	N1440n->user_data("branching_N_pi_upper", N1440p->user_data<long double>("branching_N_pi_upper"));
//	N1440n->user_data("branching_N_pi_lower", N1440p->user_data<long double>("branching_N_pi_lower"));
//	N1440n->user_data("branching_N_rho_upper", N1440p->user_data<long double>("branching_N_rho_upper"));
//	N1440n->user_data("branching_N_rho_lower", N1440p->user_data<long double>("branching_N_rho_lower"));
//	N1440n->user_data("branching_proton_photon_upper", N1440p->user_data<long double>("branching_proton_photon_upper"));
//	N1440n->user_data("branching_proton_photon_lower", N1440p->user_data<long double>("branching_proton_photon_lower"));
//	N1440n->user_data("branching_neutron_photon_upper", N1440p->user_data<long double>("branching_neutron_photon_upper"));
//	N1440n->user_data("branching_neutron_photon_lower", N1440p->user_data<long double>("branching_neutron_photon_lower"));
//
//	N1440n->width([&](double p2){ return N1440n->width();});

	/********************************************\
	|*	███    ██  ██ ███████ ██████   ██████   *|
	|*	████   ██ ███ ██           ██ ██  ████  *|
	|*	██ ██  ██  ██ ███████  █████  ██ ██ ██  *|
	|*	██  ██ ██  ██      ██ ██      ████  ██  *|
	|*	██   ████  ██ ███████ ███████  ██████   *|
	\********************************************/

//	N1520p->parity(-1);
//	N1520p->user_data("gRNpi", 1.);
//	N1520p->user_data("gRNrho", 1.);
//	N1520p->isospin(Angular_Momentum(0.5, 0.5));
//	N1520p->user_data("branching_N_pi_upper", 65._percent);
//	N1520p->user_data("branching_N_pi_lower", 55._percent);
//	N1520p->user_data("branching_N_rho_upper", 13.7_percent);
//	N1520p->user_data("branching_N_rho_lower", 9.9_percent);
//	N1520p->width([&](double p2){ return N1520p->width();});
//
//	N1520n->parity(-1);
//	N1520n->user_data("gRNpi", 1.);
//	N1520n->user_data("gRNrho", 1.);
//	N1520n->isospin(Angular_Momentum(0.5, -0.5));
//	N1520n->user_data("branching_N_pi_upper", N1520p->user_data<long double>("branching_N_pi_upper"));
//	N1520n->user_data("branching_N_pi_lower", N1520p->user_data<long double>("branching_N_pi_lower"));
//	N1520n->user_data("branching_N_rho_upper", N1520p->user_data<long double>("branching_N_rho_upper"));
//	N1520n->user_data("branching_N_rho_lower", N1520p->user_data<long double>("branching_N_rho_lower"));
//	N1520n->width([&](double p2){ return N1520n->width();});
//
	/********************************************\
	|*	███    ██  ██ ███████ ██████  ███████
	|*	████   ██ ███ ██           ██ ██
		██ ██  ██  ██ ███████  █████  ███████
		██  ██ ██  ██      ██      ██      ██
		██   ████  ██ ███████ ██████  ███████
	\********************************************/

//	N1535p->parity(-1);
//	N1535p->user_data("gRNpi", 1.);
//	N1535p->user_data("gRNrho", 1.);
//	N1535p->isospin(Angular_Momentum(0.5, 0.5));
//	N1535p->user_data("branching_N_pi_upper", 65._percent);
//	N1535p->user_data("branching_N_pi_lower", 55._percent);
//	N1535p->user_data("branching_N_rho_upper", 13.7_percent);
//	N1535p->user_data("branching_N_rho_lower", 9.9_percent);
//	N1535p->width([&](double p2){ return N1535p->width();});

//	N1535n->parity(-1);
//	N1535n->user_data("gRNpi", 1.);
//	N1535n->user_data("gRNrho", 1.);
//	N1535n->isospin(Angular_Momentum(0.5, -0.5));
//	N1535n->user_data("branching_N_pi_upper", N1535p->user_data<long double>("branching_N_pi_upper"));
//	N1535n->user_data("branching_N_pi_lower", N1535p->user_data<long double>("branching_N_pi_lower"));
//	N1535n->user_data("branching_N_rho_upper", N1535p->user_data<long double>("branching_N_rho_upper"));
//	N1535n->user_data("branching_N_rho_lower", N1535p->user_data<long double>("branching_N_rho_lower"));
//	N1535n->width([&](double p2){ return N1535n->width();});


}

Feynumeric::Vertex_Manager_Ptr VMP = std::make_shared<Feynumeric::Vertex_Manager>();

double isospin2_2(Feynumeric::Feynman_Graph::Edge_Ptr a, Feynumeric::Feynman_Graph::Edge_Ptr b)
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
	using Feynumeric::Direction;
	using namespace Feynumeric::Units;

	Feynumeric::QED::init_vertices();
	VMP->import(Feynumeric::QED::VMP);

	VMP->add(Feynumeric::Vertex(
			{
					{P["N12p"], Direction::ANY},
					{P["N"], Direction::ANY},
					{P["Pion"],   Direction::ANY}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){

				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& pi = edges[2];
				auto const g = std::any_cast<double>(R->particle()->user_data("gRNpi"));
				auto const m_pi = pi->particle()->mass();
				//auto constexpr isospin = constexpr_sqrt(2.);
				auto isospin = isospin2_2(R, N);
				return -g/m_pi * isospin * GA5 * GS(pi->four_momentum(kin));
			}
	));
	/*
	VMP->add(Feynumeric::Vertex(
			{
					{N1440p, Direction::ANY},
					{Proton, Direction::ANY},
					{Feynumeric::QED::Photon, Direction::ANY}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
				using namespace Feynumeric;
				auto const& n1440 = edges[0];
				auto const& photon = edges[2];
				auto const g = std::any_cast<double>(n1440->particle()->user_data("gRproton_photon"));
				auto const m_rho = Rho_Zero->mass();
				return -g/m_rho * dirac_sigmac(photon->four_momentum(kin), photon->lorentz_indices()[0]);
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{N1440n, Direction::ANY},
					{Neutron, Direction::ANY},
					{Feynumeric::QED::Photon, Direction::ANY}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
				using namespace Feynumeric;
				auto const& n1440 = edges[0];
				auto const& photon = edges[2];
				auto const g = std::any_cast<double>(n1440->particle()->user_data("gRneutron_photon"));
				auto const m_rho = Rho_Zero->mass();
				return -g/m_rho * dirac_sigmac(photon->four_momentum(kin), photon->lorentz_indices()[0]);
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{N1520p, Direction::ANY},
					{Proton, Direction::ANY},
					{Pi_Zero, Direction::ANY}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& pi = edges[2];
				auto mu = R->lorentz_indices()[0];
				auto const& pR = R->four_momentum(kin);
				auto const& pPi = pi->four_momentum(kin);
				auto const g = std::any_cast<double>(R->particle()->user_data("gRNpi"));
				auto const m_pi = pi->particle()->mass();
				auto const isospin = isospin2_2(R, N);
				return 1.i * isospin * g/(m_pi*m_pi) *
					CONTRACT_MATRIX(O32c(pR, mu, lambda) * pPi.contra(lambda), lambda) * GA5;

			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{N1520p, Direction::ANY},
					{Neutron, Direction::ANY},
					{Pi_Plus, Direction::ANY}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& pi = edges[2];
				auto mu = R->lorentz_indices()[0];
				auto const& pR = R->four_momentum(kin);
				auto const& pPi = pi->four_momentum(kin);
				auto const g = std::any_cast<double>(R->particle()->user_data("gRNpi"));
				auto const m_pi = pi->particle()->mass();
				auto const isospin = isospin2_2(R, N);
				return 1.i * isospin * g/(m_pi*m_pi) *
				       CONTRACT_MATRIX(O32c(pR, mu, lambda) * pPi.contra(lambda), lambda) * GA5;

			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{N1520n, Direction::ANY},
					{Proton, Direction::ANY},
					{Pi_Minus, Direction::ANY}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& pi = edges[2];
				auto mu = R->lorentz_indices()[0];
				auto const& pR = R->four_momentum(kin);
				auto const& pPi = pi->four_momentum(kin);
				auto const g = std::any_cast<double>(R->particle()->user_data("gRNpi"));
				auto const m_pi = pi->particle()->mass();
				auto const isospin = isospin2_2(R, N);
				return 1.i * isospin * g/(m_pi*m_pi) *
				       CONTRACT_MATRIX(O32c(pR, mu, lambda) * pPi.contra(lambda), lambda) * GA5;

			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{N1520n, Direction::ANY},
					{Neutron, Direction::ANY},
					{Pi_Zero, Direction::ANY}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
				using namespace Feynumeric;
				auto const& R = edges[0];
				auto const& N = edges[1];
				auto const& pi = edges[2];
				auto mu = R->lorentz_indices()[0];
				auto const& pR = R->four_momentum(kin);
				auto const& pPi = pi->four_momentum(kin);
				auto const g = std::any_cast<double>(R->particle()->user_data("gRNpi"));
				auto const m_pi = pi->particle()->mass();
				auto const isospin = isospin2_2(R, N);
				return 1.i * isospin * g/(m_pi*m_pi) *
				       CONTRACT_MATRIX(O32c(pR, mu, lambda) * pPi.contra(lambda), lambda) * GA5;

			}
	));
	 */
}