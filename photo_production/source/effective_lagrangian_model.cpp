#include <feynumeric/dirac.hpp>
#include <feynumeric/particle.hpp>
#include <feynumeric/qed.hpp>
#include <feynumeric/units.hpp>


#include "effective_lagrangian_model.hpp"

using namespace Feynumeric::Units;
using Feynumeric::Particle;
using Feynumeric::Matrix;

/**************************************************************************\
|* VALUES TAKEN FROM https://pdg.lbl.gov/2021/tables/contents_tables.html *|
\**************************************************************************/

Feynumeric::Particle_Ptr Photon   = std::make_shared<Feynumeric::Particle>("Photon", Particle::Type::NeutralBoson, 0, 0, 0, 1);

Feynumeric::Particle_Ptr Proton   = std::make_shared<Feynumeric::Particle>("Proton", Particle::Type::TrueBaryon, 938.270813_MeV, 0, +1, 0.5);
Feynumeric::Particle_Ptr Neutron  = std::make_shared<Feynumeric::Particle>("Neutron", Particle::Type::TrueBaryon, 939.5654133_MeV, 0, 0, 0.5);

Feynumeric::Particle_Ptr Pi_Zero  = std::make_shared<Feynumeric::Particle>("Pi_0", Particle::Type::NeutralMeson, 134.9768_MeV, 0, 0, 0);
Feynumeric::Particle_Ptr Pi_Plus  = std::make_shared<Feynumeric::Particle>("Pi_+", Particle::Type::AntiMeson, 139.57039_MeV, 0, +1, 0);
Feynumeric::Particle_Ptr Pi_Minus = std::make_shared<Feynumeric::Particle>("Pi_-", Particle::Type::TrueMeson, 139.57039_MeV, 0, -1, 0);

Feynumeric::Particle_Ptr Rho_Zero  = std::make_shared<Feynumeric::Particle>("Rho_0", Particle::Type::NeutralMeson, 775.26_MeV, 147.8_MeV, 0, 1);
Feynumeric::Particle_Ptr Rho_Plus  = std::make_shared<Feynumeric::Particle>("Rho_+", Particle::Type::AntiMeson, 775.11_MeV, 149.1_MeV, +1, 1);
Feynumeric::Particle_Ptr Rho_Minus = std::make_shared<Feynumeric::Particle>("Rho_-", Particle::Type::TrueMeson, 775.11_MeV, 149.1_MeV, -1, 1);

Feynumeric::Particle_Ptr N1440p   = std::make_shared<Feynumeric::Particle>("N1440_+", Particle::Type::TrueBaryon, 1.44_GeV, 350._MeV, +1, 0.5);
Feynumeric::Particle_Ptr N1440n   = std::make_shared<Feynumeric::Particle>("N1440_0", Particle::Type::TrueBaryon, 1.44_GeV, 350._MeV, 0, 0.5);

Feynumeric::Particle_Ptr N1520p   = std::make_shared<Feynumeric::Particle>("N1520_+", Particle::Type::TrueBaryon, 1.515_GeV, 110._MeV, +1, 1.5);
Feynumeric::Particle_Ptr N1520n   = std::make_shared<Feynumeric::Particle>("N1520_0", Particle::Type::TrueBaryon, 1.515_GeV, 110._MeV, 0, 1.5);

Feynumeric::Particle_Ptr N1535p   = std::make_shared<Feynumeric::Particle>("N1535_+", Particle::Type::TrueBaryon, 1.53_GeV, 150._MeV, +1, 0.5);
Feynumeric::Particle_Ptr N1535n   = std::make_shared<Feynumeric::Particle>("N1535_0", Particle::Type::TrueBaryon, 1.53_GeV, 150._MeV, 0, 0.5);

Feynumeric::Particle_Ptr N1650p   = std::make_shared<Feynumeric::Particle>("N1535_+", Particle::Type::TrueBaryon, 1.65_GeV, 125._MeV, +1, 0.5);
Feynumeric::Particle_Ptr N1650n   = std::make_shared<Feynumeric::Particle>("N1535_0", Particle::Type::TrueBaryon, 1.65_GeV, 125._MeV, 0, 0.5);


void init_particles()
{
    using namespace Feynumeric;
    //QED::init_vertices();
    QED::init_particles();
    Proton->isospin(Angular_Momentum(0.5, 0.5));
    Proton->feynman_virtual  = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
    Proton->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return u(e, kin); };
    Proton->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return ubar(e, kin); };


    Pi_Zero->isospin(Angular_Momentum(1, 0));
    Pi_Zero->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){return Matrix(1,1,1);};
    Pi_Zero->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){return Matrix(1,1,1);};
    Pi_Zero->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){return Matrix(1,1,1);};

    /********************************************\
	|*	███    ██  ██ ██   ██ ██   ██  ██████   *|
	|*	████   ██ ███ ██   ██ ██   ██ ██  ████  *|
	|*	██ ██  ██  ██ ███████ ███████ ██ ██ ██  *|
	|*	██  ██ ██  ██      ██      ██ ████  ██  *|
	|*	██   ████  ██      ██      ██  ██████   *|
    \********************************************/

	N1440p->parity(1);
    N1440p->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
    N1440p->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return u(e, kin); };
    N1440p->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return ubar(e, kin); };

    N1440p->user_data("gRNpi", 1.);
	N1440p->user_data("gRNrho", 1.);
    N1440p->isospin(Angular_Momentum(0.5, 0.5));
    N1440p->user_data("branching_N_pi_upper", 75._percent);
	N1440p->user_data("branching_N_pi_lower", 55._percent);
	N1440p->user_data("branching_N_rho_upper", 0._percent);
	N1440p->user_data("branching_N_rho_lower", 0._percent);
    N1440p->width([&](double p2){ return N1440p->width();});

	N1440n->parity(1);
	N1440n->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
	N1440n->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return u(e, kin); };
	N1440n->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return ubar(e, kin); };

	N1440n->user_data("gRNpi", 1.);
	N1440n->user_data("gRNrho", 1.);
	N1440n->isospin(Angular_Momentum(0.5, -0.5));
	N1440n->user_data("branching_N_pi_upper", N1440p->user_data<long double>("branching_N_pi_upper"));
	N1440n->user_data("branching_N_pi_lower", N1440p->user_data<long double>("branching_N_pi_lower"));
	N1440n->user_data("branching_N_rho_upper", N1440p->user_data<long double>("branching_N_rho_upper"));
	N1440n->user_data("branching_N_rho_lower", N1440p->user_data<long double>("branching_N_rho_lower"));
	N1440n->width([&](double p2){ return N1440n->width();});

	/********************************************\
	|*	███    ██  ██ ███████ ██████   ██████   *|
	|*	████   ██ ███ ██           ██ ██  ████  *|
	|*	██ ██  ██  ██ ███████  █████  ██ ██ ██  *|
	|*	██  ██ ██  ██      ██ ██      ████  ██  *|
	|*	██   ████  ██ ███████ ███████  ██████   *|
	\********************************************/

	N1520p->parity(-1);
	N1520p->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
	N1520p->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return u(e, kin); };
	N1520p->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return ubar(e, kin); };

	N1520p->user_data("gRNpi", 1.);
	N1520p->user_data("gRNrho", 1.);
	N1520p->isospin(Angular_Momentum(0.5, 0.5));
	N1520p->user_data("branching_N_pi_upper", 65._percent);
	N1520p->user_data("branching_N_pi_lower", 55._percent);
	N1520p->user_data("branching_N_rho_upper", 13.7_percent);
	N1520p->user_data("branching_N_rho_lower", 9.9_percent);
	N1520p->width([&](double p2){ return N1520p->width();});

	N1520n->parity(-1);
	N1520n->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
	N1520n->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return u(e, kin); };
	N1520n->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return ubar(e, kin); };

	N1520n->user_data("gRNpi", 1.);
	N1520n->user_data("gRNrho", 1.);
	N1520n->isospin(Angular_Momentum(0.5, -0.5));
	N1520n->user_data("branching_N_pi_upper", N1520p->user_data<long double>("branching_N_pi_upper"));
	N1520n->user_data("branching_N_pi_lower", N1520p->user_data<long double>("branching_N_pi_lower"));
	N1520n->user_data("branching_N_rho_upper", N1520p->user_data<long double>("branching_N_rho_upper"));
	N1520n->user_data("branching_N_rho_lower", N1520p->user_data<long double>("branching_N_rho_lower"));
	N1520n->width([&](double p2){ return N1520n->width();});

	/********************************************\
	|*	███    ██  ██ ███████ ██████  ███████
	|*	████   ██ ███ ██           ██ ██
		██ ██  ██  ██ ███████  █████  ███████
		██  ██ ██  ██      ██      ██      ██
		██   ████  ██ ███████ ██████  ███████
	\********************************************/

	N1535p->parity(-1);
	N1535p->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
	N1535p->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return u(e, kin); };
	N1535p->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return ubar(e, kin); };

	N1535p->user_data("gRNpi", 1.);
	N1535p->user_data("gRNrho", 1.);
	N1535p->isospin(Angular_Momentum(0.5, 0.5));
	N1535p->user_data("branching_N_pi_upper", 65._percent);
	N1535p->user_data("branching_N_pi_lower", 55._percent);
	N1535p->user_data("branching_N_rho_upper", 13.7_percent);
	N1535p->user_data("branching_N_rho_lower", 9.9_percent);
	N1535p->width([&](double p2){ return N1535p->width();});

	N1535n->parity(-1);
	N1535n->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
	N1535n->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return u(e, kin); };
	N1535n->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return ubar(e, kin); };

	N1535n->user_data("gRNpi", 1.);
	N1535n->user_data("gRNrho", 1.);
	N1535n->isospin(Angular_Momentum(0.5, -0.5));
	N1535n->user_data("branching_N_pi_upper", N1535p->user_data<long double>("branching_N_pi_upper"));
	N1535n->user_data("branching_N_pi_lower", N1535p->user_data<long double>("branching_N_pi_lower"));
	N1535n->user_data("branching_N_rho_upper", N1535p->user_data<long double>("branching_N_rho_upper"));
	N1535n->user_data("branching_N_rho_lower", N1535p->user_data<long double>("branching_N_rho_lower"));
	N1535n->width([&](double p2){ return N1535n->width();});


}

Feynumeric::Vertex_Manager_Ptr VMP = std::make_shared<Feynumeric::Vertex_Manager>();

void init_vertices()
{
	using Feynumeric::Direction;
	using namespace Feynumeric::Units;

	Feynumeric::QED::init_vertices();
	VMP->import(Feynumeric::QED::VMP);

	VMP->add(Feynumeric::Vertex(
			{
					{N1440p, Direction::BOTH},
					{Proton, Direction::BOTH},
					{Pi_Zero,   Direction::BOTH}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
				using namespace Feynumeric;
				auto const& n1440 = edges[0];
				auto const& pi = edges[2];
				auto const g = std::any_cast<double>(n1440->particle()->user_data("gRNpi"));
				auto const m_pi = pi->particle()->mass();
				return -g/m_pi * GA5 * GS(pi->four_momentum(kin));
			}
	));

	VMP->add(Feynumeric::Vertex(
			{
					{N1440p, Direction::BOTH},
					{Proton, Direction::BOTH},
					{Feynumeric::QED::Photon, Direction::BOTH}
			},
			[](Feynumeric::Kinematics const& kin, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
				using namespace Feynumeric;
				auto const& n1440 = edges[0];
				auto const& photon = edges[2];
				auto const g = std::any_cast<double>(n1440->particle()->user_data("gRNgamma"));
				auto const m_rho = Rho_Zero->mass();
				return -g/m_rho * dirac_sigmac(photon->four_momentum(kin), photon->lorentz_indices()[0]);
			}
	));
}