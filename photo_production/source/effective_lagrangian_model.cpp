#include <feynumeric/dirac.hpp>
#include <feynumeric/particle.hpp>
#include <feynumeric/qed.hpp>
#include <feynumeric/units.hpp>


#include "effective_lagrangian_model.hpp"

using namespace Feynumeric::Units;
using Feynumeric::Particle;
using Feynumeric::Matrix;

Feynumeric::Particle_Ptr Photon   = std::make_shared<Feynumeric::Particle>("Photon", Particle::Type::Majorana, 0, 0, 0, 1);

Feynumeric::Particle_Ptr Proton   = std::make_shared<Feynumeric::Particle>("Proton", Particle::Type::Particle, 938.270813_MeV, 0, +1, 0.5);
Feynumeric::Particle_Ptr Neutron  = std::make_shared<Feynumeric::Particle>("Neutron", Particle::Type::Particle, 939.5654133_MeV, 0, 0, 0.5);

Feynumeric::Particle_Ptr Pi_Zero  = std::make_shared<Feynumeric::Particle>("Pi_0", Particle::Type::Majorana, 134.9768_MeV, 0, 0, 0);
Feynumeric::Particle_Ptr Pi_Plus  = std::make_shared<Feynumeric::Particle>("Pi_+", Particle::Type::AntiParticle, 139.57039_MeV, 0, +1, 0);
Feynumeric::Particle_Ptr Pi_Minus = std::make_shared<Feynumeric::Particle>("Pi_-", Particle::Type::Particle, 139.57039_MeV, 0, -1, 0);

Feynumeric::Particle_Ptr Rho_Zero  = std::make_shared<Feynumeric::Particle>("Rho_0", Particle::Type::Majorana, 775.26_MeV, 147.8_MeV, 0, 1);
Feynumeric::Particle_Ptr Rho_Plus  = std::make_shared<Feynumeric::Particle>("Rho_+", Particle::Type::AntiParticle, 775.11_MeV, 149.1_MeV, +1, 1);
Feynumeric::Particle_Ptr Rho_Minus = std::make_shared<Feynumeric::Particle>("Rho_-", Particle::Type::Particle, 775.11_MeV, 149.1_MeV, -1, 1);

Feynumeric::Particle_Ptr N1440p   = std::make_shared<Feynumeric::Particle>("N1440_+", Particle::Type::Particle, 1.44_GeV, 350._MeV, +1, 0.5);
Feynumeric::Particle_Ptr N1440n   = std::make_shared<Feynumeric::Particle>("N1440_0", Particle::Type::Particle, 1.44_GeV, 350._MeV, 0, 0.5);


void init_particles()
{
    using namespace Feynumeric;
    //QED::init_vertices();
    QED::init_particles();
    Proton->feynman_virtual  = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
    Proton->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return u(e, kin); };
    Proton->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return ubar(e, kin); };


//    Pi_Zero->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){return Matrix(1,1,1);};
    Pi_Zero->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){return Matrix(1,1,1);};
    Pi_Zero->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){return Matrix(1,1,1);};

    N1440p->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
    N1440p->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return u(e, kin); };
    N1440p->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return ubar(e, kin); };
    N1440p->user_data("gRNpi", 1.);
    N1440p->width([&](double p2){ return N1440p->width();});
    N1440p->parity(1);
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