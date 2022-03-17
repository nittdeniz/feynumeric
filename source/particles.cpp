#include "dirac.hpp"
#include "units.hpp"
#include "particles.hpp"

namespace Feynumeric{
	namespace QED
	{
		using namespace Units;
		Particle_Ptr Photon       = std::make_shared<Particle>("Photon", Particle::Type::NeutralBoson, 0, 0, 1);

		Particle_Ptr Electron     = std::make_shared<Particle>("Electron", Particle::Type::TrueParticle, 0.5109989461_MeV, -1, 0.5);
		Particle_Ptr Positron     = std::make_shared<Particle>("Positron", Particle::Type::AntiParticle, 0.5109989461_MeV, 1, 0.5);

		Particle_Ptr Muon_Plus    = std::make_shared<Particle>("Muon_+", Particle::Type::AntiParticle, 105.6583745_MeV, 1, 0.5);
		Particle_Ptr Muon_Minus   = std::make_shared<Particle>("Muon_-", Particle::Type::TrueParticle, 105.6583745_MeV, -1, 0.5);


		void init_particles()
		{
			Photon->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Projector(e, kin); };
			Photon->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return epsilon(e, kin); };
			Photon->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){
				return epsilon_star(e, kin);
			};

			Electron->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){
				return Propagator(e, kin);
			};
			Electron->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return u(e, kin); };
			Electron->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return ubar(e, kin); };

			Positron->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){
				return Propagator(e, kin);
			};
			Positron->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return vbar(e, kin); };
			Positron->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return v(e, kin); };

			Muon_Minus->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){
				return Propagator(e, kin);
			};
			Muon_Minus->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return u(e, kin); };
			Muon_Minus->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return ubar(e, kin); };

			Muon_Plus->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){
				return Propagator(e, kin);
			};
			Muon_Plus->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return vbar(e, kin); };
			Muon_Plus->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return v(e, kin); };
		}
	}
}