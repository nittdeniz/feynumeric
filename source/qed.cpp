#include "feynman_graph.hpp"
#include "dirac.hpp"
#include "matrix.hpp"
#include "particle.hpp"
#include "qed.hpp"
#include "vertex.hpp"
#include "units.hpp"

namespace Feynumeric
{
	namespace QED
	{
		using namespace Units;
		Particle_Ptr Photon = std::make_shared<Particle>("Photon", Particle::Type::NeutralBoson, 0, 0, 0, 1);

		Particle_Ptr Electron = std::make_shared<Particle>("Electron", Particle::Type::TrueLepton, 0.5109989461_MeV, 0, -1,
		                                                   0.5);
		Particle_Ptr Positron = std::make_shared<Particle>("Positron", Particle::Type::AntiLepton, 0.5109989461_MeV, 0,
		                                                   1, 0.5);

		Particle_Ptr Muon_Plus = std::make_shared<Particle>("Muon_+", Particle::Type::AntiLepton, 105.6583745_MeV, 0, 1,
		                                                    0.5);
		Particle_Ptr Muon_Minus = std::make_shared<Particle>("Muon_-", Particle::Type::TrueLepton, 105.6583745_MeV, 0, -1,
		                                                     0.5);


		void init_particles()
		{
			Photon->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
			Photon->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return epsilon(e, kin); };
			Photon->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return epsilon_star(e, kin); };

			Electron->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
			Electron->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return u(e, kin); };
			Electron->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return ubar(e, kin); };

			Positron->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
			Positron->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return vbar(e, kin); };
			Positron->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return v(e, kin); };

			Muon_Minus->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
			Muon_Minus->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return u(e, kin); };
			Muon_Minus->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return ubar(e, kin); };

			Muon_Plus->feynman_virtual = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return Propagator(e, kin); };
			Muon_Plus->feynman_incoming = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return vbar(e, kin); };
			Muon_Plus->feynman_outgoing = [](Feynman_Graph::Edge_Ptr e, Kinematics const& kin){ return v(e, kin); };
		}

		Feynumeric::Vertex_Manager_Ptr VMP = std::make_shared<Feynumeric::Vertex_Manager>();

		void init_vertices()
		{
			/*
			using Feynumeric::Direction;
			using namespace Feynumeric::Units;
			VMP->add(Feynumeric::Vertex(
					{
							{Electron, Direction::ANY},
							{Electron, Direction::ANY},
							{Photon,   Direction::ANY}
					},
					[](Feynumeric::Kinematics const&, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
						auto const& photon = edges[2];
						return 1._e * Feynumeric::GAC[*( photon->lorentz_indices()[0] )];
					}
			));
			VMP->add(Feynumeric::Vertex(
					{
							{Electron, Direction::ANY},
							{Positron, Direction::ANY},
							{Photon,   Direction::ANY}
					},
					[](Feynumeric::Kinematics const&, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
						auto const& photon = edges[2];
						return 1._e * Feynumeric::GAC[*( photon->lorentz_indices()[0] )];
					}
			));
			VMP->add(Feynumeric::Vertex(
					{
							{Positron, Direction::ANY},
							{Positron, Direction::ANY},
							{Photon,   Direction::ANY}
					},
					[](Feynumeric::Kinematics const&, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
						auto const& photon = edges[2];
						return 1._e * Feynumeric::GAC[*( photon->lorentz_indices()[0] )];
					}
			));


			VMP->add(Feynumeric::Vertex(
					{
							{Muon_Minus, Direction::ANY},
							{Muon_Minus, Direction::ANY},
							{Photon,     Direction::ANY}
					},
					[](Feynumeric::Kinematics const&, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
						auto const& photon = edges[2];
						return 1._e * Feynumeric::GAC[*( photon->lorentz_indices()[0] )];
					})
			);
			VMP->add(Feynumeric::Vertex(
					{
							{Muon_Minus, Direction::ANY},
							{Muon_Plus,  Direction::ANY},
							{Photon,     Direction::ANY}
					},
					[](Feynumeric::Kinematics const&, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
						auto const& photon = edges[2];
						return 1._e * Feynumeric::GAC[*( photon->lorentz_indices()[0] )];
					})
			);
			 */
		}
	}
}
