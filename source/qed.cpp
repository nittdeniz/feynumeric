#include "dirac.hpp"
#include "edge_direction.hpp"
#include "feynman_graph.hpp"
#include "graph_edge.hpp"
#include "matrix.hpp"
#include "particle.hpp"
#include "qed.hpp"
#include "units.hpp"
#include "vertex.hpp"

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
			Photon->feynman_virtual = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return Propagator(e, kin); };
			Photon->feynman_incoming = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return epsilon(e, kin); };
			Photon->feynman_outgoing = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return epsilon_star(e, kin); };

			Electron->feynman_virtual = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return Propagator(e, kin); };
			Electron->feynman_incoming = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return u(e, kin); };
			Electron->feynman_outgoing = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return ubar(e, kin); };

			Positron->feynman_virtual = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return Propagator(e, kin); };
			Positron->feynman_incoming = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return vbar(e, kin); };
			Positron->feynman_outgoing = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return v(e, kin); };

			Muon_Minus->feynman_virtual = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return Propagator(e, kin); };
			Muon_Minus->feynman_incoming = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return u(e, kin); };
			Muon_Minus->feynman_outgoing = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return ubar(e, kin); };

			Muon_Plus->feynman_virtual = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return Propagator(e, kin); };
			Muon_Plus->feynman_incoming = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return vbar(e, kin); };
			Muon_Plus->feynman_outgoing = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return v(e, kin); };
		}

		Feynumeric::Vertex_Manager_Ptr VMP = std::make_shared<Feynumeric::Vertex_Manager>();

		void init_vertices()
		{
			using Feynumeric::Direction;
			using namespace Feynumeric::Units;
			VMP->add(Feynumeric::Vertex(
					{
							{Electron, Edge_Direction::ANY},
							{Electron, Edge_Direction::ANY},
							{Photon,   Edge_Direction::ANY}
					},
					[](Feynumeric::Kinematics const&, std::vector<std::shared_ptr<Graph_Edge>> const& edges){
						auto const& photon = edges[2];
						return 1._e * Feynumeric::GAC[*( photon->lorentz_indices()[0] )];
					}
			));
			VMP->add(Feynumeric::Vertex(
					{
							{Electron, Edge_Direction::ANY},
							{Positron, Edge_Direction::ANY},
							{Photon,   Edge_Direction::ANY}
					},
					[](Feynumeric::Kinematics const&, std::vector<std::shared_ptr<Graph_Edge>> const& edges){
						auto const& photon = edges[2];
						return 1._e * Feynumeric::GAC[*( photon->lorentz_indices()[0] )];
					}
			));
			VMP->add(Feynumeric::Vertex(
					{
							{Positron, Edge_Direction::ANY},
							{Positron, Edge_Direction::ANY},
							{Photon,   Edge_Direction::ANY}
					},
					[](Feynumeric::Kinematics const&, std::vector<std::shared_ptr<Graph_Edge>> const& edges){
						auto const& photon = edges[2];
						return 1._e * Feynumeric::GAC[*( photon->lorentz_indices()[0] )];
					}
			));


			VMP->add(Feynumeric::Vertex(
					{
							{Muon_Minus, Edge_Direction::ANY},
							{Muon_Minus, Edge_Direction::ANY},
							{Photon,     Edge_Direction::ANY}
					},
					[](Feynumeric::Kinematics const&, std::vector<std::shared_ptr<Graph_Edge>> const& edges){
						auto const& photon = edges[2];
						return 1._e * Feynumeric::GAC[*( photon->lorentz_indices()[0] )];
					})
			);
			VMP->add(Feynumeric::Vertex(
					{
							{Muon_Minus, Edge_Direction::ANY},
							{Muon_Plus,  Edge_Direction::ANY},
							{Photon,     Edge_Direction::ANY}
					},
					[](Feynumeric::Kinematics const&, std::vector<std::shared_ptr<Graph_Edge>> const& edges){
						auto const& photon = edges[2];
						return 1._e * Feynumeric::GAC[*( photon->lorentz_indices()[0] )];
					})
			);
		}
	}
}
