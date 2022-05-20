/**
#include "particles.hpp"
#include "vertices.hpp"
#include <feynumeric/feynman_diagram.hpp>
#include <feynumeric/direction.hpp>
#include <feynumeric/kinematics.hpp>
#include <feynumeric/vertex_manager.hpp>
#include <feynumeric/units.hpp>
#include <feynumeric/dirac.hpp>
#include <feynumeric/vertex.hpp>

Feynumeric::Vertex_Manager_Ptr VMP = std::make_shared<Feynumeric::Vertex_Manager>();

void init_vertices()
{
    using Feynumeric::Direction;
    using namespace Feynumeric::Units;
    VMP->add(Feynumeric::Topology_Vertex(
        {
				           {Electron, Direction::BOTH},
				           {Electron, Direction::BOTH},
				           {Photon, Direction::BOTH}
		           },
		           [](Feynumeric::Kinematics const&, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
			           auto const& electron_in = edges[0];
			           auto const& electron_out = edges[1];
			           auto const& photon = edges[2];
			           return 1._e * Feynumeric::GAC[*( photon->lorentz_indices()[0] )];
		           }
    ));
	VMP->add(Feynumeric::Topology_Vertex(
			{
					{Electron, Direction::BOTH},
					{Positron, Direction::BOTH},
					{Photon, Direction::BOTH}
			},
			[](Feynumeric::Kinematics const&, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges){
				auto const& electron_in = edges[0];
				auto const& electron_out = edges[1];
				auto const& photon = edges[2];
				return 1._e * Feynumeric::GAC[*( photon->lorentz_indices()[0] )];
			}
	));


    VMP->add(Feynumeric::Topology_Vertex(
		    {
				     {Muon_Minus, Direction::BOTH}
				    ,{Muon_Minus, Direction::BOTH}
				    ,{Photon, Direction::BOTH}
		    },
		    [](Feynumeric::Kinematics const&, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges)
            {
                auto const& muon_in  = edges[0];
                auto const& muon_out = edges[1];
                auto const& photon   = edges[2];
	            return 1._e * Feynumeric::GAC[*(photon->lorentz_indices()[0])];
            })
    );
	VMP->add(Feynumeric::Topology_Vertex(
			{
					{Muon_Minus, Direction::BOTH}
					,{Muon_Plus, Direction::BOTH}
					,{Photon, Direction::BOTH}
			},
			[](Feynumeric::Kinematics const&, std::vector<Feynumeric::Feynman_Graph::Edge_Ptr> const& edges)
			{
				auto const& muon_in  = edges[0];
				auto const& muon_out = edges[1];
				auto const& photon   = edges[2];
				return 1._e * Feynumeric::GAC[*(photon->lorentz_indices()[0])];
			})
	);
}
 **/