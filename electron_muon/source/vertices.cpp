#include "particles.hpp"
#include "vertices.hpp"

#include <feynumeric/diagram.hpp>
#include <feynumeric/dirac.hpp>
#include <feynumeric/edge.hpp>
#include <feynumeric/constants.hpp>

Feynumeric::Vertex_Manager VM = Feynumeric::Vertex_Manager();

void init_vertices()
{
    using Direction = Feynumeric::Vertex_Manager::Direction;
    VM.add_vertex(
            {
                 {Electron, Direction::IN}
                ,{Electron, Direction::OUT}
                ,{Photon, Direction::BOTH}
            },
            [](Feynumeric::Diagram* diagram, std::vector<Feynumeric::Edge*> const& edges)
            {
	            auto const& electron_in  = edges[0];
	            auto const& electron_out = edges[1];
	            auto const& photon       = edges[2];
	            return Feynumeric::Constants::e * Feynumeric::GA[*(photon->get_lorentz_indices()[0])];
            }
            );
    VM.add_vertex(
		    {
				     {Muon_Minus, Direction::IN}
				    ,{Muon_Minus, Direction::OUT}
				    ,{Photon, Direction::BOTH}
		    },
            [](Feynumeric::Diagram* diagram, std::vector<Feynumeric::Edge*> const& edges)
            {
                auto const& muon_in  = edges[0];
                auto const& muon_out = edges[1];
                auto const& photon   = edges[2];
	            return Feynumeric::Constants::e * Feynumeric::GA[*(photon->get_lorentz_indices()[0])];
            }
    );
}

