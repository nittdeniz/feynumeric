#include "particles.hpp"
#include "vertices.hpp"

#include <Feynumeric/edge.hpp>

Feynumeric::Vertex_Manager VM = Feynumeric::Vertex_Manager();

void init_vertices()
{
    using Direction = Feynumeric::Vertex_Manager::Direction;
    VM.add_vertex(
            [](Feynumeric::Diagram* diagram, std::vector<Feynumeric::Edge*> const& edges)
            {
                auto const& electron_in = edges[0];
                auto const& electron_out = edges[1];
                auto const& photon = edges[2];
                return Feynumeric::Matrix(4, 4, 7);
            },
            {
                {Electron, Direction::IN}
                ,{Electron, Direction::OUT}
                ,{Photon, Direction::BOTH}
            }
            );
    VM.add_vertex(
            [](Feynumeric::Diagram* diagram, std::vector<Feynumeric::Edge*> const& edges)
            {
                auto const& electron_in = edges[0];
                auto const& electron_out = edges[1];
                auto const& photon = edges[2];
                return Feynumeric::Matrix(4, 4, 11);
            },
            {
                    {Muon_Minus, Direction::IN}
                    ,{Muon_Minus, Direction::OUT}
                    ,{Photon, Direction::BOTH}
            }
    );
}

