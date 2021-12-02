#include "particles.hpp"
#include "vertices.hpp"

#include <feyncalc/edge.hpp>

Feyncalc::Vertex_Manager_Ptr VMP = std::make_shared<Feyncalc::Vertex_Manager>();

void init_vertices()
{
    using Direction = Feyncalc::Vertex_Manager::Direction;
    VMP->add_vertex(
            [](Feyncalc::Diagram* diagram, std::vector<Feyncalc::Edge> const& edges)
            {
                auto const& electron_in = edges[0];
                auto const& electron_out = edges[1];
                auto const& photon = edges[2];
                return Feyncalc::Matrix(4, 4, 7);
            },
            {
                {Electron, Direction::IN}
                ,{Electron, Direction::OUT}
                ,{Photon, Direction::BOTH}
            }
            );
    VMP->add_vertex(
            [](Feyncalc::Diagram* diagram, std::vector<Feyncalc::Edge> const& edges)
            {
                auto const& electron_in = edges[0];
                auto const& electron_out = edges[1];
                auto const& photon = edges[2];
                return Feyncalc::Matrix(4, 4, 11);
            },
            {
                    {Muon_Minus, Direction::IN}
                    ,{Muon_Minus, Direction::OUT}
                    ,{Photon, Direction::BOTH}
            }
    );
}

