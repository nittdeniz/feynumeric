#ifndef FEYNCALC_GRAPH_HPP
#define FEYNCALC_GRAPH_HPP

#include <array>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "complex.hpp"
#include "edge.hpp"
#include "kinematics.hpp"
#include "particle.hpp"
#include "vertex.hpp"

namespace Feyncalc
{
    class Graph
    {
    private:
        std::vector<Edge> _edges;
        std::vector<Edge_Id> _incoming_edge_ids;
        std::vector<Edge_Id> _outgoing_edge_ids;
        std::vector<Edge_Id> _virtual_edge_ids;

    public:
        Graph(std::vector<Edge>&& edges);

        std::vector<Edge_Id> all_edge_ids() const;
        std::map<Vertex_Id, std::vector<Edge_Id>> all_vertex_ids() const;
        std::vector<Edge_Id> edge_ids_to(Vertex_Id vertex) const;
        friend class Diagram;
    };

    using Graph_Ptr = std::shared_ptr<Graph>;
}

#endif // FEYNCALC_GRAPH_HPP