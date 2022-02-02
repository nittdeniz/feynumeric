#ifndef Feynumeric_GRAPH_HPP
#define Feynumeric_GRAPH_HPP

#include <array>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <vector>

#include "complex.hpp"
#include "edge.hpp"
#include "kinematics.hpp"
#include "particle.hpp"
#include "vertex.hpp"

namespace Feynumeric
{
    class Graph
    {
    private:
        std::list<Edge> _edges;
        std::vector<Edge*> _incoming_edges;
        std::vector<Edge*> _outgoing_edges;
        std::vector<Edge*> _virtual_edges;

        void initialize_edges();

    public:
    	Graph(Graph const& other);
        Graph(std::list<Edge>&& edges);
        Graph& operator=(Graph const& other);

        std::vector<Edge*> all_edges();
        std::map<std::size_t, std::vector<Edge*>> all_vertex_ids();
        std::vector<Edge*> edges_to(std::size_t vertex_id);
        friend class Diagram;
    };

    using Graph_Ptr = std::shared_ptr<Graph>;
}

#endif // Feynumeric_GRAPH_HPP