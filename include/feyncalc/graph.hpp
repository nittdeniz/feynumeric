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

namespace Feyncalc
{
    using std::array;
    using std::function;
    using std::map;
    using std::vector;
    using std::cerr;
    using std::abort;
    class Graph
    {
    private:
        vector<Edge> _edges;
        vector<size_t> _incoming_edges;
        vector<size_t> _outgoing_edges;
        vector<size_t> _virtual_edges;

    public:
        Graph(vector<Edge>&& edges);

        std::vector<Edge> all_edges() const;

        friend class Diagram;
    };

    using Graph_Ptr = std::shared_ptr<Graph>;
}

#endif // FEYNCALC_GRAPH_HPP