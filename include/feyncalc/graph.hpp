#ifndef FEYNCALC_GRAPH_HPP
#define FEYNCALC_GRAPH_HPP

#include <array>
#include <map>
#include <memory>
#include <vector>

#include "complex.hpp"
#include "kinematics.hpp"

#include "particle.hpp"

namespace Feyncalc
{
    using std::array;
    using std::function;
    using std::map;
    using std::vector;
    class Graph
    {
    private:
        vector<int> _incoming_vertices;
        vector<int> _virtual_vertices;
        vector<int> _outgoing_vertices;
        map<int, map<int, vector<Particle_Ptr>>> _adjacency_map;
    public:
        Graph(vector<int>&& incoming_vertices, vector<int>&& virtual_vertices, vector<int>&& outgoing_vertices, vector<array<int, 2>>&& edges);
        Graph& add_edge(int v1, int v2, Particle_Ptr particle);

        friend class Diagram;
    };

    using Graph_Ptr = std::shared_ptr<Graph>;
}

#endif // FEYNCALC_GRAPH_HPP