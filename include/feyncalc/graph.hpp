#ifndef FEYNCALC_GRAPH_HPP
#define FEYNCALC_GRAPH_HPP

#include <map>
#include <memory>
#include <vector>

#include "complex.hpp"
#include "kinematics.hpp"

#include "particle.hpp"

namespace Feyncalc
{
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
        vector<function<Complex(Kinematics const&)>> _amplitude;
    public:
        Graph(vector<int>&& incoming_vertices, vector<int>&& virtual_vertices, vector<int>&& outgoing_vertices, vector<vector<int>>&& edges);
        Graph& add_edge(int v1, int v2, Particle_Ptr particle);
        void generate_amplitude();
        Complex amplitude(Kinematics const& kinematics);


    };

    using Graph_Ptr = std::shared_ptr<Graph>;
}

#endif // FEYNCALC_GRAPH_HPP