#include "graph.hpp"

namespace Feyncalc
{
    Graph::Graph(vector<int> &&incoming_vertices, vector<int> &&virtual_vertices, vector<int> &&outgoing_vertices, vector<vector<int>>&& edges)
    : _incoming_vertices(std::move(incoming_vertices))
    , _virtual_vertices(std::move(virtual_vertices))
    , _outgoing_vertices(std::move(outgoing_vertices))
    {
        edges.size();
    }

    Graph& Graph::add_edge(int v1, int v2, Particle_Ptr particle)
    {
        _adjacency_map[v1][v2].push_back(particle);
        _adjacency_map[v2][v1].push_back(particle);
        return *this;
    }

    void Graph::generate_amplitude()
    {

    }

    Complex Graph::amplitude(const Kinematics &kinematics)
    {
        kinematics.f();
        return Feyncalc::Complex();
    }
}

