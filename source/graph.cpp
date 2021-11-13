#include "graph.hpp"

namespace Feyncalc
{
    Graph::Graph(vector<int> &&incoming_vertices, vector<int> &&virtual_vertices, vector<int> &&outgoing_vertices, vector<array<int, 2>>&& edges)
    : _incoming_vertices(std::move(incoming_vertices))
    , _virtual_vertices(std::move(virtual_vertices))
    , _outgoing_vertices(std::move(outgoing_vertices))
    {
        for( auto edge : edges )
        {
            _adjacency_map[edge[0]][edge[1]];
            _adjacency_map[edge[1]][edge[0]];
        }
    }

    Graph& Graph::add_edge(int v1, int v2, Particle_Ptr particle)
    {
        _adjacency_map[v1][v2].push_back(particle);
        _adjacency_map[v2][v1].push_back(particle);
        return *this;
    }

}

