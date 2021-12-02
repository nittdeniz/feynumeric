#include "graph.hpp"
#include "utility.hpp"

namespace Feyncalc
{
    Graph::Graph(vector<Edge>&& edges)
    : _edges(std::move(edges))
    {
        for( auto edge_it = _edges.begin(); edge_it != _edges.end(); ++edge_it ){
            if( edge_it->is_incoming()){
                _incoming_edge_ids.push_back(Edge_Id(edge_it - _edges.begin()));
            }
            else if( edge_it->is_outgoing()){
                _outgoing_edge_ids.push_back(Edge_Id(edge_it - _edges.begin()));
            }
            else if( edge_it->is_virtual()){
                _virtual_edge_ids.push_back(Edge_Id(edge_it - _edges.begin()));
            }
            else{
                std::cerr << "Undefined edge (needs to be incoming, outgoing or virtual): " << *edge_it << "\n";
                std::abort();
            }
            for( auto edge_jt = edge_it+1; edge_jt != _edges.end(); edge_jt ++ )
            {
                if( shares_vertex(*edge_it, *edge_jt) )
                {
                    edge_it->add_neighbour(Edge_Id(edge_jt - _edges.begin()));
                    edge_jt->add_neighbour(Edge_Id(edge_it - _edges.begin()));
                }
            }
        }
    }

    std::vector<Edge_Id> Graph::all_edge_ids() const
    {
        std::vector<Edge_Id> edge_ids;
        edge_ids.reserve(_edges.size());
        for( std::size_t i = 0; i < _edges.size(); ++i )
        {
            edge_ids.emplace_back(i);
        }
        return edge_ids;
    }

    std::map<Vertex_Id, std::vector<Edge_Id>> Graph::all_vertex_ids() const
    {
        std::map<Vertex_Id, std::vector<Edge_Id>> map;
        for( std::size_t i = 0; i < _edges.size(); ++i )
        {
            auto const& edge = _edges[i];
            map[Vertex_Id(edge._a)].emplace_back(i);
            map[Vertex_Id(edge._b)].emplace_back(i);
        }

        for(auto it = map.begin(); it != map.end(); ) {
            if( it->second.size() <= 1 )
            {
                it = map.erase(it);
            }
            else
            {
                ++it;
            }
        }
        return map;
    }

    std::vector<Edge_Id> Graph::edge_ids_to(Vertex_Id vertex) const
    {
        std::vector<Edge_Id> edges;
        for( std::size_t i = 0; i < _edges.size();  ++i )
        {
            auto const& edge = _edges[i];
            if( edge._a == vertex || edge._b == vertex )
            {
                edges.emplace_back(i);
            }
        }
        return edges;
    }
}

