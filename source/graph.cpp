#include "graph.hpp"
#include "utility.hpp"

namespace Feyncalc
{
    Graph::Graph(vector<Edge>&& edges)
    : _edges(std::move(edges))
    {
        for( auto edge_it = _edges.begin(); edge_it != _edges.end(); ++edge_it ){
            if( edge_it->is_incoming()){
                _incoming_edges.push_back(edge_it - _edges.begin());
            }
            else if( edge_it->is_outgoing()){
                _outgoing_edges.push_back(edge_it - _edges.begin());
            }
            else if( edge_it->is_virtual()){
                _virtual_edges.push_back(edge_it - _edges.begin());
            }
            else{
                cerr << "Undefined edge (needs to be incoming, outgoing or virtual): " << *edge_it << "\n";
                abort();
            }
            for( auto edge_jt = edge_it+1; edge_jt != _edges.end(); edge_jt ++ )
            {
                if( shares_vertex(*edge_it, *edge_jt) )
                {
                    edge_it->add_neighbour(edge_jt - _edges.begin());
                    edge_jt->add_neighbour(edge_it - _edges.begin());
                }
            }
        }
    }

    vector<Edge> Graph::all_edges() const
    {
        return _edges;
    }
}

