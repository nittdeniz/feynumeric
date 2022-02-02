#include "graph.hpp"
#include "utility.hpp"

namespace Feynumeric
{
	Graph::Graph(const Graph& other)
	: _edges(other._edges)
	{
		initialize_edges();
	}

	Graph::Graph(std::list<Edge>&& edges)
    : _edges(std::move(edges))
    {
		initialize_edges();
    }

	Graph& Graph::operator=(const Graph& other){
		_edges = other._edges;
		initialize_edges();
		return *this;
	}

	void Graph::initialize_edges()
    {
		for( auto& edge : _edges )
		{
			edge.clear_neighbours();
		}
		for( auto edge_it = _edges.begin(); edge_it != _edges.end(); ++edge_it )
	    {
			auto* edge_ptr = &(*edge_it);
		    if( edge_ptr->is_incoming() )
		    {
			    _incoming_edges.push_back(edge_ptr);
		    }
		    else if( edge_ptr->is_outgoing() )
		    {
			    _outgoing_edges.push_back(edge_ptr);
		    }
		    else if( edge_ptr->is_virtual() )
		    {
			    _virtual_edges.push_back(edge_ptr);
		    }
		    else
		    {
			    critical_error("Undefined edge (needs to be incoming, outgoing, or virtual): " + edge_ptr->to_string());
		    }
		    auto temp = edge_it;
		    std::advance(temp, 1);
		    for( auto edge_jt = temp; edge_jt != _edges.end(); ++edge_jt )
		    {
		    	auto* edge2_ptr = &(*edge_jt);
			    if( shares_vertex(edge_ptr, edge2_ptr) )
			    {
				    edge_ptr->add_neighbour(edge2_ptr);
				    edge2_ptr->add_neighbour(edge_ptr);
			    }
		    }
	    }
    }

    std::vector<Edge*> Graph::all_edges()
    {
        std::vector<Edge*> edges;
        edges.reserve(_edges.size());
        for( auto& edge : _edges )
        {
        	edges.emplace_back(&edge);
        }
        return edges;
    }

    std::map<std::size_t, std::vector<Edge*>> Graph::all_vertex_ids()
    {
        std::map<std::size_t, std::vector<Edge*>> map;
        for( auto& edge : _edges )
        {
        	map[edge._a].emplace_back(&edge);
        	map[edge._b].emplace_back(&edge);
        }

        std::erase_if(map, [](auto const& item){
        	auto const& [key, value] = item;
        	return value.size() <= 1;
        });

        return map;
    }

    std::vector<Edge*> Graph::edges_to(std::size_t vertex_id)
    {
        std::vector<Edge*> edges;
        for( auto& edge : _edges )
        {
        	if( edge._a == vertex_id || edge._b == vertex_id )
	        {
        		edges.push_back(&edge);
	        }
        }
        return edges;
    }
}

