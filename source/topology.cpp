#include "format.hpp"
#include "messages.hpp"
#include "topology.hpp"


namespace Feynumeric
{
	Topology::Topology(std::vector<Edge> const& edge_list)
	: _edge_list(edge_list)
	{
		create_adjacency_map();
		validate();
	}

	Topology::Topology(Topology const& copy)
	: _adjacency_map(copy._adjacency_map)
	, _edge_list(copy._edge_list)
	, _incoming_edges(copy._incoming_edges)
	, _outgoing_edges(copy._outgoing_edges)
	, _virtual_edges(copy._virtual_edges)
	{
	}

	Topology& Topology::operator=(Topology const& copy){
		_adjacency_map = copy._adjacency_map;
		_edge_list = copy._edge_list;
		_incoming_edges = copy._incoming_edges;
		_outgoing_edges = copy._outgoing_edges;
		_virtual_edges = copy._outgoing_edges;
		return *this;
	}

	void Topology::validate() const{
		for( auto const& edge : _edge_list )
		{
			if( edge.direction == Direction::VIRTUAL )
			{
				if( !(_adjacency_map.at(edge.from).size() > 1 && _adjacency_map.at(edge.to).size() > 1) )
				{
					critical_error(FORMAT("Edge({}, {}) does not appear to be virtual.", edge.from, edge.to));
				}
			}
			else
			{
				if( (_adjacency_map.at(edge.from).size() > 1 && _adjacency_map.at(edge.to).size() > 1) )
				{
					critical_error(FORMAT("Edge({}, {}) does not appear to be external.", edge.from, edge.to));
				}
			}
		}
	}

	void Topology::create_adjacency_map(){
		std::size_t k = 0;
		for( auto const& edge : _edge_list )
		{
			_adjacency_map[edge.from][edge.to].push_back(k);
			_adjacency_map[edge.to][edge.from].push_back(k);
			switch( edge.direction )
			{
				case Direction::INCOMING:
					_incoming_edges.push_back(k);
					break;
				case Direction::OUTGOING:
					_outgoing_edges.push_back(k);
					break;
				case Direction::VIRTUAL:
					_virtual_edges.push_back(k);
					break;
				default:
					critical_error("Invalid direction in topology.");
			}
			k++;
		}
	}

}