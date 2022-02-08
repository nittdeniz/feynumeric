#ifndef FEYNUMERIC_TOPOLOGY_HPP
#define FEYNUMERIC_TOPOLOGY_HPP

#include <map>
#include <vector>

#include "direction.hpp"

namespace Feynumeric
{
	class Topology{
	public:
		using Vertex_ID = std::size_t;
		struct Edge
		{
			Vertex_ID from, to;
			Direction direction;
		};
		using Adjacency_Map = std::map<Vertex_ID, std::map<Vertex_ID, std::vector<std::size_t>>>;
	private:
		Adjacency_Map _adjacency_map;
		std::vector<Edge> _edge_list;
		std::vector<std::size_t> _incoming_edges, _outgoing_edges, _virtual_edges;

		void create_adjacency_map();
		void validate() const;
	public:
		Topology(std::vector<Edge> const& edge_list);

		Topology(Topology const& copy);
		Topology& operator=(Topology const& copy);

		friend class Feynman_Graph;
		friend class Feynman_Diagram;
	};
}
#endif // FEYNUMERIC_TOPOLOGY_HPP