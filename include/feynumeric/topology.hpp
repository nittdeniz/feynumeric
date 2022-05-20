#ifndef FEYNUMERIC_TOPOLOGY_HPP
#define FEYNUMERIC_TOPOLOGY_HPP

#include <map>
#include <vector>

#include "direction.hpp"

namespace Feynumeric
{
    struct Vertex{
        enum class Type{
            INCOMING,
            VIRTUAL,
            OUTGOING
        };
        using ID = std::size_t;
        Type type;
        ID id;
        Vertex(std::string const& id);
        Vertex(Type t, ID i);
        Vertex(Vertex const& vertex);
        Vertex& operator=(Vertex const& vertex);
    };
    struct Edge
    {
        Vertex from, to;
        Edge(Vertex const& f, Vertex const& t);
        Edge(std::string const& from, std::string const& to);
    };
	class Topology{
	public:
	private:
	    std::vector<Edge> _edges;
	    std::vector<std::size_t> _incoming_edges;
	    std::vector<std::size_t> _outgoing_edges;
	    std::vector<std::size_t> _virtual_edges;
	    std::map<Vertex::ID, std::map<Vertex::ID, std::vector<std::size_t>>> _adjacency_map;
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