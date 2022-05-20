#ifndef FEYNUMERIC_TOPOLOGY_HPP
#define FEYNUMERIC_TOPOLOGY_HPP

#include <map>
#include <vector>

#include "direction.hpp"

namespace Feynumeric
{
    enum class Type{
        INCOMING,
        VIRTUAL,
        OUTGOING
    };
    class Topology_Vertex{
    public:
    private:
        Type _type;
        std::string _label;
        std::size_t _id;
        std::size_t _num;
    public:
        Topology_Vertex(std::string const& label);
        Topology_Vertex(Topology_Vertex const& vertex);
        Topology_Vertex& operator=(Topology_Vertex const& vertex);
        std::size_t id() const;
        void id(std::size_t i);
        std::size_t num() const;
        std::string label() const;
        Type type() const;
        friend class Topology_Edge;
        friend class Topology;
    };
    class Topology_Edge
    {
    private:
        Topology_Vertex _from, _to;
        Type _type;
    public:
        Topology_Edge(std::string const& from, std::string const& to);
        Topology_Edge(Topology_Vertex const& f, Topology_Vertex const& t);
        Topology_Edge(Topology_Edge const& copy);
        Topology_Edge& operator=(Topology_Edge const& copy);
        Type type() const;
        Topology_Vertex& from();
        Topology_Vertex& to();

    };
	class Topology{
	public:
	private:
	    std::vector<Topology_Edge> _edges;
	    std::vector<std::size_t> _incoming_edges;
	    std::vector<std::size_t> _outgoing_edges;
	    std::vector<std::size_t> _virtual_edges;
	    std::map<std::size_t, std::map<std::size_t, std::vector<std::size_t>>> _adjacency_map;
		void validate();
	public:
		Topology(std::vector<Topology_Edge> edge_list);

		Topology(Topology const& copy);
		Topology& operator=(Topology const& copy);

		friend class Feynman_Graph;
		friend class Feynman_Diagram;
	};
}
#endif // FEYNUMERIC_TOPOLOGY_HPP