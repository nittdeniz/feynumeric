#include "format.hpp"
#include "messages.hpp"
#include "topology.hpp"

namespace Feynumeric{
    Vertex::Vertex(std::string const& id)
    {
        if( id.size() < 2 ){
            critical_error("Topology::Vertex must be at least 2 characters long.");
        }
        switch( id[0] ){
            case 'i':{
                type = Type::INCOMING;
                break;
            }
            case 'o':{
                type = Type::OUTGOING;
                break;
            }
            case 'v':{
                type = Type::VIRTUAL;
                break;
            }
            default:{
                critical_error(FORMAT("Topology::Vertex first character must be 'i', 'o', or 'v'. Given: {}.", id[0]));
            }
        }
        id = std::stoull(id.substr(1));
    }

    Vertex::Vertex(Vertex::Type t, Vertex::ID i)
    : type(t)
    , id(i)
    {
    }

    Vertex::Vertex(const Vertex &vertex)
    : type(vertex.type)
    , id(vertex.id)
    {

    }

    Vertex& Vertex::operator=(const Vertex &vertex)
    {
        type = vertex.type;
        id = vertex.id;
        return *this;
    }

    Edge::Edge(const Vertex &f, const Vertex &t)
    : from(f)
    , to(t)
    {

    }

    Edge::Edge(std::string const& from, std::string const& to)
    : from(Vertex(from))
    , to(Vertex(to))
    {

    }

    void Topology::validate() const
    {
        if( _edges.empty() ){
            critical_error("Topology contains no edges.");
        }
        if( _incoming_edges.empty() ){
            critical_error("Topology contains no incoming edges.");
        }
        if( _outgoing_edges.empty() ){
            critical_error("Topology contains no outgoing edges.");
        }
    }

    Topology::Topology(const std::vector<Edge> &edge_list)
    {
        for( auto const& edge : edge_list ){
            auto pos = _edges.size();
            if( edge.from.type == Vertex::Type::INCOMING && edge.to.type == Vertex::Type::VIRTUAL ){
                _incoming_edges.push_back(pos);
            }
            else if( edge.from.type == Vertex::Type::VIRTUAL && edge.to.type == Vertex::Type::OUTGOING ){
                _outgoing_edges.push_back(pos);
            }
            else if( edge.from.type == Vertex::Type::VIRTUAL && edge.to.type == Vertex::Type::VIRTUAL ){
                _virtual_edges.push_back(pos);
            }
            else{
                critical_error("Topology::Edge must be pair (i, v), (v, v), or (v, o).");
            }
            _adjacency_map[edge.from.id][edge.to.id].push_back(pos);
            _edges.push_back(edge);
        }
    }

    Topology::Topology(const Topology &copy)
    {

    }

    Topology &Topology::operator=(const Topology &copy)
    {
    }
}