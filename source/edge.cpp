#include <iostream>
#include "messages.hpp"
#include "edge.hpp"
#include "particle.hpp"

namespace Feynumeric
{
    Edge::Edge(Vertex_Id a, Vertex_Id b, Edge::Type type, Particle_Ptr particle)
    : _a(a)
    , _b(b)
    , _type(type)
    , _particle(particle)
    {
    }

    Edge::Edge(const Edge &edge)
    : _a(edge._a)
    , _b(edge._b)
    , _type(edge._type)
    , _particle(edge._particle)
    , _momentum(edge._momentum)
    , _neighbour_ids(edge._neighbour_ids)
    {

    }

    Edge &Edge::operator=(const Edge &edge)
    {
        _a = edge._a;
        _b = edge._b;
        _type = edge._type;
        _particle = edge._particle;
        _momentum = edge._momentum;
        _neighbour_ids = edge._neighbour_ids;
        return *this;
    }

    Vertex_Id Edge::a() const
    {
        return _a;
    }

    Vertex_Id Edge::b() const
    {
        return _b;
    }

    bool Edge::is_incoming() const
    {
        return _type == Type::INCOMING;
    }

    bool Edge::is_outgoing() const
    {
        return _type == Type::OUTGOING;
    }

    bool Edge::is_virtual() const
    {
        return _type == Type::VIRTUAL;
    }

    bool Edge::is_undefined() const
    {
        return _type == Type::UNDEFINED;
    }

    void Edge::add_neighbour(Edge_Id neighbour)
    {
        _neighbour_ids.push_back(neighbour);
    }

    vector<Edge_Id> Edge::neighbour_ids() const
    {
        return _neighbour_ids;
    }

    Matrix Edge::momentum() const
    {
        return _momentum;
    }

    std::ostream& operator<<(std::ostream &out, const Edge &edge)
    {
        return out << "(" << static_cast<std::size_t>(edge._a) << ", " << static_cast<std::size_t>(edge._b) << ")";
    }

    std::vector<std::size_t> Edge::get_lorentz_indices() const
    {
        return _assigned_indices;
    }

    std::function<Matrix()> Edge::feynman_rule() const
    {
        if( is_incoming() )
        {
            return _particle->feynman_incoming;
        }
        else if( is_outgoing() )
        {
            return _particle->feynman_outgoing;
        }
        else if( is_virtual() )
        {
            return _particle->feynman_virtual;
        }
        else
        {
            critical_error("Can't assign feynman rule to edge: " + to_string() + ".");
        }
    }

    Edge::Edge(std::size_t a, std::size_t b, Edge::Type type)
    : _a(Vertex_Id{a})
    , _b(Vertex_Id{b})
    , _type(type)
    {
        if( _type == Type::UNDEFINED )
        {
            warning("Edge Type specified as UNDEFINED at " + to_string() + ".");
        }
    }

    void Edge::momentum(const Matrix &momentum)
    {
        _momentum = momentum;
    }

    bool operator==(const Edge &lhs, const Edge &rhs)
    {
        return lhs._particle == rhs._particle &&
               lhs._momentum == rhs._momentum &&
                (
                    ( lhs._a == rhs._a && lhs._b == rhs._b ) ||
                    ( lhs._b == rhs._a && lhs._a == rhs._b )
                );
    }

    bool operator!=(const Edge &lhs, const Edge &rhs)
    {
        return !(lhs==rhs);
    }

    bool shares_vertex(const Edge &lhs, const Edge &rhs)
    {
        return lhs._a == rhs._a || lhs._a == rhs._b || lhs._b == rhs._a || lhs._b == rhs._b;
    }

    std::string Edge::to_string() const
    {
        std::stringstream in;
        in << *this;
        return in.str();
    }

    void Edge::assign_lorentz_index(std::size_t id)
    {
        _assigned_indices.emplace_back(id);
    }

    void Edge::clear_lorentz_indices()
    {
        _assigned_indices.clear();
    }


    void Edge::assign_angular_momentum(std::size_t id)
    {
        _angular_momentum_index = id;
    }

    std::optional<Vertex_Id> shared_vertex(Edge const& lhs, Edge const& rhs)
    {
        if( lhs._a == rhs._a || lhs._a == rhs._b )
        {
            return lhs._a;
        }
        if( lhs._b == rhs._a || lhs._b == rhs._b )
        {
            return lhs._b;
        }
        return std::nullopt;
    }

    void Edge::particle(Particle_Ptr new_particle)
    {
        _particle = new_particle;
    }

    Particle_Ptr Edge::particle() const
    {
        return _particle;
    }
}