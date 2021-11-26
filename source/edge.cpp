#include <iostream>
#include "edge.hpp"

namespace Feyncalc
{
    Edge::Edge(size_t a, size_t b, Edge::Type type, Particle_Ptr particle, vector<int> momentum)
    : _a(a)
    , _b(b)
    , _type(type)
    , _particle(particle)
    , _momentum(momentum)
    {

    }

    Edge::Edge(const Edge &edge)
    : _a(edge._a)
    , _b(edge._b)
    , _type(edge._type)
    , _particle(edge._particle)
    , _momentum(edge._momentum)
    , _neighbours(edge._neighbours)
    {

    }

    Edge &Edge::operator=(const Edge &edge)
    {
        _a = edge._a;
        _b = edge._b;
        _type = edge._type;
        _particle = edge._particle;
        _momentum = edge._momentum;
        _neighbours = edge._neighbours;
        return *this;
    }

    size_t Edge::a() const
    {
        return _a;
    }

    size_t Edge::b() const
    {
        return _b;
    }

    void Edge::a(size_t new_a)
    {
        _a = new_a;
    }

    void Edge::b(size_t new_b)
    {
        _b = new_b;
    }

    Edge Edge::sorted() const
    {
        return Edge(std::min(_a, _b), std::max(_a, _b), _type, _particle, _momentum);
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

    Matrix Edge::momentum(const vector<Matrix> &momenta) const
    {
        if( momenta.size() != _momentum.size() )
        {
            std::cerr << "Momenta size does not match size of edge_momenta vector.\n";
            abort();
        }
        Matrix result;
        for( size_t i = 0; i < _momentum.size(); ++i )
        {
            result += _momentum.at(i) * momenta.at(i);
        }
        return result;
    }

    void Edge::type(Edge::Type new_type)
    {
        _type = new_type;
    }

    void Edge::add_neighbour(size_t neighbour)
    {
        _neighbours.push_back(neighbour);
    }

    vector<size_t> Edge::neighbours() const
    {
        return _neighbours;
    }

    std::ostream& operator<<(std::ostream &out, const Edge &edge)
    {
        return out << "(" << edge._a << ", " << edge._b << ")";
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

    void Edge::particle(Particle_Ptr new_particle)
    {
        _particle = new_particle;
    }

    Particle_Ptr Edge::particle() const
    {
        return _particle;
    }
}