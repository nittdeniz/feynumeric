/*
#include <iostream>
#include "diagram.hpp"
#include "edge.hpp"
#include "messages.hpp"
#include "particle.hpp"

namespace Feynumeric
{
    Edge::Edge(std::size_t a, std::size_t b, Edge::Type type, Particle_Ptr particle)
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
    , _neighbours(edge._neighbours)
    , _assigned_indices(edge._assigned_indices)
    , _angular_momentum(edge._angular_momentum)
    {
		std::cerr << "Edge(Edge const& edge)\n";
    }

    Edge &Edge::operator=(const Edge &edge)
    {
    	std::cerr << "operator=(Edge const& edge)\n";
        _a = edge._a;
        _b = edge._b;
        _type = edge._type;
        _particle = edge._particle;
        _momentum = edge._momentum;
        _neighbours = edge._neighbours;
        _assigned_indices = edge._assigned_indices;
        _angular_momentum = edge._angular_momentum;
        return *this;
    }



    void Edge::set_diagram(Diagram* diagram)
    {
    	_diagram = diagram;
    }

    std::size_t Edge::a() const
    {
        return _a;
    }

    std::size_t Edge::b() const
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

    void Edge::add_neighbour(Edge* neighbour)
    {
        _neighbours.push_back(neighbour);
    }

	void Edge::clear_neighbours(){
		_neighbours.clear();
	}

	std::vector<Edge*> Edge::neighbours()
    {
        return _neighbours;
    }

    Matrix Edge::momentum() const
    {
        return _momentum;
    }

    std::ostream& operator<<(std::ostream &out, const Edge &edge)
    {
        return out << "(" << static_cast<std::size_t>(edge._a) << ", " << static_cast<std::size_t>(edge._b) << ")";
    }

    std::vector<Lorentz_Index_Ptr> Edge::get_lorentz_indices() const
    {
    	return _assigned_indices;
    }

    Momentum_Func Edge::four_momentum() const
    {
    	if( _diagram == nullptr )
	    {
    		critical_error("Edge::four_momentum() _diagram == nullptr.");
	    }
    	return _diagram->four_momentum(_momentum);
    }

    Angular_Momentum_Ptr Edge::spin() const
    {
    	return _angular_momentum;
    }

    std::function<Matrix(Kinematics const&)> Edge::feynman_rule()
    {
	    using namespace std::placeholders;
        if( is_incoming() )
        {
            return std::bind(_particle->feynman_incoming, this, _1);
        }
        else if( is_outgoing() )
        {
            return std::bind(_particle->feynman_outgoing, this, _1);
        }
        else if( is_virtual() )
        {
            return std::bind(_particle->feynman_virtual, this, _1);
        }
        else
        {
            critical_error("Can't assign feynman rule to edge: " + to_string() + ".");
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

    bool shares_vertex(Edge* lhs, Edge* rhs)
    {
        return lhs->_a == rhs->_a || lhs->_a == rhs->_b
            || lhs->_b == rhs->_a || lhs->_b == rhs->_b;
    }

    std::string Edge::to_string() const
    {
        std::stringstream in;
        in << *this;
        return in.str();
    }

    void Edge::assign_lorentz_index(Lorentz_Index_Ptr const& id)
    {
        _assigned_indices.emplace_back(id);
    }

    void Edge::clear_lorentz_indices()
    {
        _assigned_indices.clear();
    }


    void Edge::assign_angular_momentum(Angular_Momentum_Ptr const& id)
    {
	    _angular_momentum = id;
    }

    std::optional<std::size_t> shared_vertex(Edge* lhs, Edge* rhs)
    {
        if( lhs->_a == rhs->_a || lhs->_a == rhs->_b )
        {
            return lhs->_a;
        }
        if( lhs->_b == rhs->_a || lhs->_b == rhs->_b )
        {
            return lhs->_b;
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
*/