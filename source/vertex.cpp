#include <algorithm>

#include "vertex.hpp"

namespace Feynumeric
{
	Vertex::Vertex(std::vector<Particle_Direction> const& particle_directions,
                                     std::function<Matrix(Kinematics const&, Particle_List const&)>&& f)
    : _particle_directions(particle_directions)
    , _vertex_function(f)
	{
	}

	std::function<Matrix(Kinematics const&, Vertex::Particle_List const&)> Vertex::vertex_function() const
	{
		return _vertex_function;
	}

	Vertex::Vertex(Vertex const& v)
	: _particle_directions(v._particle_directions)
	, _vertex_function(v._vertex_function)
	{

	}

	Vertex& Vertex::operator=(Vertex const& v)
    {
        _particle_directions = v._particle_directions;
        _vertex_function = v._vertex_function;
        return *this;
    }
}