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

	std::function<Matrix(Kinematics const&, Particle_List const&)> Vertex::vertex_function() const
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

	Vertex::Particle_Direction::Particle_Direction(Particle_Ptr ptr, Direction dir)
	: particle(ptr)
	, direction(dir)
	{

	}

	Vertex::Particle_Direction::Particle_Direction(Vertex::Particle_Direction const& pd)
	: particle(pd.particle)
	, direction(pd.direction)
	{

	}

	Vertex::Particle_Direction& Vertex::Particle_Direction::operator=(Vertex::Particle_Direction const& pd)
	{
		particle = pd.particle;
		direction = pd.direction;
		return *this;
	}

	bool operator<(Vertex::Particle_Direction const& a, Vertex::Particle_Direction const& b)
	{
		if( a.particle->name() == b.particle->name() )
		{
			return a.direction < b.direction;
		}
		return a.particle->name() < b.particle->name();
	}
}