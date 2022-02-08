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

	bool canonical_sort(Vertex::Particle_Direction const& a, Vertex::Particle_Direction const& b)
	{
		if( a.particle->name() == b.particle->name() )
		{
			return a.direction < b.direction;
		}
		return a.particle->name() < b.particle->name();
	}

	std::size_t canonical_hash(std::vector<Vertex::Particle_Direction> copy)
	{
		std::sort(copy.begin(), copy.end(), canonical_sort);
		std::stringstream stream;
		for( auto const& item : copy )
		{
			stream << item.particle->name() << item.direction;
		}
		auto str = stream.str();
		return std::hash<std::string>{}(str);
	}

	Particle_List Vertex::sort(Particle_List list, Feynman_Graph::Vertex_Ptr vertex)
	{
		Particle_List result;
		result.reserve(list.size());
		for( auto const& item : _particle_directions )
		{
			bool found = false;
			for( auto it = list.begin(); it != list.end(); it++ ){
				auto direction = (*it)->front() == vertex ? Direction::INCOMING : Direction::OUTGOING;
				if( item.particle == (*it)->particle() && item.direction == direction )
				{
					found = true;
					result.push_back(*it);
					list.erase(it);
					break;
				}
			}
			if( !found )
			{
				critical_error("Could not find corresponding edge at vertex.");
			}
		}
		return result;
	}
}