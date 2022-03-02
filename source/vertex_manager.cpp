#include <iostream>
#include <optional>
#include <sstream>

#include "vertex_manager.hpp"


namespace Feynumeric
{
	Vertex_Manager::Vertex_Manager(std::shared_ptr<Vertex_Manager> const& copy)
	: _vertices(copy->_vertices)
	{

	}

	void Vertex_Manager::add(Vertex const& vertex)
	{
		for( std::size_t i = 0; i < vertex._particle_directions.size(); ++i )
		{
			if( vertex._particle_directions[i].direction == Direction::BOTH )
			{
				Vertex v1(vertex);
				v1._particle_directions[i].direction = Direction::INCOMING;
				Vertex v2(vertex);
				v2._particle_directions[i].direction = Direction::OUTGOING;
				add(v1);
				add(v2);
				return;
			}
		}
		auto hash = canonical_hash(vertex._particle_directions);
		_vertices[hash] = std::make_shared<Vertex>(vertex);
	}

	std::optional<Vertex_Ptr> Vertex_Manager::find_vertex(std::size_t hash)
	{
		if( _vertices.contains(hash) )
		{
			return _vertices.at(hash);
		}
		return std::nullopt;
	}
}


