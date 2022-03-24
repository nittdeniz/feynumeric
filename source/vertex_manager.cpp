#include <algorithm>
#include <iostream>
#include <optional>
#include <sstream>

#include "vertex_manager.hpp"


namespace Feynumeric
{
	Vertex_Manager::Vertex_Manager()
	{
	}

	std::map<std::size_t, std::size_t>
	Vertex_Manager::permutation_map(std::vector<Vertex::Particle_Direction> const& lst){
		std::map<std::size_t, std::size_t> result;
		for( std::size_t i = 0; i < lst.size(); ++i ){
			std::size_t min_pos{0};
			while( result.contains(min_pos) ) min_pos++;
			for( std::size_t j = 0; j < lst.size(); ++j ){
				if( result.contains(j) ){
					continue;
				}
				if( lst[j] < lst[min_pos] ){
					min_pos = j;
				}
			}
			result[min_pos] = i;
		}
		return result;
	}

	std::map<std::size_t, std::size_t> Vertex_Manager::invert_map(std::map<std::size_t, std::size_t> const& map)
	{
		std::map<std::size_t, std::size_t> inverted;
		for( auto const& [key, value] : map )
		{
			inverted[value] = key;
		}
		return inverted;
	}

	int Vertex_Manager::directions_to_int(std::vector<Vertex::Particle_Direction> const& list, std::map<std::size_t, std::size_t> const& permutation){
		int result{0};
		for( std::size_t i = 0; i < list.size(); ++i )
		{
			if( list[permutation.at(i)].direction == Direction::INCOMING ){
				result += std::pow(10, i);
			} else if( list[permutation.at(i)].direction == Direction::OUTGOING ){
				result -= std::pow(10, i);
			} else{
				critical_error("directions_to_int has invalid direction.");
			}
		}
		return result;
	}

	void Vertex_Manager::add(Vertex const& vertex){
		for( std::size_t i = 0; i < vertex._particle_directions.size(); ++i )
		{
			if( vertex._particle_directions[i].direction == Direction::ANY )
			{
				Vertex a(vertex);
				Vertex b(vertex);
				a._particle_directions[i].direction = Direction::INCOMING;
				b._particle_directions[i].direction = Direction::OUTGOING;
				add(a);
				add(b);
				return;
			}
		}
		auto permutation = permutation_map(vertex._particle_directions);
		Node* current = &_root;
		for( std::size_t i = 0; i < vertex._particle_directions.size(); ++i ){
			auto pd = vertex._particle_directions[permutation[i]];
			current = &(current->children[pd.particle->name()]);
		}

		auto key = directions_to_int(vertex._particle_directions, permutation);
		current->data.insert({key, vertex});
		current->permutation = invert_map(permutation);
	}

	std::optional<Vertex> Vertex_Manager::find(std::vector<Vertex::Particle_Direction> lst)
	{
		std::optional<Vertex> result;
		do{
			Node* current = &_root;
			for( auto const& item : lst )
			{
				if( current->children.contains(item.particle->name()) )
				{
					current = &current->children[item.particle->name()];
				}
				else
				{
					// advance particles to parents
				}
			}
		}while(true);


		std::sort(lst.begin(), lst.end());

		return result;
	}
}


