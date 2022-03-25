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
			if( list[permutation.at(i)].direction == Edge_Direction::IN ){
				result += std::pow(10, i);
			} else if( list[permutation.at(i)].direction == Edge_Direction::OUT ){
				result -= std::pow(10, i);
			} else{
				critical_error("directions_to_int has invalid direction.");
			}
		}
		return result;
	}

	int Vertex_Manager::directions_to_int(std::vector<Vertex::Particle_Direction> const& list){
		int result{0};
		for( std::size_t i = 0; i < list.size(); ++i )
		{
			if( list[i].direction == Edge_Direction::IN ){
				result += std::pow(10, i);
			} else if( list[i].direction == Edge_Direction::OUT ){
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
			if( vertex._particle_directions[i].direction == Edge_Direction::ANY )
			{
				Vertex a(vertex);
				Vertex b(vertex);
				a._particle_directions[i].direction = Edge_Direction::IN;
				b._particle_directions[i].direction = Edge_Direction::OUT;
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

	void Vertex_Manager::import(Vertex_Manager const& v)
	{
		for( auto const& child : v._root.children ){
			_root.children[child.first] = child.second;
		}
	}

	std::vector<std::vector<Vertex::Particle_Direction>>
	Vertex_Manager::advance(std::vector<Vertex::Particle_Direction> const& list) const
	{
		std::vector<std::vector<Vertex::Particle_Direction>> result;
		for( std::size_t i = 0; i < list.size(); ++i ){
			if( list[i].particle->parent() ){
				auto temp = list;
				temp[i].particle = list[i].particle->parent();
				std::sort(temp.begin(), temp.end());
				result.push_back(temp);
				//result.back()[i].particle = list[i].particle->parent();
			}
		}
		return result;
	}



	std::optional<Vertex_Manager::Vertex_Permutation> Vertex_Manager::find(std::vector<Vertex::Particle_Direction> const& particle_directions)
	{
		std::optional<Vertex_Manager::Vertex_Permutation> result;
		std::vector<std::vector<Vertex::Particle_Direction>> search_list = {particle_directions};
		do{
			std::vector<std::vector<Vertex::Particle_Direction>> updated_list;
			std::vector<std::size_t> found_indices;
			std::size_t depth{0};
			for( std::size_t i = 0; i < search_list.size(); ++i ){
				Node* current = &_root;
				auto& list = search_list[i];
				for( auto& entry : list ){
					if( current->children.contains(entry.particle->name()) ){
						current = &current->children[entry.particle->name()];
						depth++;
					}
					else{
						auto advanced_list = advance(list);
						updated_list.insert(updated_list.end(), advanced_list.begin(), advanced_list.end());
						break;
					}
				}
				if( depth == list.size() ){
					auto key = directions_to_int(list);
					if( current->data.contains(key) ){
						result = {current->data[key], current->permutation};
						found_indices.push_back(i);
					}
				}
			}
			if( found_indices.size() == 1 ){
				break;
			} else if( found_indices.size() > 1 ){
				critical_error("Found more than one suitable vertex.");
			}
			std::sort(updated_list.begin(), updated_list.end());
			auto last = std::unique(updated_list.begin(), updated_list.end(),
						   [](std::vector<Vertex::Particle_Direction> const& a, std::vector<Vertex::Particle_Direction> const& b)
						   {
								for( std::size_t i = 0; i < a.size(); ++i ){
									if( a[i].particle != b[i].particle ) return false;
								}
								return true;
						   }
						   );
			updated_list.erase(last, updated_list.end());
			search_list = updated_list;
		}while( !search_list.empty() );
		return result;
	}
}


