#ifndef Feynumeric_VERTEX_MANAGER_HPP
#define Feynumeric_VERTEX_MANAGER_HPP

#include <functional>
#include <map>
#include <memory>
#include <utility>
#include <string>


#include "vertex.hpp"

namespace Feynumeric
{
    class Vertex_Manager
    {
    private:
		struct Node
		{
			std::map<int, Vertex> data;
			std::map<std::size_t, std::size_t> permutation;
			std::map<std::string, Node> children;
		};

		int directions_to_int(std::vector<Vertex::Particle_Direction> const& list, std::map<std::size_t, std::size_t> const& permutation);

		Node _root;
		std::map<std::size_t, std::size_t> permutation_map(std::vector<Vertex::Particle_Direction> const& lst);
		std::map<std::size_t, std::size_t> invert_map(std::map<std::size_t, std::size_t> const& map);

    public:
    	Vertex_Manager();
    	void add(Vertex const& vertex);

    	std::optional<Vertex> find(std::vector<Vertex::Particle_Direction> lst);
    };
    using Vertex_Manager_Ptr = std::shared_ptr<Vertex_Manager>;
}
#endif // Feynumeric_VERTEX_MANAGER_HPP