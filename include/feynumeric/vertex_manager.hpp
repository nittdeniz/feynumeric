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
    	std::map<std::size_t, Vertex_Ptr> _vertices;
    public:
    	Vertex_Manager() = default;
    	Vertex_Manager(std::shared_ptr<Vertex_Manager> const& copy);
    	void add(Vertex const& vertex);
    	std::optional<Vertex_Ptr> find_vertex(std::size_t hash);
    };
    using Vertex_Manager_Ptr = std::shared_ptr<Vertex_Manager>;
}
#endif // Feynumeric_VERTEX_MANAGER_HPP