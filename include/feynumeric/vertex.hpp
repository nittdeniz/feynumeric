#ifndef Feynumeric_VERTEX_HPP
#define Feynumeric_VERTEX_HPP

#include <functional>
#include <memory>
#include <vector>

#include "direction.hpp"
#include "feynman_graph.hpp"
#include "kinematics.hpp"
#include "particle.hpp"

namespace Feynumeric
{
	using Particle_List = std::vector<Feynman_Graph::Edge_Ptr>;
    class Vertex
    {
    public:
    	struct Particle_Direction
	    {
    		Particle_Ptr particle;
    		Direction direction;
    		Particle_Direction(Particle_Ptr ptr, Direction dir);
    		Particle_Direction(Particle_Direction const& pd);
    		Particle_Direction& operator=(Particle_Direction const& pd);
	    };
    private:
	    std::vector<Particle_Direction> _particle_directions;
	    std::function<Matrix(Kinematics const&, Particle_List const&)> _vertex_function;
    public:
		Vertex(std::vector<Particle_Direction> const& particle_directions, std::function<Matrix(Kinematics const&, Particle_List const&)>&&);
        Vertex(Vertex const& v);
        Vertex& operator=(Vertex const& v);

        Particle_List sort(Particle_List list, Feynman_Graph::Vertex_Ptr vertex);

        std::function<Matrix(Kinematics const&, Particle_List const&)> vertex_function() const;


        friend class Vertex_Manager;
    };
	bool canonical_sort(Vertex::Particle_Direction const& a, Vertex::Particle_Direction const& b);
	std::size_t canonical_hash(std::vector<Vertex::Particle_Direction> copy);
    using Vertex_Ptr = std::shared_ptr<Vertex>;
}
#endif // Feynumeric_VERTEX_HPP