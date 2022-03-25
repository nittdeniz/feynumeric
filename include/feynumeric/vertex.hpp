#ifndef Feynumeric_VERTEX_HPP
#define Feynumeric_VERTEX_HPP

#include <functional>
#include <memory>
#include <vector>

#include "edge_direction.hpp"
#include "kinematics.hpp"
#include "particle.hpp"

namespace Feynumeric
{
	class Graph_Edge;
	class Graph_Vertex;
    class Vertex
    {
    public:
	    using Edge_Ptr = std::shared_ptr<Graph_Edge>;
	    using Vertex_Ptr = std::shared_ptr<Graph_Vertex>;
	    using Particle_List = std::vector<Edge_Ptr>;

    	struct Particle_Direction
	    {
    		Particle_Ptr particle;
    		Edge_Direction direction;
    		Particle_Direction(Particle_Ptr ptr, Edge_Direction dir);
    		Particle_Direction(Particle_Direction const& pd);
    		Particle_Direction& operator=(Particle_Direction const& pd);
	    };
    private:
	    std::vector<Particle_Direction> _particle_directions;
	    std::function<Matrix(Kinematics const&, Particle_List const&)> _vertex_function;
    public:
    	Vertex() = default;
		Vertex(std::vector<Particle_Direction> const& particle_directions, std::function<Matrix(Kinematics const&, Particle_List const&)>&&);
        Vertex(Vertex const& v);
        Vertex& operator=(Vertex const& v);

        Particle_List sort(Particle_List list, Vertex_Ptr vertex);
        Particle_List sort(Particle_List const& list, std::map<std::size_t, std::size_t> const& permutation);

        std::function<Matrix(Kinematics const&, Particle_List const&)> vertex_function() const;


        friend class Graph_Vertex;
        friend class Vertex_Manager;
        friend bool operator<(Vertex::Particle_Direction const& a, Vertex::Particle_Direction const& b);
    };
	bool operator<(Vertex::Particle_Direction const& a, Vertex::Particle_Direction const& b);
    using Vertex_Ptr = std::shared_ptr<Vertex>;
}
#endif // Feynumeric_VERTEX_HPP