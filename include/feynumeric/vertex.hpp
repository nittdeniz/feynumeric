#ifndef Feynumeric_VERTEX_HPP
#define Feynumeric_VERTEX_HPP

#include <functional>
#include <memory>
#include <vector>

#include "edge_direction.hpp"
#include "kinematics.hpp"
#include "particle.hpp"
#include "particle_direction.hpp"

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
        friend bool operator<(Particle_Direction const& a, Particle_Direction const& b);
    };
    using Vertex_Ptr = std::shared_ptr<Vertex>;
}
#endif // Feynumeric_VERTEX_HPP