#ifndef FEYNUMERIC_FEYNMAN_GRAPH_HPP
#define FEYNUMERIC_FEYNMAN_GRAPH_HPP

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <vector>

#include "four_vector.hpp"
#include "kinematics.hpp"
#include "matrix.hpp"
#include "particle_direction.hpp"
#include "topology.hpp"
#include "types.hpp"
#include "vertex.hpp"


namespace Feynumeric{
	class Angular_Momentum;
	class Feynman_Diagram;
	class Particle;
	class Lorentz_Index;
	class Graph_Vertex;

	class Feynman_Graph{
	private:
        using Edge_Ptr = std::shared_ptr<Graph_Edge>;
		using Angular_Momentum_Ptr = std::shared_ptr<Angular_Momentum>;
		using Lorentz_Index_Ptr = std::shared_ptr<Lorentz_Index>;
		using Vertex_Ptr = std::shared_ptr<Graph_Vertex>;
		using Edge_Map = std::map<std::size_t, std::map<std::size_t, std::vector<Edge_Ptr>>>;
		using Vertex_Map = std::map<std::size_t, Vertex_Ptr>;

		Feynman_Diagram* _diagram;
		Topology _topology;

		std::vector<Edge_Ptr> _incoming;
		std::vector<Edge_Ptr> _outgoing;
		std::vector<Edge_Ptr> _virtual;
		std::vector<Edge_Ptr> _dummies;

		Edge_Map _edges;
		Vertex_Map _vertices;

		void create_graph(std::vector<Particle_Ptr> const& incoming_list, std::vector<Particle_Ptr> const& virtual_list, std::vector<Particle_Ptr> const& outgoing_list);
		void validate_input(std::vector<Particle_Ptr> const& incoming_list, std::vector<Particle_Ptr> const& virtual_list, std::vector<Particle_Ptr> const& outgoing_list);
	public:
		Feynman_Graph(Feynman_Diagram* diagram, Topology const& topology, std::vector<Particle_Ptr> const& incoming_list, std::vector<Particle_Ptr> const& virtual_list, std::vector<Particle_Ptr> const& outgoing_list);

		friend class Feynman_Diagram;
		friend class Feynman_Process;
		friend class Graph_Edge;
		friend class Graph_Vertex;
	};
}
#endif /// FEYNUMERIC_FEYNMAN_GRAPH_HPP