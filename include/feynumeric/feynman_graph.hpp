#ifndef FEYNUMERIC_FEYNMAN_GRAPH_HPP
#define FEYNUMERIC_FEYNMAN_GRAPH_HPP

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <vector>

#include "kinematics.hpp"
#include "matrix.hpp"
#include "four_vector.hpp"
#include "topology.hpp"
#include "types.hpp"
#include "vertex.hpp"


namespace Feynumeric{
	class Angular_Momentum;
	class Feynman_Diagram;
	class Particle;
	class Lorentz_Index;
	class Graph_Vertex;
	class Graph_Edge : public std::enable_shared_from_this<Graph_Edge>
	{
	private:
		std::size_t _eid;
		using Angular_Momentum_Ptr = std::shared_ptr<Angular_Momentum>;
		using Lorentz_Index_Ptr = std::shared_ptr<Lorentz_Index>;
		using Vertex_Ptr = std::shared_ptr<Graph_Vertex>;
		Feynman_Diagram* _diagram;
		Particle_Ptr _particle;
		Vertex_Ptr _front;
		Vertex_Ptr _back;
		Matrix _relative_momentum;
		Angular_Momentum_Ptr _spin;
		std::vector<Lorentz_Index_Ptr> _lorentz_indices;

		std::function<Matrix(Kinematics const&)> _feynman_rule;
	public:
		Graph_Edge(std::size_t id, Feynman_Diagram* diagram, Particle_Ptr const& P);
		Graph_Edge(Graph_Edge const& edge);
		Graph_Edge& operator=(Graph_Edge const& edge);

		bool is_incoming() const;
		bool is_outgoing() const;
		bool is_virtual() const;

		Vertex_Ptr front() const;
		Vertex_Ptr back() const;

		std::size_t id() const;

		Particle_Ptr particle() const;

		void spin(Angular_Momentum_Ptr const& spin);
		Angular_Momentum_Ptr spin() const;

		void add_lorentz_index(Lorentz_Index_Ptr const& index);
		void lorentz_indices(std::vector<Lorentz_Index_Ptr> const& list);
		std::vector<Lorentz_Index_Ptr> lorentz_indices() const;
		std::vector<Lorentz_Index_Ptr> lorentz_indices(Vertex_Ptr const& ptr) const;

		Four_Vector four_momentum(Kinematics const&) const;

		void front(Vertex_Ptr const& v);
		void back(Vertex_Ptr const& v);

		Matrix relative_momentum() const;
		void relative_momentum(Matrix const& momentum);

		std::function<Matrix(Kinematics const&)> feynman_rule();

		friend class Graph_Vertex;
	};

	using Edge_Ptr = std::shared_ptr<Graph_Edge>;

	class Graph_Vertex : public std::enable_shared_from_this<Graph_Vertex>
	{
	private:
		std::size_t _vid;
		Feynman_Diagram* _diagram;
		std::vector<Edge_Ptr> _front;
		std::vector<Edge_Ptr> _back;
		std::function<Matrix()> _feynman_rule;
		std::vector<Vertex::Particle_Direction> particle_directions();

	public:
		explicit Graph_Vertex(std::size_t id, Feynman_Diagram* diagram);
		Graph_Vertex(Graph_Vertex const& vertex);
		Graph_Vertex& operator=(Graph_Vertex const& vertex);

		std::vector<Edge_Ptr> front() const;
		std::vector<Edge_Ptr> back() const;
		std::vector<Edge_Ptr> all() const;

		std::size_t id() const;

		void front(Edge_Ptr const& e);
		void back(Edge_Ptr const& e);

		std::string particles_to_string() const;
		std::function<Matrix(Kinematics const& kin)> feynman_rule();
	};
	class Feynman_Graph{
	private:
		using Angular_Momentum_Ptr = std::shared_ptr<Angular_Momentum>;
		using Lorentz_Index_Ptr = std::shared_ptr<Lorentz_Index>;
		using Vertex_Ptr = std::shared_ptr<Graph_Vertex>;
		using Edge_Map = std::map<Vertex::ID, std::map<Vertex::ID, std::vector<Edge_Ptr>>>;
		using Vertex_Map = std::map<Vertex::ID, Vertex_Ptr>;

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