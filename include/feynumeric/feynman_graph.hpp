#ifndef FEYNUMERIC_FEYNMAN_GRAPH_HPP
#define FEYNUMERIC_FEYNMAN_GRAPH_HPP

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <vector>

#include "angular_momentum.hpp"
#include "kinematics.hpp"
#include "lorentz_index.hpp"
#include "matrix.hpp"
#include "four_vector.hpp"
#include "topology.hpp"
#include "types.hpp"


namespace Feynumeric{
	class Feynman_Diagram;
	class Feynman_Graph{
	public:
		class Edge;
		class Vertex;
	public:
		using Edge_Ptr = std::shared_ptr<Edge>;
		using Vertex_Ptr = std::shared_ptr<Vertex>;
	public:
	class Edge : public std::enable_shared_from_this<Feynman_Graph::Edge>
		{
		private:
			Feynman_Diagram* _diagram;
			Particle_Ptr _particle;
			Vertex_Ptr _front;
			Vertex_Ptr _back;
			Matrix _relative_momentum;
			Angular_Momentum_Ptr _spin;
			std::vector<Lorentz_Index_Ptr> _lorentz_indices;

			std::function<Matrix(Kinematics const&)> _feynman_rule;
		public:
			Edge(Feynman_Diagram* diagram, Particle_Ptr const& P);
			Edge(Edge const& edge);
			Edge& operator=(Edge const& edge);

			bool is_incoming() const;
			bool is_outgoing() const;
			bool is_virtual() const;

			Vertex_Ptr front() const;
			Vertex_Ptr back() const;

			Particle_Ptr particle() const;

			void spin(Angular_Momentum_Ptr const& spin);
			Angular_Momentum_Ptr spin() const;

			void add_lorentz_index(Lorentz_Index_Ptr const& index);
			std::vector<Lorentz_Index_Ptr> lorentz_indices() const;

			Four_Vector four_momentum(Kinematics const&) const;

			void front(Feynman_Graph::Vertex_Ptr const& v);
			void back(Feynman_Graph::Vertex_Ptr const& v);

			Matrix relative_momentum() const;
			void relative_momentum(Matrix const& momentum);

			std::function<Matrix(Kinematics const&)> feynman_rule();
		};

		class Vertex : public std::enable_shared_from_this<Feynman_Graph::Vertex>
		{
		private:
			Feynman_Diagram* _diagram;
			std::vector<Edge_Ptr> _front;
			std::vector<Edge_Ptr> _back;
			std::function<Matrix()> _feynman_rule;

		public:
			Vertex(Feynman_Diagram* diagram);
			Vertex(Vertex const& vertex);
			Vertex& operator=(Vertex const& vertex);

			std::vector<Edge_Ptr> front() const;
			std::vector<Edge_Ptr> back() const;
			std::vector<Edge_Ptr> all() const;

			void front(Edge_Ptr const& e);
			void back(Edge_Ptr const& e);
			std::size_t hash() const;


			std::function<Matrix(Kinematics const& kin)> feynman_rule();
		};
	private:
		enum class Vertex_Direction : unsigned char
		{
			IN = 0x01,
			OUT = 0x02
		};

		using Edge_Map = std::map<Topology::Vertex_ID, std::map<Topology::Vertex_ID, std::vector<Edge_Ptr>>>;
		using Vertex_Map = std::map<Topology::Vertex_ID, Vertex_Ptr>;

		Feynman_Diagram* _diagram;
		Topology _topology;

		std::vector<Edge_Ptr> _incoming;
		std::vector<Edge_Ptr> _outgoing;
		std::vector<Edge_Ptr> _virtual;

		Edge_Map _edges;
		Vertex_Map _vertices;

		void create_graph(std::vector<Particle_Ptr> const& incoming_list, std::vector<Particle_Ptr> const& virtual_list, std::vector<Particle_Ptr> const& outgoing_list);
		void validate_input(std::vector<Particle_Ptr> const& incoming_list, std::vector<Particle_Ptr> const& virtual_list, std::vector<Particle_Ptr> const& outgoing_list);
	public:
		Feynman_Graph(Feynman_Diagram* diagram, Topology const& topology, std::vector<Particle_Ptr> const& incoming_list, std::vector<Particle_Ptr> const& virtual_list, std::vector<Particle_Ptr> const& outgoing_list);

		friend class Feynman_Diagram;
		friend class Feynman_Process;
	};
}
#endif /// FEYNUMERIC_FEYNMAN_GRAPH_HPP