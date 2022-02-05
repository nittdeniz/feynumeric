#ifndef FEYNUMERIC_FEYNMAN_GRAPH_HPP
#define FEYNUMERIC_FEYNMAN_GRAPH_HPP
/*
#include <functional>
#include <list>
#include <memory>
#include <vector>


namespace Feynumeric{
	class Feynman_Graph{
		class Edge;
		class Vertex;

		using Edge_Ptr = std::shared_ptr<Edge>;
		using Vertex_Ptr = std::shared_ptr<Vertex>;

		class Edge
		{
		private:
			Feynman_Diagram* _diagram;
			Particle const& P;
			Vertex_Ptr _front;
			Vertex_Ptr _back;

			std::function<Matrix()> _feynman_rule;
		public:
			Edge(Feynman_Diagram* diagram, Particle const& P);
			Edge(Edge const& edge);
			Edge& operator=(Edge const& edge);

			Vertex_Ptr front() const;
			Vertex_Ptr back() const;

			Particle const& particle();

			void front(Vertex_Ptr const& v);
			void back(Vertex_Ptr const& v);

			std::function<Matrix()> feynman_rule() const;
		};

		class Vertex
		{
		private:
			Feynman_Diagram* _diagram;
			std::vector<Edge_Ptr> _front;
			std::vector<Edge_Ptr> _back;
			std::function<Matrix()> _feynman_rule;

		public:
			Vertex(Feynman_Diagram* diagram);
			Vertex(Vertex const& Vertex);
			Vertex& operator=(Vertex const& Vertex);

			std::vector<Edge_Ptr> front() const;
			std::vector<Edge_Ptr> back() const;

			void front(Edge_Ptr const& e);
			void back(Edge_Ptr const& e);

			std::function<Matrix()> feynman_rule() const;


		};
	private:
		std::list<Edge> _incoming_edges;
		std::list<Edge> _outgoing_edges;
		std::list<Vertex> _vertices;

		void create_graph();
	public:
		Feynman_Graph(Feynman_Diagram* diagram, Topology const& topology);
	};
}
 */
#endif /// FEYNUMERIC_FEYNMAN_GRAPH_HPP