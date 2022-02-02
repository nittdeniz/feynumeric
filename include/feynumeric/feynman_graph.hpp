#ifndef FEYNUMERIC_FEYNMAN_GRAPH_HPP
#define FEYNUMERIC_FEYNMAN_GRAPH_HPP

#include <memory>
#include <vector>

namespace Feynumeric{
	class Feynman_Graph{
	private:
		class Edge;
		class Vertex;

		using Edge_Ptr = std::shared_ptr<Edge>;
		using Vertex_Ptr = std::shared_ptr<Vertex>;

		class Edge
		{
		private:
			std::shared_ptr<Vertex> front;
			std::shared_ptr<Vertex> back;
		};
		class Vertex
		{
		private:
			std::vector<Edge_Ptr> front;
			std::vector<Edge_Ptr> back;
		};
	public:
	};
}
#endif /// FEYNUMERIC_FEYNMAN_GRAPH_HPP