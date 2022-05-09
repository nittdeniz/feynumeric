#ifndef FEYNUMERIC_ISOSPIN_HPP
#define FEYNUMERIC_ISOSPIN_HPP

#include <memory>

namespace Feynumeric{
	class Graph_Edge;
	double isospin(std::shared_ptr<Graph_Edge> out, std::shared_ptr<Graph_Edge> in, std::shared_ptr<Graph_Edge> pi);

	double isospin_tau(std::shared_ptr<Graph_Edge> a, std::shared_ptr<Graph_Edge> b);
	double isospin_T(std::shared_ptr<Graph_Edge> a, std::shared_ptr<Graph_Edge> b);
}

#endif