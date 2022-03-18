#ifndef FEYNUMERIC_FEYNMAN_RULE_HPP
#define FEYNUMERIC_FEYNMAN_RULE_HPP

#include <iostream>
#include <string>

#include "feynman_graph.hpp"
#include "kinematics.hpp"

namespace Feynumeric
{
	class Feynman_Rule
	{
	public:
		Feynman_Rule(std::string&& str);
		Feynman_Rule(std::string const& str);
		Feynman_Rule(Feynman_Rule const& other);
		Feynman_Rule(Feynman_Rule&& other);
		Feynman_Rule& operator=(Feynman_Rule const& other);
		Feynman_Rule& operator=(Feynman_Rule&& other);

		virtual Matrix operator()() = 0;
		virtual void print(Feynman_Diagram* diagram);
	};

	class Feynman_Edge_Rule : public Feynman_Rule
	{
	private:
		Feynman_Graph::Edge_Ptr _edge;
		std::function<Matrix(Feynman_Graph::Edge_Ptr const&, Kinematics const&)> _rule;
	public:
		virtual Matrix operator()(Kinematics const& kin)
		{
			return _rule(_edge, kin);
		}
		virtual void print(Feynman_Diagram*)
		{
			std::cout << "feynman_edge\n";
		}
	};

	class Feynman_Vertex_Rule : public Feynman_Rule
	{
	private:
		Feynman_Graph::Edge_Ptr _edge;
		std::function<Matrix(Feynman_Graph::Edge_Ptr const&, Kinematics const&)> _rule;
	public:
		virtual Matrix operator()(Kinematics const& kin)
		{
			return _rule(_edge, kin);
		}
		virtual void print(Feynman_Diagram*)
		{
			std::cout << "feynman_vertex\n";
		}
	};
}
#endif