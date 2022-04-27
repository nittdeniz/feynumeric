#include "feynman_graph.hpp"
#include "isospin.hpp"
#include "matrix.hpp"
#include "messages.hpp"

namespace Feynumeric{
	double isospin_tau(std::shared_ptr<Graph_Edge> a, std::shared_ptr<Graph_Edge> b){
		static const Matrix tau(2, 2, {1, std::sqrt(2.), std::sqrt(2.), -1});
		if( a->particle()->isospin().j() != 0.5 || b->particle()->isospin().j() != 0.5 )
		{
			Feynumeric::critical_error(FORMAT("Both {} and {} must have isospin_j = 0.5\n", a->particle()->name(), b->particle()->name()));
		}
		std::size_t in, out;
		if( a->front() == b->back() )
		{
			in = 1-(a->particle()->isospin().m() + 0.5);
			out = 1-(b->particle()->isospin().m() + 0.5);
		}
		else if( a->back() == b->front() )
		{
			in = 1-(b->particle()->isospin().m() + 0.5);
			out = 1-(a->particle()->isospin().m() + 0.5);
		}
		else
		{
			Feynumeric::critical_error("Unsupported isospin for meeting fermions.");
		}
		return tau.at(out, in).real();
	}

	double isospin_T(std::shared_ptr <Graph_Edge> isospin32, std::shared_ptr <Graph_Edge> isospin12){
		static const Matrix T(2, 4, {
				-1, std::sqrt(2./3.), std::sqrt(1./3.), 0,
				0, -std::sqrt(1./3.), std::sqrt(2./3.), 1
		});
		std::size_t indexDelta = static_cast<std::size_t>(1.5 - isospin32->particle()->isospin().m());
		std::size_t indexN = static_cast<std::size_t>(0.5 - isospin12->particle()->isospin().m());
		return T.at(indexN, indexDelta).real();
	}
}

