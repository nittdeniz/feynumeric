#include "feynman_graph.hpp"
#include "isospin.hpp"
#include "matrix.hpp"
#include "messages.hpp"

namespace Feynumeric{
	double
	isospin(std::shared_ptr<Graph_Edge> out, std::shared_ptr<Graph_Edge> in, std::shared_ptr<Graph_Edge> pi){
		static const Matrix T_12_12_pi_in(2, 2, {
			1, -std::sqrt(2),
			std::sqrt(2), -1
		});
		static const Matrix T_12_12_pi_out(T_12_12_pi_in.T());
		static const Matrix T_32_12_pi_in(4, 2, {
			1, 0,
			std::sqrt(2/3.), std::sqrt(1/3.),
			std::sqrt(1/3.), std::sqrt(2/3.),
			0, 1
		});
		static const Matrix T_32_12_pi_out(4, 2, {
			1, 0,
			-std::sqrt(2/3.), std::sqrt(1/3.),
			std::sqrt(1/3.), -std::sqrt(2/3.),
			0, 1
		});
		static const Matrix T_12_32_pi_in(T_32_12_pi_in.T());
		static const Matrix T_12_32_pi_out(T_32_12_pi_out.T());
		static const Matrix T_32_32_pi_in(4, 4, {
			1, -std::sqrt(2/3.), 0, 0,
			std::sqrt(2/3.), 1/3., -std::sqrt(8/3.), 0,
			0, std::sqrt(8/3.), -1/3., -std::sqrt(2/3.),
			0, 0, std::sqrt(2/3.), -1
		});
		static const Matrix T_32_32_pi_out(T_32_32_pi_in.T());
		auto const iso_out = out->particle()->isospin();
		auto const iso_in  = in->particle()->isospin();
		bool const incoming_pion = pi->front() == in->front();

		auto const index_out = static_cast<std::size_t>(iso_out.j() - iso_out.m());
		auto const index_in  = static_cast<std::size_t>(iso_in.j() - iso_in.m());

		if( iso_out.j() == 1.5 && iso_in.j() == 1.5 ){
			if( incoming_pion ){
				return T_32_32_pi_in.at(index_out, index_in).real();
			}else{
				return T_32_32_pi_out.at(index_out, index_in).real();
			}
		}
		if( iso_out.j() == 1.5 && iso_in.j() == 0.5 ){
			if( incoming_pion ){
				return T_32_12_pi_in.at(index_out, index_in).real();
			}else{
				return T_32_12_pi_out.at(index_out, index_in).real();
			}
		}
		if( iso_out.j() == 0.5 && iso_in.j() == 1.5 ){
			if( incoming_pion ){
				return T_12_32_pi_in.at(index_out, index_in).real();
			}else{
				return T_12_32_pi_out.at(index_out, index_in).real();
			}
		}
		if( iso_out.j() == 0.5 && iso_in.j() == 0.5 ){
			if( incoming_pion ){
				return T_12_12_pi_in.at(index_out, index_in).real();
			}else{
				return T_12_12_pi_out.at(index_out, index_in).real();
			}
		}
		critical_error(FORMAT("Unsupported isospin matrix for particles {} {} {}", out->particle()->name(), in->particle()->name(), pi->particle()->name()));
	}

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

