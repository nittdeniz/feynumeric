#include <algorithm>
#include <iostream>
#include <utility.hpp>

#include "feynman_diagram.hpp"
#include "format.hpp"
#include "messages.hpp"
#include "particle.hpp"

namespace Feynumeric
{
	Feynman_Diagram::Feynman_Diagram(Topology const& topology, Vertex_Manager_Ptr const& VMP, std::vector<Particle_Ptr> const& incoming_particles,
	                                 std::vector<Particle_Ptr> const& virtual_particles,
	                                 std::vector<Particle_Ptr> const& outgoing_particles)
			: _graph(this, topology, incoming_particles, virtual_particles, outgoing_particles)
			, _VMP(VMP)
	{
	}

	void Feynman_Diagram::add_spin(Feynman_Graph::Edge_Ptr const& edge_ptr)
	{
		auto const& spin = edge_ptr->particle()->spin();
		if( spin.j() > 0 )
		{
			_spins.push_back(std::make_shared<Angular_Momentum>(spin.j(), spin.j()));
			edge_ptr->spin(_spins.back());
		}
	}

	void Feynman_Diagram::add_lorentz_indices(Feynman_Graph::Edge_Ptr const& edge_ptr)
	{
		double j = edge_ptr->particle()->spin().j();
		while( j --> 0.5 )
		{
			_lorentz_indices.push_back(std::make_shared<Lorentz_Index>());
			edge_ptr->add_lorentz_index(_lorentz_indices.back());
		}
	}

	void Feynman_Diagram::generate_amplitude()
	{
		fix_momenta();
		_spins.reserve(_graph._outgoing.size() + _graph._incoming.size());
		for( auto& edge_ptr : _graph._outgoing )
		{
			add_spin(edge_ptr);
			add_lorentz_indices(edge_ptr);

			if( edge_ptr->particle()->is_fermion() )
			{
//				std::cerr << "Amplitude-Feynman: " << edge_ptr->particle()->name() << "\n";
//				for( auto const& ind : edge_ptr->lorentz_indices() )
//				{
//					std::cerr << ind << "\n";
//				}
				_amplitude.push_back(edge_ptr->feynman_rule());
				trace_fermion_line(edge_ptr, edge_ptr->back(),Direction::OUTGOING);
			}
			else if( !edge_ptr->particle()->is_anti_fermion() )
			{
//				std::cerr << "Amplitude-Feynman: " << edge_ptr->particle()->name() << "\n";
//				for( auto const& ind : edge_ptr->lorentz_indices() )
//				{
//					std::cerr << ind << "\n";
//				}
				_amplitude.push_back(edge_ptr->feynman_rule());
			}
		}
		for( auto& edge_ptr : _graph._incoming )
		{
			add_spin(edge_ptr);
			add_lorentz_indices(edge_ptr);
			if( edge_ptr->particle()->is_anti_fermion() )
			{
//				std::cerr << "Amplitude-Feynman: " << edge_ptr->particle()->name() << "\n";
//				for( auto const& ind : edge_ptr->lorentz_indices() )
//				{
//					std::cerr << ind << "\n";
//				}
				_amplitude.push_back(edge_ptr->feynman_rule());
				trace_fermion_line(edge_ptr, edge_ptr->front(), Direction::INCOMING);
			}
			else if( !edge_ptr->particle()->is_fermion() )
			{
//				std::cerr << "Amplitude-Feynman: " << edge_ptr->particle()->name() << "\n";
//				for( auto const& ind : edge_ptr->lorentz_indices() )
//				{
//					std::cerr << ind << "\n";
//				}
				_amplitude.push_back(edge_ptr->feynman_rule());
			}
		}
		for( auto& edge_ptr : _graph._virtual )
		{
			// do this twice since it's a virtual particle
			add_lorentz_indices(edge_ptr);
			add_lorentz_indices(edge_ptr);
			if( !(edge_ptr->particle()->is_fermion() || edge_ptr->particle()->is_anti_fermion()) )
			{
//				std::cerr << "Amplitude-Feynman: " << edge_ptr->particle()->name() << "\n";
//				for( auto const& ind : edge_ptr->lorentz_indices() )
//				{
//					std::cerr << ind << "\n";
//				}
				_amplitude.push_back(edge_ptr->feynman_rule());
			}
		}
	}

	void Feynman_Diagram::fix_momenta()
	{
		std::size_t const n_external = _graph._topology._incoming_edges.size() + _graph._topology._outgoing_edges.size();
		Matrix zeroes(n_external, 1);
		std::size_t position{0};

		for( auto const& edge_ptr : _graph._incoming )
		{
			Matrix p = zeroes;
			p[position] = 1;
			edge_ptr->relative_momentum(p);
//			_four_momenta.push_back(generate_four_momentum(Direction::INCOMING, position));
			position++;
		}
		for( auto const& edge_ptr : _graph._outgoing )
		{
			Matrix p = zeroes;
			p[position] = -1;
			edge_ptr->relative_momentum(p);
//			_four_momenta.push_back(generate_four_momentum(Direction::OUTGOING, position));
			position++;
		}

		auto virtuals_left = _graph._virtual;

		auto n_undetermined_momenta = [&](Feynman_Graph::Vertex_Ptr const& vertex)
		{
			std::size_t n{0};
			std::size_t n_edges{0};
			for( auto const& e : vertex->all() )
			{
				n_edges++;
				if( e->relative_momentum().n_cols() == 1 )
				{
					n++;
				}
			}
			return n_edges - n;
		};

		while( !virtuals_left.empty() )
		{
			bool changed = false;
			for( auto const& edge_ptr : virtuals_left )
			{
				if( edge_ptr->relative_momentum().n_cols() == 1 )
				{
					virtuals_left.erase(std::remove(virtuals_left.begin(), virtuals_left.end(), edge_ptr), virtuals_left.end());
					changed = true;
					break;
				}
				if( n_undetermined_momenta(edge_ptr->back()) == 1 )
				{
					Matrix sum_outgoing = zeroes;
					Matrix sum_incoming = zeroes;
					for( auto const& incoming_ptr : edge_ptr->back()->back() )
					{
						sum_incoming += incoming_ptr->relative_momentum();
					}
					for( auto const& outgoing_ptr : edge_ptr->back()->front() )
					{
						if( outgoing_ptr == edge_ptr )
							continue;
						sum_outgoing += outgoing_ptr->relative_momentum();
					}
					edge_ptr->relative_momentum(sum_outgoing + sum_incoming);
					changed = true;
					virtuals_left.erase(std::remove(virtuals_left.begin(), virtuals_left.end(), edge_ptr), virtuals_left.end());
					break;
				}
				if( n_undetermined_momenta(edge_ptr->front()) == 1 )
				{
					Matrix sum_outgoing = zeroes;
					Matrix sum_incoming = zeroes;
					for( auto const& incoming_ptr : edge_ptr->back()->back() )
					{
						if( incoming_ptr == edge_ptr )
							continue;
						sum_incoming += incoming_ptr->relative_momentum();
					}
					for( auto const& outgoing_ptr : edge_ptr->front()->front() )
					{
						sum_outgoing += outgoing_ptr->relative_momentum();
					}
					changed = true;
					virtuals_left.erase(std::remove(virtuals_left.begin(), virtuals_left.end(), edge_ptr), virtuals_left.end());
					break;
				}
			}
			if( !changed )
			{
				critical_error("Can not fix momenta.");
			}
		}
		for( auto const& [vertex_a, map] : _graph._edges )
		{
			for( auto const& [vertex_b, edge_list] : map )
			{
				if( vertex_a > vertex_b ) continue;
//				for( auto const& edge_ptr : edge_list )
//				{
//					std::cerr << "Edge(" << vertex_a << ", " << vertex_b << "): " << edge_ptr->relative_momentum() << "\n";
//				}
			}
		}

	}

/*
	Feynman_Diagram::Momentum_Func Feynman_Diagram::generate_four_momentum(Direction const& direction, std::size_t pos) const
	{
		if( direction == Direction::INCOMING )
		{
			std::size_t n_particles = _graph._incoming.size();
			switch( n_particles )
			{
				case 1:
					return [](Particle_Ptr const& P, Kinematics const& kin){
						return Four_Momentum(kin.momenta[0], P->mass(), 1);
					};
				case 2:
					switch( pos )
					{
						case 0:
							return [](Particle_Ptr const& P, Kinematics const& kin){
								return Four_Momentum(kin.momenta[0], P->mass(), 1);
							};
						case 1:
							return [](Particle_Ptr const& P, Kinematics const& kin){
								return Four_Momentum(-kin.momenta[0], P->mass(), 1);
							};
						default:
							critical_error(FORMAT("generate_four_momentum: position argument ({}) out of bounds ({}).", pos, n_particles));
					};
				default:
					break;
			}
		}
		if( direction == Direction::OUTGOING )
		{
			std::size_t n_particles = _graph._outgoing.size();
			std::size_t offset = _graph._incoming.size();
			switch( n_particles )
			{
				case 1:
					return [](Particle_Ptr const& P, Kinematics const& kin){
						return Four_Momentum(kin.momenta[1], P->mass(), 1);
					};
				case 2:
					switch( pos - offset )
					{
						case 0:
							return [offset](Particle_Ptr const& P, Kinematics const& kin){
								return Four_Momentum(kin.momenta[1], P->mass(), kin.cosines[0]);
							};
						case 1:
							return [offset](Particle_Ptr const& P, Kinematics const& kin){
								return Four_Momentum(-kin.momenta[1], P->mass(), kin.cosines[0]);
							};
						default:
							critical_error(FORMAT("generate_four_momentum: position argument ({}) out of bounds ({}).", pos, n_particles));
					};
				default:
					break;
			}
		}
		critical_error("Particle configuration not supported.");
	}
*/
//	Four_Momentum Feynman_Diagram::four_momentum(std::size_t index, Particle_Ptr const& P, Kinematics const& kin)
//	{
//		return _four_momenta.at(index)(P, kin);
//	}

	void Feynman_Diagram::trace_fermion_line(Feynman_Graph::Edge_Ptr const& ptr, Direction const& start_direction)
	{
		switch( start_direction )
		{
			case Direction::OUTGOING:
				if( ptr->is_incoming() )
				{
					if( ptr->particle()->is_fermion() )
					{
//						std::cerr << "Amplitude-Feynman: " << ptr->particle()->name() << "\n";
//						for( auto const& ind : ptr->lorentz_indices() )
//						{
//							std::cerr << ind << "\n";
//						}
						_amplitude.push_back(ptr->feynman_rule());
						return;
					}
					critical_error(FORMAT("Outgoing fermion line ends in incoming line, which is not a fermion ({}).", ptr->particle()->name()));
				}
				if( ptr->is_outgoing() )
				{
					if( ptr->particle()->is_anti_fermion() )
					{
//						std::cerr << "Amplitude-Feynman: " << ptr->particle()->name() << "\n";
//						for( auto const& ind : ptr->lorentz_indices() )
//						{
//							std::cerr << ind << "\n";
//						}
						_amplitude.push_back(ptr->feynman_rule());
						return;
					}
					critical_error(FORMAT("Outgoing fermion line ends in outgoing line, which is not an antifermion ({}).", ptr->particle()->name()));
				}
				trace_fermion_line(ptr, ptr->back(), start_direction);
				break;
			case Direction::INCOMING:
				if( ptr->is_outgoing() )
				{
					if( ptr->particle()->is_anti_fermion() )
					{
//						std::cerr << "Amplitude-Feynman: " << ptr->particle()->name() << "\n";
//						for( auto const& ind : ptr->lorentz_indices() )
//						{
//							std::cerr << ind << "\n";
//						}
						_amplitude.push_back(ptr->feynman_rule());
						return;
					}
					critical_error(FORMAT("Incoming antifermion line ends in outgoing line, which is not an antifermion ({}).", ptr->particle()->name()));
				}
				if( ptr->is_incoming() )
				{
					if( ptr->particle()->is_fermion() )
					{
//						std::cerr << "Amplitude-Feynman: " << ptr->particle()->name() << "\n";
//						for( auto const& ind : ptr->lorentz_indices() )
//						{
//							std::cerr << ind << "\n";
//						}
						_amplitude.push_back(ptr->feynman_rule());
						return;
					}
					critical_error(FORMAT("Incoming antifermion line ends in incoming line, which is not a fermion ({}).", ptr->particle()->name()));
				}
				trace_fermion_line(ptr, ptr->front(), start_direction);
				break;
			default:
				critical_error("trace_fermion_line starting_direction must be external.");
		}
	}


	void Feynman_Diagram::trace_fermion_line(Feynman_Graph::Edge_Ptr const& origin, Feynman_Graph::Vertex_Ptr const& vertex_ptr, Direction const& start_direction)
	{
		if( vertex_ptr == nullptr )
		{
			critical_error("trace_fermion_line: Vertex_Ptr is null.");
		}
		Feynman_Graph::Edge_Ptr next;
		for( auto const& e : vertex_ptr->all() ){
			if( e == origin ) continue;
			if( e->particle()->is_anti_fermion() || e->particle()->is_fermion()){
				if( next == nullptr ){
					next = e;
				} else{
					critical_error("Vertex has more than two fermions.");
				}
			}
		}

//		std::cerr << "Amplitude-Feynman: ";
//		for( auto const& edge_ptr : vertex_ptr->all() )
//		{
//			std::cerr << edge_ptr->particle()->name() << " ";
//			for( auto const& ind : edge_ptr->lorentz_indices() )
//			{
//				std::cerr << ind << "\n";
//			}
//		}
//		std::cerr << "\n";
		_amplitude.push_back(vertex_ptr->feynman_rule());
		trace_fermion_line(next, start_direction);
	}

	void Feynman_Diagram::iterate_spins()
	{

		for( auto& J : _spins )
		{
			auto j = J->j();
			if( J->m() > -J->j() )
			{
				J->m(J->m()-1);
				return;
			}
			else if( J->m() == -j )
			{
				J->m(j);
			}
		}
	}

	void Feynman_Diagram::iterate_indices()
	{
		for( auto& mu_ptr : _lorentz_indices )
		{
			Lorentz_Index& mu = *mu_ptr;
			++mu;
			if( mu > 0 )
			{
				return;
			}
		}
	}

	Complex Feynman_Diagram::evaluate_amplitude(Kinematics const& kin)
	{
		Complex M{0., 0.};
		std::size_t N_indices = std::pow(4, _lorentz_indices.size());
		for( std::size_t i = 0; i < N_indices; ++i )
		{
			Matrix interim(_amplitude[0](kin));
			for( std::size_t k = 1; k < _amplitude.size(); ++k )
			{
//					std::cerr << FORMAT("k: {}\n", k);
//					std::cerr << interim << "\n";
				interim *= _amplitude[k](kin);
//					std::cerr << interim << "\n\n";
			}
			try{
				M += interim.try_as_complex();
			}
			catch( Matrix::dimension_exception const& e ){
				critical_error("Invariant amplitude does not evaluate to a scalar.");
			}
			iterate_indices();
		}
		return _phase * M;
	}

/*
	double Feynman_Diagram::dsigma_dcos(Kinematics& kin)
	{
		kin.momenta.push_back(momentum(kin.sqrt_s, _graph._incoming[0]->particle()->mass(), _graph._incoming[1]->particle()->mass()));
		kin.momenta.push_back(momentum(kin.sqrt_s, _graph._outgoing[0]->particle()->mass(), _graph._outgoing[1]->particle()->mass()));
		std::cerr << "q_in: " << kin.momenta[0] << "\n";
		std::cerr << "q_out: " << kin.momenta[1] << "\n";
		Complex M{0., 0.};
		std::size_t N_indices = std::pow(4, _lorentz_indices.size());
		std::size_t N_spins = [&](){
			std::size_t n = 1;
			for( auto const& j : _spins )
			{
				n *= 2*j->j() + 1;
			}
			return n;
		}();
		for( std::size_t i = 0; i < N_indices; ++i )
		{

			for( std::size_t j = 0; j < N_spins; ++j )
			{
				Matrix interim(_amplitude[0](kin));
				for( std::size_t k = 1; k < _amplitude.size(); ++k )
				{
					interim *= _amplitude[k](kin);
				}
				try{
					M += interim.try_as_complex();
				}
				catch( Matrix::dimension_exception const& e ){
					critical_error("Invariant amplitude does not evaluate to a scalar.");
				}
				iterate_spins();
			}
			iterate_indices();
		}
		return (M * std::conj(M)).real();
	}
*/
	void Feynman_Diagram::phase(Complex phi)
	{
		if( !almost_identical(1., std::abs(phi)) )
		{
			warning("Phase angle phi has modulus != 1.");
		}
		_phase = phi;
	}

	Vertex_Manager_Ptr Feynman_Diagram::Vertex_Manager()
	{
		return _VMP;
	}

	Feynman_Diagram_Ptr operator*(Complex phi, Feynman_Diagram_Ptr& p)
	{
		p->phase(phi);
		return p;
	}

	std::vector<Particle_Ptr> Feynman_Diagram::incoming_particles() const
	{
		std::vector<Particle_Ptr> result(_graph._incoming.size());
		for( std::size_t i = 0; i < result.size(); ++i )
		{
			result[i] = _graph._incoming.at(i)->particle();
		}
		return result;
	}

	std::vector<Particle_Ptr> Feynman_Diagram::outgoing_particles() const
	{
		std::vector<Particle_Ptr> result(_graph._outgoing.size());
		for( std::size_t i = 0; i < result.size(); ++i )
		{
			result[i] = _graph._outgoing.at(i)->particle();
		}
		return result;
	}
}