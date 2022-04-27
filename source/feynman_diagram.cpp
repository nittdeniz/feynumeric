#include <algorithm>
#include <iostream>
#include <utility.hpp>

#include "feynman_diagram.hpp"
#include "format.hpp"
#include "messages.hpp"
#include "particle.hpp"

namespace Feynumeric
{
	Feynman_Diagram::Feynman_Diagram(std::string&& name, Topology const& topology, Vertex_Manager_Ptr const& VMP,
	                                 std::initializer_list<Particle_Ptr> const& incoming_particles,
	                                 std::initializer_list<Particle_Ptr> const& virtual_particles,
	                                 std::initializer_list<Particle_Ptr> const& outgoing_particles)
			: _graph(this, topology, incoming_particles, virtual_particles, outgoing_particles)
			  , _VMP(VMP)
			  , _name(name){
		initialize();
		validate();
	}

	void Feynman_Diagram::add_spin(std::shared_ptr<Graph_Edge> const& edge_ptr){
		auto const spin = edge_ptr->particle()->spin();
		_spins.push_back(std::make_shared<Angular_Momentum>(spin));
		edge_ptr->spin(_spins.back());
	}

	std::string const& Feynman_Diagram::name() const{
		return _name;
	}

	void Feynman_Diagram::add_lorentz_indices(std::shared_ptr<Graph_Edge> const& edge_ptr){
		double j = edge_ptr->particle()->spin().j();
		while( j-- > 0.5 ){
			_lorentz_indices.push_back(std::make_shared<Lorentz_Index>());
			edge_ptr->add_lorentz_index(_lorentz_indices.back());
		}
	}

	Feynman_Diagram_Ptr create_diagram(std::string&& name, Topology const& topology, Vertex_Manager_Ptr const& VMP,
	                                   std::initializer_list<Particle_Ptr> const& incoming_particles,
	                                   std::initializer_list<Particle_Ptr> const& virtual_particles,
	                                   std::initializer_list<Particle_Ptr> const& outgoing_particles){
		return std::make_shared<Feynman_Diagram>(std::move(name), topology, VMP, incoming_particles, virtual_particles,
		                                         outgoing_particles);
	}

	void Feynman_Diagram::cross_outgoing(std::size_t a, std::size_t b){
		auto& A = _graph._outgoing[a];
		auto& B = _graph._outgoing[b];
		auto A_spin = A->spin();
		auto A_LI = A->lorentz_indices();
		auto A_p = A->relative_momentum();
		A->spin(B->spin());
		A->lorentz_indices(B->lorentz_indices());
		A->relative_momentum(B->relative_momentum());
		B->spin(A_spin);
		B->lorentz_indices(A_LI);
		B->relative_momentum(A_p);
		fix_internal_momenta();
	}

	void Feynman_Diagram::initialize(){
		fix_momenta();
		_spins.reserve(_graph._outgoing.size() + _graph._incoming.size());
		for( auto& edge_ptr : _graph._outgoing ){
			add_spin(edge_ptr);
			add_lorentz_indices(edge_ptr);
		}

		for( auto& edge_ptr : _graph._incoming ){
			add_spin(edge_ptr);
			add_lorentz_indices(edge_ptr);
		}
		for( auto& edge_ptr : _graph._virtual ){
			add_lorentz_indices(edge_ptr);
			add_lorentz_indices(edge_ptr);
		}
	}

	std::string Feynman_Diagram::index_to_string(Lorentz_Index_Ptr const& ptr){
		if( _indices_to_string.contains(ptr)){
			return _indices_to_string[ptr];
		}
		_indices_to_string[ptr] = FORMAT("mu{}", _indices_to_string.size());
		return _indices_to_string[ptr];
	}

	std::string Feynman_Diagram::pretty_momentum(Matrix const& relative) const{
		std::string result;
		for( std::size_t i = 0; i < relative.elements(); ++i ){
			if( relative.at(i).real() > 0 ){
				result += FORMAT("+p{}", i + 1);
			} else if( relative.at(i).real() < 0 ){
				result += FORMAT("-p{}", i + 1);
			}
		}
		return result;
	}

	void Feynman_Diagram::generate_amplitude(){
		_amplitude.clear();
		for( auto& edge_ptr : _graph._outgoing ){
			if( edge_ptr->particle()->is_true_fermion()){
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
				print_feynman_edge_rule("_u", edge_ptr);
#endif
				_amplitude.push_back(edge_ptr->feynman_rule());
				trace_fermion_line(edge_ptr, edge_ptr->back(), Direction::OUTGOING);
			} else if( !edge_ptr->particle()->is_anti_fermion()){
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
				print_feynman_edge_rule("e*", edge_ptr);
#endif
				_amplitude.push_back(edge_ptr->feynman_rule());
			}
		}
		for( auto& edge_ptr : _graph._incoming ){
			if( edge_ptr->particle()->is_anti_fermion()){
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
				print_feynman_edge_rule("_v", edge_ptr);
#endif
				_amplitude.push_back(edge_ptr->feynman_rule());
				trace_fermion_line(edge_ptr, edge_ptr->front(), Direction::INCOMING);
			} else if( !edge_ptr->particle()->is_true_fermion()){
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
				print_feynman_edge_rule("e", edge_ptr);
#endif
				_amplitude.push_back(edge_ptr->feynman_rule());
			}
		}
		for( auto& edge_ptr : _graph._virtual ){
			if( !( edge_ptr->particle()->is_true_fermion() || edge_ptr->particle()->is_anti_fermion())){
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
				print_feynman_edge_rule("D", edge_ptr);
#endif
				_amplitude.push_back(edge_ptr->feynman_rule());
			}
		}
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
		std::cout << "\n";
#endif
	}

	void Feynman_Diagram::fix_momenta(){
		fix_external_momenta();
		fix_internal_momenta();
	}

	void Feynman_Diagram::fix_external_momenta(){
		std::size_t const n_external =
				_graph._topology._incoming_edges.size() + _graph._topology._outgoing_edges.size();
		Matrix zeroes(n_external, 1);
		std::size_t position{0};

		for( auto const& edge_ptr : _graph._incoming ){
			Matrix p = zeroes;
			p[position] = 1;
			edge_ptr->relative_momentum(p);
//			_four_momenta.push_back(generate_four_momentum(Direction::INCOMING, position));
			position++;
		}
		for( auto const& edge_ptr : _graph._outgoing ){
			Matrix p = zeroes;
			p[position] = -1;
			edge_ptr->relative_momentum(p);
//			_four_momenta.push_back(generate_four_momentum(Direction::OUTGOING, position));
			position++;
		}
	}

	void Feynman_Diagram::fix_internal_momenta(){
		// reset momenta
		for( auto& e : _graph._virtual ){
			e->relative_momentum(Matrix());
		}

		auto virtuals_left = _graph._virtual;

		auto n_undetermined_momenta = [&](Feynman_Graph::Vertex_Ptr const& vertex){
			std::size_t n{0};
			std::size_t n_edges{0};
			for( auto const& e : vertex->all()){
				n_edges++;
				if( e->relative_momentum().n_cols() == 1 ){
					n++;
				}
			}
			return n_edges - n;
		};

		std::size_t const n_external =
				_graph._topology._incoming_edges.size() + _graph._topology._outgoing_edges.size();
		Matrix zeroes(n_external, 1);

		while( !virtuals_left.empty() ){
			bool changed = false;
			for( auto const& edge_ptr : virtuals_left ){
				if( edge_ptr->relative_momentum().n_cols() == 1 ){
					virtuals_left.erase(std::remove(virtuals_left.begin(), virtuals_left.end(), edge_ptr),
					                    virtuals_left.end());
					changed = true;
					break;
				}
				if( n_undetermined_momenta(edge_ptr->back()) == 1 ){
					Matrix sum_outgoing = zeroes;
					Matrix sum_incoming = zeroes;
					for( auto const& incoming_ptr : edge_ptr->back()->back()){
						sum_incoming += incoming_ptr->relative_momentum();
					}
					for( auto const& outgoing_ptr : edge_ptr->back()->front()){
						if( outgoing_ptr == edge_ptr )
							continue;
						sum_outgoing += outgoing_ptr->relative_momentum();
					}
					edge_ptr->relative_momentum(sum_outgoing + sum_incoming);
					changed = true;
					virtuals_left.erase(std::remove(virtuals_left.begin(), virtuals_left.end(), edge_ptr),
					                    virtuals_left.end());
					break;
				}
				if( n_undetermined_momenta(edge_ptr->front()) == 1 ){
					Matrix sum_outgoing = zeroes;
					Matrix sum_incoming = zeroes;
					for( auto const& incoming_ptr : edge_ptr->back()->back()){
						if( incoming_ptr == edge_ptr )
							continue;
						sum_incoming += incoming_ptr->relative_momentum();
					}
					for( auto const& outgoing_ptr : edge_ptr->front()->front()){
						sum_outgoing += outgoing_ptr->relative_momentum();
					}
					edge_ptr->relative_momentum(sum_outgoing + sum_incoming);
					changed = true;
					virtuals_left.erase(std::remove(virtuals_left.begin(), virtuals_left.end(), edge_ptr),
					                    virtuals_left.end());
					break;
				}
			}
			if( !changed ){
				critical_error("Can not fix momenta.");
			}
		}
		/*
		for( auto const& [vertex_a, map] : _graph._edges )
		{
			for( auto const& [vertex_b, edge_list] : map )
			{
				if( vertex_a > vertex_b ) continue;
				#ifdef DEBUG_AMPLITUDE
				for( auto const& edge_ptr : edge_list )
				{
					std::cerr << "Edge(" << vertex_a << ", " << vertex_b << "): " << edge_ptr->relative_momentum() << "\n";
				}
				#endif
			}
		}
		*/
	}

	void Feynman_Diagram::print_feynman_edge_rule(std::string const& id, std::shared_ptr<Graph_Edge> const& ptr, bool reverse){
		std::string momentum = pretty_momentum(ptr->relative_momentum());
		if( !ptr->is_virtual()){
			momentum.erase(0, 1);
		}
		std::cout << id << "(" << ptr->particle()->name() << ", " << momentum;
		auto lorentz_indices = ptr->lorentz_indices();
		auto print = [&](std::vector<Lorentz_Index_Ptr> const& list){
			for( auto const& item : list ){
				std::cout << ", " << index_to_string(item);
			}
		};
		if( reverse ){
			auto lhs = std::vector<Lorentz_Index_Ptr>(lorentz_indices.begin(), lorentz_indices.begin() + lorentz_indices.size()/2);
			auto rhs = std::vector<Lorentz_Index_Ptr>(lorentz_indices.begin() + lorentz_indices.size()/2, lorentz_indices.end());
			print(rhs);
			print(lhs);
		}else{
			print(lorentz_indices);
		}
		std::cout << ") ";
	}

	void Feynman_Diagram::print_feynman_vertex_rule(Feynman_Graph::Vertex_Ptr const& ptr){
		bool first = true;
		std::cout << "V(";
		for( auto const& edge_ptr : ptr->all()){
			if( !first ){
				std::cout << " // ";
			}
			std::cout << edge_ptr->particle()->name() << ", ";
			if( edge_ptr->is_virtual() && contains(ptr->front(), edge_ptr) ){
				std::cout << "-(" << pretty_momentum(edge_ptr->relative_momentum()) << ")";
			}else{
				std::cout << pretty_momentum(edge_ptr->relative_momentum());
			}

			auto lorentz_indices = edge_ptr->lorentz_indices(ptr);
			auto print = [&](std::vector<Lorentz_Index_Ptr> const& list){
				for( auto const& item : list ){
					std::cout << ", " << index_to_string(item);
				}
			};
			print(lorentz_indices);
			first = false;
		}
		std::cout << ") ";
	}

	void Feynman_Diagram::trace_fermion_line(std::shared_ptr<Graph_Edge> const& ptr, Direction const& start_direction,
	                                         Edge_Direction vertex_direction){
		if( ptr == nullptr ){
			return;
		}
		switch( start_direction ){
			case Direction::OUTGOING:
				if( ptr->is_incoming()){
					if( ptr->particle()->is_true_fermion()){
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
						print_feynman_edge_rule("u", ptr);
#endif
						_amplitude.push_back(ptr->feynman_rule());
						return;
					}
					critical_error(FORMAT("Outgoing fermion line ends in incoming line, which is not a fermion ({}).",
					                      ptr->particle()->name()));
				} else if( ptr->is_outgoing()){
					if( ptr->particle()->is_anti_fermion()){
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
						print_feynman_edge_rule("v", ptr);
#endif
						_amplitude.push_back(ptr->feynman_rule());
						return;
					}
					critical_error(
							FORMAT("Outgoing fermion line ends in outgoing line, which is not an anti fermion ({}).",
							       ptr->particle()->name()));
				} else{
					if( ptr->particle()->is_true_fermion()){
						if( vertex_direction == Edge_Direction::IN ){
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
							print_feynman_edge_rule("D", ptr, true);
#endif
							_amplitude.push_back(ptr->feynman_rule());
							trace_fermion_line(ptr, ptr->back(), start_direction);
							return;
						} else{
							critical_error(
									FORMAT("Outgoing fermion line joins intermediate fermion in wrong direction. Change to an anti fermion or reverse edge in topography."));
						}
					} else if( ptr->particle()->is_anti_fermion()){
						if( vertex_direction == Edge_Direction::OUT ){
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
							print_feynman_edge_rule("D", ptr);
#endif
							_amplitude.push_back(ptr->feynman_rule());
							trace_fermion_line(ptr, ptr->front(), start_direction);
							return;
						} else{
							critical_error(
									FORMAT("Outgoing fermion line joins intermediate anti fermion in wrong direction. Change to a fermion or reverse edge in topography."));
						}
					}
				}
				critical_error("Reached end of control flow which should not be reachable.");
			case Direction::INCOMING:
				if( ptr->is_outgoing()){
					if( ptr->particle()->is_anti_fermion()){
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
						print_feynman_edge_rule("v", ptr);
#endif
						_amplitude.push_back(ptr->feynman_rule());
						return;
					}
					critical_error(
							FORMAT("Incoming anti fermion line ends in outgoing line, which is not an anti fermion ({}).",
							       ptr->particle()->name()));
				} else if( ptr->is_incoming()){
					if( ptr->particle()->is_true_fermion()){
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
						print_feynman_edge_rule("u", ptr);
#endif
						_amplitude.push_back(ptr->feynman_rule());
						return;
					}
					critical_error(
							FORMAT("Incoming antifermion line ends in incoming line, which is not a fermion ({}).",
							       ptr->particle()->name()));
				} else{
					if( ptr->particle()->is_true_fermion()){
						if( vertex_direction == Edge_Direction::OUT ){
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
							print_feynman_edge_rule("D", ptr);
#endif
							_amplitude.push_back(ptr->feynman_rule());
							trace_fermion_line(ptr, ptr->front(), start_direction);
							return;
						} else{
							critical_error(
									FORMAT("Incoming anti fermion line joins intermediate fermion in wrong direction. Change to an anti fermion or reverse edge in topography."));
						}
					} else if( ptr->particle()->is_anti_fermion()){
						if( vertex_direction == Edge_Direction::IN ){
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
							print_feynman_edge_rule("D", ptr);
#endif
							_amplitude.push_back(ptr->feynman_rule());
							trace_fermion_line(ptr, ptr->back(), start_direction);
							return;
						} else{
							critical_error(
									FORMAT("Incoming anti fermion line joins intermediate anti fermion in wrong direction. Change to a fermion or reverse edge in topography."));
						}
					}
				}
				critical_error("Reached end of control flow which should not be reachable.");
			default:
				critical_error("trace_fermion_line starting_direction must be external.");
		}
	}


	void Feynman_Diagram::trace_fermion_line(std::shared_ptr<Graph_Edge> const& origin,
	                                         Feynman_Graph::Vertex_Ptr const& vertex_ptr,
	                                         Direction const& start_direction){
		std::shared_ptr<Graph_Edge> next;
		Edge_Direction vertex_direction;
		for( auto const& e : vertex_ptr->front()){
			if( e == origin ) continue;
			if( e->particle()->is_fermion()){
				if( next == nullptr ){
					next = e;
					vertex_direction = Edge_Direction::OUT;
				} else{
					critical_error("Vertex has more than two fermions.");
				}
			}
		}
		for( auto const& e : vertex_ptr->back()){
			if( e == origin ) continue;
			if( e->particle()->is_fermion()){
				if( next == nullptr ){
					next = e;
					vertex_direction = Edge_Direction::IN;
				} else{
					critical_error("Vertex has more than two fermions.");
				}
			}
		}
#if DEBUG_AMPLITUDE == 1 || PRINT_AMPLITUDE == 1
		print_feynman_vertex_rule(vertex_ptr);
#endif
		if( next == nullptr ){
			CRITICAL_ERROR("Vertex has no connecting fermion.");
		}
		_amplitude.push_back(vertex_ptr->feynman_rule());
		trace_fermion_line(next, start_direction, vertex_direction);
	}

	void Feynman_Diagram::iterate_spins(){
		for( auto& J : _spins ){
			++( *J );
			if( J->m() < J->j()){
				break;
			}
		}
	}

	void Feynman_Diagram::iterate_indices(){
		for( auto& mu_ptr : _lorentz_indices ){
			Lorentz_Index& mu = *mu_ptr;
			++mu;
			if( mu > 0 ){
				return;
			}
		}
	}

	void Feynman_Diagram::reset_spins(){
		for( auto& spin : _spins ){
			spin->reset();
		}
	}

	void Feynman_Diagram::reset_indices(){
		for( auto& index : _lorentz_indices ){
			index->reset();
		}
	}

	Complex Feynman_Diagram::evaluate_amplitude(Kinematics const& kin){
		Complex M{0., 0.};
		std::size_t N_indices = std::pow(4, _lorentz_indices.size());
		for( std::size_t i = 0; i < N_indices; ++i ){
#if DEBUG_AMPLITUDE == 1
			std::cerr << "indices: " << i << "\n";
			Matrix interim(_amplitude[0](kin));
			std::cerr << "@@@0:\n" << interim << "\n";
#else
			Matrix interim(_amplitude[0](kin));
#endif
			for( std::size_t k = 1; k < _amplitude.size(); ++k ){
#if DEBUG_AMPLITUDE == 1
				auto temp = _amplitude[k](kin);
				std::cerr << "@@@"<<k<<"\n" << temp << "\n";
				interim *= temp;
#else
				interim *= _amplitude[k](kin);
#endif
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

	void Feynman_Diagram::phase(Complex phi){
		if( !almost_identical(1., std::abs(phi))){
			warning("Phase angle phi has modulus != 1.");
		}
		_phase = phi;
	}

	Vertex_Manager_Ptr Feynman_Diagram::Vertex_Manager(){
		return _VMP;
	}

	Feynman_Diagram_Ptr operator*(Complex phi, Feynman_Diagram_Ptr& p){
		p->phase(phi);
		return p;
	}

	std::vector<Particle_Ptr> Feynman_Diagram::incoming_particles() const{
		std::vector<Particle_Ptr> result(_graph._incoming.size());
		for( std::size_t i = 0; i < result.size(); ++i ){
			result[i] = _graph._incoming.at(i)->particle();
		}
		return result;
	}

	std::vector<Particle_Ptr> Feynman_Diagram::outgoing_particles() const{
		std::vector<Particle_Ptr> result(_graph._outgoing.size());
		for( std::size_t i = 0; i < result.size(); ++i ){
			result[i] = _graph._outgoing.at(i)->particle();
		}
		return result;
	}

	void Feynman_Diagram::validate(){
		// incoming quantum numbers == outgoing quantum numbers
		double charge = 0;
		double isospin = 0;
		double baryon_number = 0;
		double lepton_number = 0;

		for( auto const& e : _graph._incoming ){
			charge += e->particle()->charge();
			isospin += e->particle()->isospin().m();
			baryon_number += e->particle()->baryon_number();
			lepton_number += e->particle()->lepton_number();
		}

		for( auto const& e : _graph._outgoing ){
			charge -= e->particle()->charge();
			isospin -= e->particle()->isospin().m();
			baryon_number -= e->particle()->baryon_number();
			lepton_number -= e->particle()->lepton_number();
		}

		if( charge != 0 ){
			warning("Charge is not conserved between incoming and outgoing particles.");
		}
		if( isospin != 0 ){
			warning("Isospin is not conserved between incoming and outgoing particles.");
		}
		if( baryon_number != 0 ){
			warning("Baryon number is not conserved between incoming and outgoing particles.");
		}
		if( lepton_number != 0 ){
			warning("Lepton number is not conserved between incoming and outgoing particles.");
		}

		for( auto const&[id, v] : _graph._vertices ){
			charge = 0;
			isospin = 0;
			baryon_number = 0;
			lepton_number = 0;
			for( auto const& e : v->back()){
				charge += e->particle()->charge();
				isospin += e->particle()->isospin().m();
				baryon_number += e->particle()->baryon_number();
				lepton_number += e->particle()->lepton_number();
			}
			for( auto const& e : v->front()){
				charge -= e->particle()->charge();
				isospin -= e->particle()->isospin().m();
				baryon_number -= e->particle()->baryon_number();
				lepton_number -= e->particle()->lepton_number();
			}

			std::string particles = v->particles_to_string();

			if( charge != 0 ){
				warning(FORMAT("Charge is not conserved at vertex: {}", particles));
			}
			if( isospin != 0 ){
				warning(FORMAT("Isospin is not conserved at vertex: {}", particles));
			}
			if( baryon_number != 0 ){
				warning(FORMAT("Baryon number is not conserved at vertex: {}", particles));
			}
			if( lepton_number != 0 ){
				warning(FORMAT("Lepton number is not conserved at vertex: {}", particles));
			}
		}
	}
}