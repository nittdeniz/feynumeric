/*
#include "diagram.hpp"
#include "messages.hpp"
#include "momentum.hpp"
#include "particle.hpp"
#include "vertex_manager.hpp"

#include <functional>
#include <iostream>
#include <utility.hpp>

namespace Feynumeric
{
    Diagram::Diagram(Vertex_Manager* vertex_manager,  Graph const& graph,  std::vector<Particle_Ptr> &&incoming_list,  std::vector<Particle_Ptr> &&virtual_list,  std::vector<Particle_Ptr> &&outgoing_list)
        : _vertex_manager(vertex_manager)
        , _graph(std::make_shared<Graph>(graph))
        , _incoming_particles(std::move(incoming_list))
        , _virtual_particles(std::move(virtual_list))
        , _outgoing_particles(std::move(outgoing_list))
    {
        if( _incoming_particles.size() != _graph->_incoming_edges.size() ){
            critical_error("Number of incoming particles [" + std::to_string(_incoming_particles.size()) + "] does not match number of incoming edges [" + std::to_string(_graph->_incoming_edges.size()) + "] of graph.");
        }
        if( _outgoing_particles.size() != _graph->_outgoing_edges.size() ){
	        critical_error("Number of outgoing particles [" + std::to_string(_outgoing_particles.size()) + "] does not match number of outgoing edges [" + std::to_string(_graph->_outgoing_edges.size()) + "] of graph.");
        }
        if( _virtual_particles.size() != _graph->_virtual_edges.size() ){
	        critical_error("Number of virtual particles [" + std::to_string(_virtual_particles.size()) + "] does not match number of virtual edges [" + std::to_string(_graph->_virtual_edges.size()) + "] of graph.");
        }

        std::size_t momentum_index = 0;
        Matrix zeroes(n_total_external(), 1);
        for( size_t i = 0; i < _incoming_particles.size(); ++i ){
            Matrix momentum = zeroes;
            momentum(0, momentum_index) = 1;
            momentum_index++;
            _graph->_incoming_edges[i]->particle(_incoming_particles[i]);
            _graph->_incoming_edges[i]->momentum(momentum);
        }
        for( size_t i = 0; i < _outgoing_particles.size(); ++i )
        {
            Matrix momentum = zeroes;
            momentum(0, momentum_index) = -1;
            momentum_index++;
            _graph->_outgoing_edges[i]->particle(_outgoing_particles[i]);
            _graph->_outgoing_edges[i]->momentum(momentum);
        }
        for( size_t i = 0; i < _virtual_particles.size(); ++i )
        {
            _graph->_virtual_edges[i]->particle(_virtual_particles[i]);
            _graph->_virtual_edges[i]->momentum(zeroes);
        }

        fix_momenta();
        generate_momentum_functions();
        link_edges_to_this();
    }

    void Diagram::register_lorentz_indices()
    {
        for( auto& edge : _graph->_edges )
        {
            edge.clear_lorentz_indices();
            int n_indices = edge.particle()->n_lorentz_indices();
            if( edge.is_virtual() )
            {
                n_indices *= 2;
            }
            while( n_indices --> 0 )
            {
	            _lorentz_indices.push_back(std::make_shared<Lorentz_Index>());
	            edge.assign_lorentz_index(_lorentz_indices.back());
            }
        }
        std::cerr << "registered indices\n";
    }

    std::size_t Diagram::n_total_external() const
    {
        return _incoming_particles.size() + _outgoing_particles.size();
    }

    void Diagram::generate_momentum_functions(){
    	if( _incoming_particles.size() == 1 )
	    {
    		if( _outgoing_particles.size() == 2 )
		    {
			    critical_error("Combination of incoming+outgoing particles is not supported [generate_momentum_functions]");
		    }
    		else
		    {
			    critical_error("Combination of incoming+outgoing particles is not supported [generate_momentum_functions]");
		    }
	    }
    	else if( _incoming_particles.size() == 2 )
	    {
    		_momenta.push_back([&](Kinematics const& kin)
		                       {
									double q = momentum(kin.sqrt_s, _incoming_particles[0]->mass(), _incoming_particles[1]->mass());
									return Four_Momentum(q, _incoming_particles[0]->mass(), 0, 0);
		                       });
		    _momenta.push_back([&](Kinematics const& kin)
		                       {
			                       double q = momentum(kin.sqrt_s, _incoming_particles[0]->mass(), _incoming_particles[1]->mass());
			                       return Four_Momentum(-q, _incoming_particles[1]->mass(),0, 0);
		                       });
			if( _outgoing_particles.size() == 2 )
			{
				_momenta.push_back([&](Kinematics const& kin)
				                   {
					                   double q = momentum(kin.sqrt_s, _outgoing_particles[0]->mass(), _outgoing_particles[1]->mass());
					                   return Four_Momentum(q, _outgoing_particles[0]->mass(), kin.cosines[0], 0);
				                   });
				_momenta.push_back([&](Kinematics const& kin)
				                   {
					                   double q = momentum(kin.sqrt_s, _outgoing_particles[0]->mass(), _outgoing_particles[1]->mass());
					                   return Four_Momentum(q, _outgoing_particles[1]->mass(), kin.cosines[0], 0);
				                   });
			}
			else
			{
				critical_error("Combination of incoming+outgoing particles is not supported [generate_momentum_functions]");
			}
	    }
    	else
	    {
    		critical_error("Number of incoming particles is not supported [generate_momentum_functions]");
	    }
    }

    void Diagram::fix_momenta()
    {
        const Matrix zeroes(n_total_external(), 1);
        const auto all_vertices = _graph->all_vertex_ids();

        auto is_incoming = [](Edge* e, std::size_t vertex_id)
        {
            return e->b() == vertex_id;
        };
        auto is_outgoing = [](Edge* e, std::size_t vertex_id)
        {
            return e->a() == vertex_id;
        };
        auto n_unknown_momenta = [&](std::size_t vertex_id)
        {
            std::size_t known_momenta{0};
            auto edges = all_vertices.at(vertex_id);
            for( auto edge : edges ){
                if( edge->momentum() != zeroes )
                {
                    known_momenta++;
                }
            }
            return edges.size() - known_momenta;
        };

        auto remaining_vertices = all_vertices;

        while( !remaining_vertices.empty() )
        {
            bool changed = false;
            for( auto const& [vertex_id, edge_ptr_list] : remaining_vertices )
            {
                auto const n = n_unknown_momenta(vertex_id);
                if( n == 0 )
                {
                    changed = true;
                    remaining_vertices.erase(vertex_id);
                    break;
                }
                if( n == 1 )
                {
                    changed = true;
                    Matrix sum_outgoing = zeroes;
                    Matrix sum_incoming = zeroes;
                    Edge* unknown_edge_ptr;
                    for( auto& edge_ptr : edge_ptr_list )
                    {
                        if( edge_ptr->momentum() == zeroes )
                        {
	                        unknown_edge_ptr = edge_ptr;
                        }
                        else if( is_incoming(edge_ptr, vertex_id) )
                        {
                            sum_incoming += edge_ptr->momentum();
                        }
                        else if( is_outgoing(edge_ptr, vertex_id) )
                        {
                            sum_outgoing += edge_ptr->momentum();
                        }
                    }
                    // p_in1 + p_in2 + ... == p_out1 + p_out2 + ...

                    if( is_incoming(unknown_edge_ptr, vertex_id) )
                    {
                        unknown_edge_ptr->momentum(-( sum_outgoing + sum_incoming));
                    }
                    else
                    {
                        unknown_edge_ptr->momentum(( sum_incoming + sum_outgoing));
                    }
                    remaining_vertices.erase(vertex_id);
                    break;
                }
            }
            if( !changed )
            {
                critical_error("Can not fix momenta.");
            }
        }
        std::cerr << "---------------------------------------\n";
        for( auto const& edge : _graph->_edges )
        {
            std::cerr << edge << ": " << edge.momentum() << "\n";
        }
    }

	void Diagram::link_edges_to_this(){
		for( auto& edge : _graph->_edges )
		{
			edge.set_diagram(this);
		}
	}

	Momentum_Func Diagram::four_momentum(Matrix const& M)
    {
		if( M.n_rows() != _incoming_particles.size() + _outgoing_particles.size() )
		{
			critical_error("Edge momentum matrix size differs from number of external particles.");
		}
		return [&](Kinematics const& kin)
		{
			Four_Momentum result;
			for( std::size_t i = 0; i < _momenta.size(); ++i )
			{
				result += M.at(i) * _momenta[i](kin);
			}
			return result;
		};
    }

    void Diagram::register_angular_momenta()
    {
    	for( auto& edge : _graph->_edges )
	    {
    		edge.set_diagram(this);
	    }

        _angular_momenta.clear();
        for( auto& edge : _graph->_outgoing_edges )
        {
        	std::cerr << "register_angular_momenta: " << *edge << "\n";
            auto const& spin = edge->particle()->spin();
            if( spin.j() > 0 )
            {
	            _angular_momenta.push_back(std::make_shared<Angular_Momentum>(edge->particle()->spin()));
	            edge->assign_angular_momentum(_angular_momenta.back());
            }
        }
	    for( auto& edge : _graph->_incoming_edges )
	    {
		    auto const spin = edge->particle()->spin();
		    if( spin.j() > 0 )
		    {
			    _angular_momenta.push_back(std::make_shared<Angular_Momentum>(edge->particle()->spin()));
			    edge->assign_angular_momentum(_angular_momenta.back());
		    }
	    }
	    std::cerr << "REGISTERED: " << _angular_momenta.size() << "\n";
    }

    void Diagram::generate_amplitude()
    {
        _remaining_edges = _graph->all_edges();
        _remaining_vertices = _graph->all_vertex_ids();

        register_lorentz_indices();
        register_angular_momenta();

        for( auto const& edge_ptr : _graph->_outgoing_edges )
        {
            if( edge_ptr->particle()->is_fermion() )
            {
	            _starting_edge_ptr = edge_ptr;
                trace_fermion_line(edge_ptr);
            }
        }
        for( auto const& edge_ptr : _graph->_incoming_edges )
        {
            if( edge_ptr->particle()->is_anti_fermion() )
            {
	            _starting_edge_ptr = edge_ptr;
                trace_fermion_line(edge_ptr);
            }
        }
        // TODO: Loops
        for( auto const& edge_ptr : _remaining_edges )
        {
            if( edge_ptr->particle()->is_fermion() || edge_ptr->particle()->is_anti_fermion() )
            {
                critical_error("There are still fermions in the remaining particles.\n");
            }
            std::cerr << "adding edge: " << *edge_ptr << " with particle " << edge_ptr->particle()->name() << "\n";
            std::cerr << "Lorentz Indices: ";
            for( auto lorentz : edge_ptr->get_lorentz_indices() )
            {
                std::cerr << lorentz << " ";
            }
            std::cerr << "\n";
            std::cerr << "_amplitude.push_back: " << *edge_ptr << "\n";
            _amplitude.push_back(edge_ptr->feynman_rule());

        }
        for( auto const& [vertex_id, edge_list] : _remaining_vertices )
        {
            add_vertex(vertex_id);
        }
    }

    Diagram::Diagram(const Diagram &diagram)
    : std::enable_shared_from_this<Diagram>(diagram)
    , _vertex_manager(diagram._vertex_manager)
    , _graph(diagram._graph)
    , _incoming_particles(diagram._incoming_particles)
    , _virtual_particles(diagram._virtual_particles)
    , _outgoing_particles(diagram._outgoing_particles)
    , _amplitude(diagram._amplitude)
    , _lorentz_indices(diagram._lorentz_indices)
    , _angular_momenta(diagram._angular_momenta)
    {
	    link_edges_to_this();
    }

    Diagram &Diagram::operator=(const Diagram &diagram)
    {
        _vertex_manager = diagram._vertex_manager;
        _graph = diagram._graph;
        _incoming_particles = diagram._incoming_particles;
        _virtual_particles = diagram._virtual_particles;
        _outgoing_particles = diagram._outgoing_particles;
	    _angular_momenta = diagram._angular_momenta;
	    _lorentz_indices = diagram._lorentz_indices;
	    link_edges_to_this();
        return *this;
    }

    void Diagram::trace_fermion_line(Edge* current_edge_ptr)
    {
        using std::abort;
        using std::cerr;
        if( !current_edge_ptr->particle() )
        {
            critical_error("Particle is not initialized. Edge: " + current_edge_ptr->to_string());
        }
        const auto edge_it = std::find(_remaining_edges.begin(), _remaining_edges.end(), current_edge_ptr);
        if( edge_it == _remaining_edges.end() )
        {
            return;
        }
        _remaining_edges.erase(edge_it);

        std::function<bool(Particle const&)> is_same_type;
        std::function<bool(Particle const&)> is_opposite_type;
        std::function<bool(Edge const&)> is_same_direction;
        std::function<bool(Edge const&)> is_opposite_direction;


        if( _starting_edge_ptr->is_outgoing() )
        {
            is_same_type = is_fermion;
            is_opposite_type = is_anti_fermion;
            is_same_direction = [&](Edge const& edge){return edge.is_outgoing();};
            is_opposite_direction = [&](Edge const& edge){return edge.is_incoming();};
        }
        else if( _starting_edge_ptr->is_incoming() )
        {
            is_same_type = is_anti_fermion;
            is_opposite_type = is_fermion;
            is_same_direction = [&](Edge const& edge){return edge.is_incoming();};
            is_opposite_direction = [&](Edge const& edge){return edge.is_outgoing();};
        }

        if( current_edge_ptr != _starting_edge_ptr && is_same_direction(*current_edge_ptr) && is_same_type(*current_edge_ptr->particle()) )
        {
            critical_error("Starting Edge " + _starting_edge_ptr->to_string() + " ends in same direction and same type as " + current_edge_ptr->to_string() + ".");
        }

        std::cerr << "adding edge: " << *current_edge_ptr << " with particle " << current_edge_ptr->particle()->name() << "\n";
        std::cerr << "Lorentz Indices: ";
        for( auto lorentz : current_edge_ptr->get_lorentz_indices() )
        {
            std::cerr << lorentz << " ";
        }
        std::cerr << "_amplitude.push_back: " << *current_edge_ptr << "\n";
        _amplitude.push_back(current_edge_ptr->feynman_rule());

        auto neighbouring_edges = current_edge_ptr->neighbours();
        for( auto& neighbour_ptr : neighbouring_edges )
        {
            if( contains(_remaining_edges, neighbour_ptr) )
            {
                if( is_same_type(*neighbour_ptr->particle()) && neighbour_ptr->is_virtual() )
                {
                    add_vertex(current_edge_ptr, neighbour_ptr);
                    trace_fermion_line(neighbour_ptr);
                    return;
                }
                else if(   ( is_same_type(*neighbour_ptr->particle()) && is_opposite_direction(*neighbour_ptr) )
                        || ( is_opposite_type(*neighbour_ptr->particle()) && is_same_direction(*neighbour_ptr) )
                )
                {
                    if( !neighbour_ptr->is_outgoing() && !neighbour_ptr->is_incoming() )
                    {
                        critical_error("Inconsistent Edge along the path " + _starting_edge_ptr->to_string() + " to " + neighbour_ptr->to_string() + ".");
                    }
                    add_vertex(current_edge_ptr, neighbour_ptr);
                    std::cerr << "_amplitude.push_back: " << *neighbour_ptr << "\n";
                    _amplitude.push_back(neighbour_ptr->feynman_rule());
                    std::cerr << "adding edge: " << *neighbour_ptr << " with particle " << neighbour_ptr->particle()->name() << "\n";
                    std::cerr << "Lorentz Indices: ";
                    for( auto lorentz : neighbour_ptr->get_lorentz_indices() )
                    {
                        std::cerr << lorentz << " ";
                    }
                    std::erase(_remaining_edges, neighbour_ptr);
                    return;
                }
                else
                {
                    critical_error("Inconsistent Edge along the path " + _starting_edge_ptr->to_string() + " to " + neighbour_ptr->to_string() + ".");
                }
            }
        }
    }

    void Diagram::add_vertex(Edge* a, Edge* b)
    {
        auto optional_vertex = shared_vertex(a, b);
        if( !optional_vertex.has_value() ){
            critical_error("Edges " + a->to_string() + " and " + b->to_string() + " do not share a vertex.");
        }
        add_vertex(optional_vertex.value());
    }

    void Diagram::add_vertex(std::size_t vertex_id)
    {
        if( !_remaining_vertices.contains(vertex_id) )
        {
            return;

        }
        std::cerr << "Adding Vertex\nID: " << vertex_id << "\n";

        auto const& vertex = _remaining_vertices[vertex_id];

        std::vector<Edge*> vertex_edges;
        vertex_edges.reserve(vertex.size());
        for( auto const& edge_ptr : vertex )
        {
            vertex_edges.push_back(edge_ptr);
        }

        auto vf = _vertex_manager->get_vertex_function(vertex_id, vertex_edges);
        auto vertex_func = std::bind(vf, this, vertex_edges);

        std::cerr << "_amplitude.push_back: " << vertex_id << "\n";
        _amplitude.push_back(vertex_func);
        _remaining_vertices.erase(vertex_id);
    }

    void Diagram::iterate_spins()
    {

		for( auto& J : _angular_momenta )
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

    void Diagram::iterate_indices()
    {
    	for( auto& mu_ptr : _lorentz_indices )
	    {
    		Lorentz_Index& mu = *mu_ptr;
		    ++mu;
    		if( mu > 0 )
		    {
				break;
		    }
	    }
    }

    Complex Diagram::calculate_amplitude(const double sqrt_s, const double cos_theta)
    {
        if( _amplitude.size() == 0 )
        {
            critical_error("Amplitude is empty.");
        }
        std::cerr << "Amplitude Size: " << _amplitude.size() << "\n";
        const double sin_theta = std::sqrt(1 - cos_theta * cos_theta)+1;

        auto spin_configurations = [](std::list<Angular_Momentum_Ptr> const& lst)
        {
			std::size_t result = 1;
			for( auto const& item : lst )
			{
				result *= (2*item->j()+1);
			}
			return result;
        };

        auto index_configurations = [](std::list<Lorentz_Index_Ptr> const& lst)
        {
        	return std::pow(4, lst.size());
        };

        std::size_t N_spins = spin_configurations(_angular_momenta);
        double temp_indices = index_configurations(_lorentz_indices);
        std::size_t N_indices = static_cast<std::size_t>(temp_indices);
        if( N_indices != temp_indices )
        {
        	critical_error("Overflow: Too many lorentz indices.\n");
        }

        Complex result{0};

        Kinematics kin;

        for( std::size_t index_it = 0; index_it < N_indices; ++index_it ){
	        for( std::size_t spin_it = 0; spin_it < N_spins; ++spin_it ){
		        Matrix interim_result = _amplitude[0](kin);
		        for( size_t i = 1; i < _amplitude.size(); i++ ){
			        interim_result *= _amplitude[i](kin);
		        }
		        try{
			        result += interim_result.try_as_complex();
		        }
		        catch( Matrix::dimension_exception const& ex ){
			        critical_error("Matrix Amplitude does not evaluate to a scalar.");
		        }
		        iterate_spins();
	        }
	        iterate_indices();
        }
        return sqrt_s * sin_theta * result;
    }

    #ifdef CATCH2_TESTING_ENABLED
    bool Diagram::assert_diagram_validity() const
    #else
    void Diagram::assert_diagram_validity() const
    #endif
    {
        #ifdef CATCH2_TESTING_ENABLED
        return true;
        #endif
    }

}
*/