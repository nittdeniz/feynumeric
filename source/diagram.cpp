#include "diagram.hpp"
#include "messages.hpp"
#include "particle.hpp"
#include "vertex_manager.hpp"

#include <functional>
#include <iostream>
#include <utility.hpp>

namespace Feynumeric
{
    Diagram::Diagram(Vertex_Manager_Ptr const& vertex_manager, Graph const& graph, vector<Particle_Ptr> &&incoming_list, vector<Particle_Ptr> &&virtual_list,
                               vector<Particle_Ptr> &&outgoing_list)
                               : _vertex_manager(vertex_manager)
                               , _graph(graph)
                               , _incoming_particles(std::move(incoming_list))
                               , _virtual_particles(std::move(virtual_list))
                               , _outgoing_particles(std::move(outgoing_list))
    {
        using std::abort;
        using std::cerr;
        if( _incoming_particles.size() != _graph._incoming_edge_ids.size() )
        {
            cerr << "Number of incoming particles [" << _incoming_particles.size() << "] does not match number of incoming edges [" << _graph._incoming_edge_ids.size() << "] of graph.\n";
            abort();
        }
        if( _outgoing_particles.size() != _graph._outgoing_edge_ids.size() )
        {
            cerr << "Number of outgoing particles [" << _outgoing_particles.size() << "] does not match number of outgoing edges [" << _graph._outgoing_edge_ids.size() << "] of graph.\n";
            abort();
        }
        if( _virtual_particles.size() != _graph._virtual_edge_ids.size() )
        {
            cerr << "Number of virtual particles [" << _virtual_particles.size() << "] does not match number of virtual edges [" << _graph._virtual_edge_ids.size() << "] of graph.\n";
            abort();
        }

        std::size_t momentum_index = 0;
        Matrix zeroes(n_total_external(), 1);
        for( size_t i = 0; i < _incoming_particles.size(); ++i ){
            Matrix momentum = zeroes;
            momentum(0, momentum_index) = 1;
            momentum_index++;
            _graph._edges[_graph._incoming_edge_ids[i]].particle(_incoming_particles[i]);
            _graph._edges[_graph._incoming_edge_ids[i]].momentum(momentum);
        }
        for( size_t i = 0; i < _outgoing_particles.size(); ++i )
        {
            Matrix momentum = zeroes;
            momentum(0, momentum_index) = -1;
            momentum_index++;
            _graph._edges[_graph._outgoing_edge_ids[i]].particle(_outgoing_particles[i]);
            _graph._edges[_graph._outgoing_edge_ids[i]].momentum(momentum);
        }
        for( size_t i = 0; i < _virtual_particles.size(); ++i )
        {
            _graph._edges[_graph._virtual_edge_ids[i]].particle(_virtual_particles[i]);
            _graph._edges[_graph._virtual_edge_ids[i]].momentum(zeroes);
        }

        fix_momenta();
    }

    void Diagram::register_lorentz_indices()
    {
        auto n_lorentz_index_pos = _lorentz_indices.size();
        for( auto& edge : _graph._edges )
        {
            edge.clear_lorentz_indices();
            int n_indices = edge.particle()->n_lorentz_indices();
            if( edge.is_virtual() )
            {
                n_indices *= 2;
            }
            while( n_indices --> 0 )
            {
               _lorentz_indices.emplace_back();
               edge.assign_lorentz_index(n_lorentz_index_pos);
               ++n_lorentz_index_pos;
            }
        }
    }

    std::size_t Diagram::n_total_external() const
    {
        return _incoming_particles.size() + _outgoing_particles.size();
    }

    void Diagram::fix_momenta()
    {
        const Matrix zeroes(n_total_external(), 1);
        const auto all_vertices = _graph.all_vertex_ids();

        auto is_incoming = [](Edge e, Vertex_Id vid)
        {
            return e.b() == vid;
        };
        auto is_outgoing = [](Edge e, Vertex_Id vid)
        {
            return e.a() == vid;
        };
        auto n_unknown_momenta = [&](Vertex_Id vid)
        {
            std::size_t known_momenta{0};
            auto edge_ids = all_vertices.at(vid);
            for( auto edge_id : edge_ids ){
                auto const &edge = _graph._edges[edge_id];
                if( edge.momentum() != zeroes )
                {
                    known_momenta++;
                }
            }
            return edge_ids.size() - known_momenta;
        };

        auto remaining_vertices = all_vertices;

        while( !remaining_vertices.empty() )
        {
            bool changed = false;
            for( auto pair : remaining_vertices )
            {
                auto const& vertex_id = pair.first;
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
                    Edge_Id unknown_id;
                    auto const& edge_id_list = pair.second;
                    for( auto edge_id : edge_id_list )
                    {
                        auto const& edge = _graph._edges[edge_id];
                        if( edge.momentum() == zeroes )
                        {
                            unknown_id = edge_id;
                        }
                        else if( is_incoming(edge, vertex_id) )
                        {
                            sum_incoming += edge.momentum();
                        }
                        else if( is_outgoing(edge, vertex_id) )
                        {
                            sum_outgoing += edge.momentum();
                        }
                    }
                    // p_in1 + p_in2 + ... == p_out1 + p_out2 + ...
                    auto& unknown_edge = _graph._edges[unknown_id];
                    if( is_incoming(unknown_edge, vertex_id) )
                    {
                        unknown_edge.momentum(-(sum_outgoing + sum_incoming));
                    }
                    else
                    {
                        unknown_edge.momentum((sum_incoming + sum_outgoing));
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
        for( auto const& edge : _graph._edges )
        {
            std::cerr << edge << ": " << edge.momentum() << "\n";
        }
    }

    void Diagram::register_angular_momenta()
    {
        _angular_momenta.clear();
        std::size_t angular_momentum_index = 0;
        for( auto const& edge_id : _graph._outgoing_edge_ids )
        {
            auto& edge = _graph._edges[edge_id];
            auto const spin = edge.particle()->spin();
            if( spin.j() > 0 )
            {
                _angular_momenta.emplace_back(edge.particle()->spin());
                edge.assign_angular_momentum(angular_momentum_index);
                ++angular_momentum_index;
            }
        }
    }

    void Diagram::generate_amplitude()
    {
        _remaining_edge_ids = _graph.all_edge_ids();
        _remaining_vertices = _graph.all_vertex_ids();

        register_lorentz_indices();
        register_angular_momenta();

        for( auto const& edge_id : _graph._outgoing_edge_ids )
        {
            auto const& edge = _graph._edges[edge_id];
            if( edge.particle()->is_fermion() )
            {
                _starting_edge_id = edge_id;
                trace_fermion_line(edge_id);
            }
        }
        for( auto const& edge_id : _graph._incoming_edge_ids )
        {
            auto const& edge = _graph._edges[edge_id];
            if( edge.particle()->is_anti_fermion() )
            {
                _starting_edge_id = edge_id;
                trace_fermion_line(edge_id);
            }
        }
        // TODO: Loops
        for( auto const& edge_id : _remaining_edge_ids )
        {
            auto const& edge = _graph._edges[edge_id];
            if( edge.particle()->is_fermion() || edge.particle()->is_anti_fermion() )
            {
                critical_error("There are still fermions in the remaining particles.\n");
            }
            std::cerr << "adding edge: " << edge << " with particle " << edge.particle()->name() << "\n";
            std::cerr << "Lorentz Indices: ";
            for( auto lorentz : edge.get_lorentz_indices() )
            {
                std::cerr << lorentz << " ";
            }
            std::cerr << "\n";
            _amplitude.push_back(edge.feynman_rule());
        }
        for( auto const& vertex_id : _remaining_vertices )
        {
            add_vertex(vertex_id.first);
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
    {

    }

    Diagram &Diagram::operator=(const Diagram &diagram)
    {
        _vertex_manager = diagram._vertex_manager;
        _graph = diagram._graph;
        _incoming_particles = diagram._incoming_particles;
        _virtual_particles = diagram._virtual_particles;
        _outgoing_particles = diagram._outgoing_particles;
        return *this;
    }

    void Diagram::trace_fermion_line(Edge_Id current_edge_id)
    {
        using std::abort;
        using std::cerr;
        auto const& current_edge = _graph._edges[current_edge_id];
        if( !current_edge.particle() )
        {
            critical_error("Particle is not initialized. Edge: " + current_edge.to_string());
        }
        const auto edge_it = std::find(_remaining_edge_ids.begin(), _remaining_edge_ids.end(), current_edge_id);
        if( edge_it == _remaining_edge_ids.end() )
        {
            return;
        }
        _remaining_edge_ids.erase(edge_it);

        std::function<bool(Particle const&)> is_same_type;
        std::function<bool(Particle const&)> is_opposite_type;
        std::function<bool(Edge const&)> is_same_direction;
        std::function<bool(Edge const&)> is_opposite_direction;

        auto const& starting_edge = _graph._edges[_starting_edge_id];

        if( starting_edge.is_outgoing() )
        {
            is_same_type = is_fermion;
            is_opposite_type = is_anti_fermion;
            is_same_direction = [&](Edge const& edge){return edge.is_outgoing();};
            is_opposite_direction = [&](Edge const& edge){return edge.is_incoming();};
        }
        else if( starting_edge.is_incoming() )
        {
            is_same_type = is_anti_fermion;
            is_opposite_type = is_fermion;
            is_same_direction = [&](Edge const& edge){return edge.is_incoming();};
            is_opposite_direction = [&](Edge const& edge){return edge.is_outgoing();};
        }

        if( current_edge != starting_edge && is_same_direction(current_edge) && is_same_type(*current_edge.particle()) )
        {
            critical_error("Starting Edge " + starting_edge.to_string() + " ends in same direction and same type as " + current_edge.to_string() + ".");
        }

        std::cerr << "adding edge: " << current_edge << " with particle " << current_edge.particle()->name() << "\n";
        std::cerr << "Lorentz Indices: ";
        for( auto lorentz : current_edge.get_lorentz_indices() )
        {
            std::cerr << lorentz << " ";
        }
        _amplitude.push_back(current_edge.feynman_rule());

        const auto neighbouring_edges = current_edge.neighbour_ids();
        for( auto const& neighbour_id : neighbouring_edges )
        {
            auto const& neighbour = _graph._edges[neighbour_id];

            if( contains(_remaining_edge_ids, neighbour_id) )
            {
                if( is_same_type(*neighbour.particle())  &&  neighbour.is_virtual() )
                {
                    add_vertex(current_edge_id, neighbour_id);
                    trace_fermion_line(neighbour_id);
                    return;
                }
                else if(   ( is_same_type(*neighbour.particle()) &&  is_opposite_direction(neighbour) )
                        || ( is_opposite_type(*neighbour.particle()) &&  is_same_direction(neighbour) )
                )
                {
                    if( !neighbour.is_outgoing() && !neighbour.is_incoming() )
                    {
                        critical_error("Inconsistent Edge along the path " + starting_edge.to_string() + " to " + neighbour.to_string() + ".");
                    }
                    add_vertex(current_edge_id, neighbour_id);
                    _amplitude.push_back(neighbour.feynman_rule());
                    std::cerr << "adding edge: " << neighbour << " with particle " << neighbour.particle()->name() << "\n";
                    std::cerr << "Lorentz Indices: ";
                    for( auto lorentz : neighbour.get_lorentz_indices() )
                    {
                        std::cerr << lorentz << " ";
                    }
                    std::erase(_remaining_edge_ids, neighbour_id);
                    return;
                }
                else
                {
                    critical_error("Inconsistent Edge along the path " + starting_edge.to_string() + " to " + neighbour.to_string() + ".");
                }
            }
        }
    }

    void Diagram::add_vertex(Edge_Id a_id, Edge_Id b_id)
    {
        auto const &a = _graph._edges[a_id];
        auto const &b = _graph._edges[b_id];

        auto optional_vertex = shared_vertex(a, b);
        if( !optional_vertex.has_value() ){
            critical_error("Edges " + a.to_string() + " and " + b.to_string() + " do not share a vertex.");
        }
        add_vertex(optional_vertex.value());
    }

    void Diagram::add_vertex(Vertex_Id vertex_id)
    {
        if( !_remaining_vertices.contains(vertex_id) )
        {
            return;

        }
        std::cerr << "Adding Vertex\nID: " << static_cast<std::size_t>(vertex_id) << "\n";

        auto const& vertex = _remaining_vertices[vertex_id];

        std::vector<Edge> vertex_edges;
        vertex_edges.reserve(vertex.size());
        for( auto const& edge_id : vertex )
        {
            vertex_edges.push_back(_graph._edges[edge_id]);
        }

        auto vertex_func = std::bind(_vertex_manager->get_vertex_function(vertex_id, vertex_edges),
                                     this, vertex_edges);

        _amplitude.push_back(vertex_func);
        _remaining_vertices.erase(vertex_id);
    }


    Complex Diagram::calculate_amplitude(const double sqrt_s, const double cos_theta) const
    {
        if( _amplitude.size() == 0 )
        {
            critical_error("Amplitude is empty.");
        }
        const double sin_theta = std::sqrt(1 - cos_theta * cos_theta)+1;
        Matrix result = _amplitude.at(0)();
        for( size_t i = 1; i < _amplitude.size(); i++ )
        {
            result *= _amplitude.at(i)();
        }
        try
        {
            std::cerr << "sqrt_s: " << sqrt_s << "\n";
            std::cerr << "sin_theta: " << sin_theta << "\n";
            std::cerr << "result: " << result.try_as_complex() << "\n";
            return sqrt_s * sin_theta*result.try_as_complex();
        }
        catch( Matrix::dimension_exception const& ex )
        {
            critical_error("Matrix Amplitude does not evaluate to a scalar.");
        }
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
