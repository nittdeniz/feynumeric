#include "diagram.hpp"
#include "particle.hpp"

#include <functional>
#include <iostream>
#include <utility.hpp>

namespace Feyncalc
{
    Diagram::Diagram(Graph const& graph, vector<Particle_Ptr> &&incoming_list, vector<Particle_Ptr> &&virtual_list,
                               vector<Particle_Ptr> &&outgoing_list)
                               : _graph(graph)
                               , _incoming_particles(std::move(incoming_list))
                               , _virtual_particles(std::move(virtual_list))
                               , _outgoing_particles(std::move(outgoing_list))
    {
        if( _incoming_particles.size() != _graph._incoming_edges.size() )
        {
            cerr << "Number of incoming particles [" << _incoming_particles.size() << "] does not match number of incoming edges [" << _graph._incoming_edges.size() << "] of graph.\n";
            abort();
        }
        if( _outgoing_particles.size() != _graph._outgoing_edges.size() )
        {
            cerr << "Number of outgoing particles [" << _outgoing_particles.size() << "] does not match number of outgoing edges [" << _graph._outgoing_edges.size() << "] of graph.\n";
            abort();
        }
        if( _virtual_particles.size() != _graph._virtual_edges.size() )
        {
            cerr << "Number of virtual particles [" << _virtual_particles.size() << "] does not match number of virtual edges [" << _graph._virtual_edges.size() << "] of graph.\n";
            abort();
        }
        
        for( size_t i = 0; i < _incoming_particles.size(); ++i ){
//            cerr << "Setting edge: " << _graph._edges[_graph._incoming_edges[i]] << "\n";
            _graph._edges[_graph._incoming_edges[i]].particle(_incoming_particles[i]);
//            cerr << "Particle: "  << _graph._edges[_graph._incoming_edges[i]].particle()->name() << "\n";
        }
        for( size_t i = 0; i < _outgoing_particles.size(); ++i )
        {
//            cerr << "Setting edge: " << _graph._edges[_graph._outgoing_edges[i]] << "\n";
            _graph._edges[_graph._outgoing_edges[i]].particle(_outgoing_particles[i]);
//            cerr << "Particle: "  << _graph._edges[_graph._outgoing_edges[i]].particle()->name() << "\n";
        }
        for( size_t i = 0; i < _virtual_particles.size(); ++i )
        {
//            cerr << "Setting edge: " << _graph._edges[_graph._virtual_edges[i]] << "\n";
            _graph._edges[_graph._virtual_edges[i]].particle(_virtual_particles[i]);
//            cerr << "Particle: "  << _graph._edges[_graph._virtual_edges[i]].particle()->name() << "\n";
        }

    }


    Diagram::Diagram(const Diagram &diagram)
    : _graph(diagram._graph)
    , _incoming_particles(diagram._incoming_particles)
    , _virtual_particles(diagram._virtual_particles)
    , _outgoing_particles(diagram._outgoing_particles)
    , _amplitude(diagram._amplitude)
    {

    }

    Diagram &Diagram::operator=(const Diagram &diagram)
    {
        _graph = diagram._graph;
        _incoming_particles = diagram._incoming_particles;
        _virtual_particles = diagram._virtual_particles;
        _outgoing_particles = diagram._outgoing_particles;
        return *this;
    }

    void Diagram::trace_fermion_line(vector<Edge>& remaining_edges, Edge const& starting_edge, Edge const& current_edge)
    {
        if( !current_edge.particle() )
        {
            cerr << "Particle not initialised " << current_edge << "\n";
            abort();
        }

//        cerr << "Tracing: " << current_edge << "\n";

        bool is_remaining = false;
        for( auto const& edge : remaining_edges )
        {
            if( edge == current_edge )
            {
//                cerr << "Erasing edge: "  << current_edge << "\n";
//                cerr << "Before:  " << remaining_edges.size() << "\n";
                is_remaining = true;
                std::erase(remaining_edges, current_edge);
//                cerr << "After: " << remaining_edges.size() << "\n";
                break;
            }

        }
        if( !is_remaining )
        {
            return;
        }
        if( current_edge.particle() == nullptr )
        {
            cerr << "Edge has no particle assigned. " << current_edge << "\n";
            abort();
        }

        std::function<bool(Particle const&)> is_same_type;
        std::function<bool(Particle const&)> is_opposite_type;
        std::function<bool(Edge const&)> is_same_direction;
        std::function<bool(Edge const&)> is_opposite_direction;

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
            cerr << "Starting Edge " << starting_edge << " ends in same direction and same type " << current_edge << ".\n";
            abort();
        }

        if( current_edge.is_incoming() )
        {
            _amplitude.push_back(current_edge.particle()->feynman_incoming);
        }
        else if( current_edge.is_virtual() )
        {
            _amplitude.push_back(current_edge.particle()->feynman_virtual);
        }
        else
        {
            _amplitude.push_back(current_edge.particle()->feynman_outgoing);
        }

        const auto neighbouring_edges = current_edge.neighbours();
//        cerr << "Neighbours: " << neighbouring_edges.size() << "\n";
        for( auto const& neighbour_id : neighbouring_edges )
        {
            if( contains(remaining_edges, _graph._edges[neighbour_id]) )
            {
                if( is_same_type(*_graph._edges[neighbour_id].particle()) )
                {
                    if( _graph._edges[neighbour_id].is_virtual() )
                    {
                        trace_fermion_line(remaining_edges, starting_edge, _graph._edges[neighbour_id]);
                        return;
                    }
                    else if((is_same_type(*_graph._edges[neighbour_id].particle()) && is_opposite_direction(_graph._edges[neighbour_id])) || (is_opposite_type(*_graph._edges[neighbour_id].particle()) && is_same_direction(_graph._edges[neighbour_id])) )
                    {
                        if( _graph._edges[neighbour_id].is_outgoing() )
                        {
                            _amplitude.push_back(_graph._edges[neighbour_id].particle()->feynman_outgoing);
                        }
                        else if( _graph._edges[neighbour_id].is_incoming() )
                        {
                            _amplitude.push_back(_graph._edges[neighbour_id].particle()->feynman_incoming );
                        }
                        else
                        {
                            cerr << "Inconsistent Edge along the path " << starting_edge << " to " << _graph._edges[neighbour_id] << " .\n";
                            abort();
                        }
                        std::erase(remaining_edges, _graph._edges[neighbour_id]);
                        return;
                    }
                    else
                    {
                        cerr << "Inconsistent Edge along the path " << starting_edge << " to " << _graph._edges[neighbour_id] << ".\n";
                        abort();
                    }
                }
            }

        }
    }

    void Diagram::generate_amplitude()
    {
        using std::abort;
        using std::cerr;

        auto remaining_edges = _graph.all_edges();
        _amplitude.clear();

        for( auto const& edge_id : _graph._outgoing_edges )
        {
            auto const& edge = _graph._edges[edge_id];
            if( is_fermion(*edge.particle()) && contains(remaining_edges, edge) )
            {
                trace_fermion_line(remaining_edges, edge, edge);
            }
        }

        for( auto const& edge_id : _graph._incoming_edges )
        {
            auto const& edge = _graph._edges[edge_id];
            if( is_anti_fermion(*edge.particle()) && contains(remaining_edges, edge) )
            {
                trace_fermion_line(remaining_edges, edge, edge);
            }
        }

        for( auto const& edge_id : _graph._virtual_edges )
        {
            auto const& edge = _graph._edges[edge_id];
            if( (is_fermion(*edge.particle()) || is_anti_fermion(*edge.particle())) && contains(remaining_edges, edge) )
            {
                cerr << "Loops are not included yet.\n";
                abort();
            }
        }

        while( remaining_edges.size() > 0 )
        {
            auto edge = remaining_edges[0];
            if( !edge.particle() )
            {
                cerr << "Particle_Ptr not initialised " << edge << "\n";
                abort();
            }
            if( is_fermion(*edge.particle()) || is_anti_fermion(*edge.particle()) )
            {
                cerr << "Remaining fermion is not part of incoming/outgoing fermion line.\n";
                cerr << *edge.particle() << " " << edge << "\n";
                abort();
            }
            if( edge.is_outgoing() )
            {
                _amplitude.push_back(edge.particle()->feynman_outgoing);
            }
            else if( edge.is_incoming() )
            {
                _amplitude.push_back(edge.particle()->feynman_incoming);
            }
            else
            {
                _amplitude.push_back(edge.particle()->feynman_virtual);
            }
            std::erase(remaining_edges, edge);
        }
    }

    void Diagram::try_generate_amplitude()
    {
        /*
        assert_diagram_validity();

        using std::cerr;
        using std::abort;


        auto remaining_edges = _graph.all_edges();

        _amplitude.clear();

        // check outgoing particles for fermions
        int index = 0;
        for( const auto& particle : _outgoing_particles )
        {
            if( is_fermion(*particle) )
            {
                const int graph_index = _graph._outgoing_vertices[index];
                const vector<Graph::Edge> edges = _graph.edges_to(graph_index);
                if( edges.size() != 1 )
                {
                    cerr << "Outgoing vertex " << graph_index << " is connected via more than one edge.\n";
                    abort();
                }
                trace_fermion_line(remaining_edges, edges.at(0), edges.at(0));
            }
            index++;
        }
        // check incoming particles for anti fermions

        index = 0;
        for( const auto& particle : _incoming_particles )
        {
            if( is_anti_fermion(*particle) )
            {
                const int graph_index = _graph._incoming_vertices[index];
                const vector<Graph::Edge> edges = _graph.edges_to(graph_index);
                if( edges.size() != 1 )
                {
                    cerr << "Incoming vertex " << graph_index << " is connected via more than one edge.\n";
                    abort();
                }
                trace_fermion_line(remaining_edges, edges.at(0), edges.at(0));
            }
            index++;
        }

        // check virtual particles for fermion loops // @TODO

        // remaining rules
        while( remaining_edges.size() > 0 )
        {
            auto edge = remaining_edges[0];
            if( is_fermion(*edge.particle) || is_anti_fermion(*edge.particle) )
            {
                cerr << "Remaining fermion is not part of incoming/outgoing fermion line.\n";
                cerr << *edge.particle << " (" << edge.a << ", " << edge.b <<")\n";
                abort();
            }
            if( _graph.is_outgoing(edge) )
            {
                _amplitude.push_back(edge.particle->feynman_outgoing);
            }
            else if( _graph.is_incoming(edge) )
            {
                _amplitude.push_back(edge.particle->feynman_incoming);
            }
            else
            {
                _amplitude.push_back(edge.particle->feynman_virtual);
            }
            std::erase(remaining_edges, edge);
            std::erase(remaining_edges, edge);
        }
         */
    }

    Complex Diagram::calculate_amplitude(const double sqrt_s, const double cos_theta) const
    {
        if( _amplitude.size() == 0 )
        {
            cerr << "Amplitude is empty.\n";
            abort();
        }
        const double sin_theta = std::sqrt(1 - cos_theta * cos_theta)+1;
        Matrix result = _amplitude.at(0)();
        for( size_t i = 1; i < _amplitude.size(); i++ )
        {
            result *= _amplitude.at(i)();
        }
        try
        {
            return sqrt_s * sin_theta*result.try_as_complex();
        }
        catch( Matrix::dimension_exception const& ex )
        {
            cerr << "Matrix Amplitude does not evaluate to a scalar.\n";
            abort();
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
