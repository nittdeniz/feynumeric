#include "feynman_diagram.hpp"
#include "feynman_graph.hpp"
#include "graph_edge.hpp"
#include "graph_vertex.hpp"
#include "format.hpp"
#include "four_vector.hpp"
#include "messages.hpp"
#include "particle.hpp"
#include "utility.hpp"
#include "vertex.hpp"

#include <functional>
#include <memory>
#include <optional>

namespace Feynumeric
{
	Feynman_Graph::Feynman_Graph(Feynman_Diagram* diagram, Topology const& topology, std::vector<Particle_Ptr> const& incoming_list, std::vector<Particle_Ptr> const& virtual_list, std::vector<Particle_Ptr> const& outgoing_list)
	: _diagram(diagram)
	, _topology(topology)
	{
		validate_input(incoming_list, virtual_list, outgoing_list);
		create_graph(incoming_list, virtual_list, outgoing_list);
	}

	void Feynman_Graph::validate_input(std::vector<Particle_Ptr> const& incoming_list,
	                                   std::vector<Particle_Ptr> const& virtual_list,
	                                   std::vector<Particle_Ptr> const& outgoing_list)
	{
		if( _topology._incoming_edges.size () != incoming_list.size() )
		{
			critical_error("Feynman_Graph: Number of incoming particles does not match number of incoming edges.");
		}
		if( _topology._outgoing_edges.size() != outgoing_list.size() )
		{
			critical_error("Feynman_Graph: Number of outgoing particles does not match number of outgoing edges.");
		}
		if( _topology._virtual_edges.size() != virtual_list.size() )
		{
			critical_error("Feynman_Graph: Number of virtual particles does not match number of virtual edges.");
		}
	}

	void Feynman_Graph::create_graph(std::vector<Particle_Ptr> const& incoming_list, std::vector<Particle_Ptr> const& virtual_list, std::vector<Particle_Ptr> const& outgoing_list){
        auto get_particle_ptr = [&](std::size_t edge_id)
        {
            auto& edge = _topology._edges[edge_id];
            switch( edge.type() )
            {
                case Type::INCOMING:
                    return incoming_list[edge.from().num()];
                case Type::OUTGOING:
                    return outgoing_list[edge.to().num()];
                case Type::VIRTUAL:{
                    auto it = std::find(_topology._virtual_edges.begin(), _topology._virtual_edges.end(), edge_id);
                    return virtual_list[it - _topology._virtual_edges.begin()];
                }
                default:
                    critical_error("Invalid enum value for Direction in Feynman_Graph::create_graph.");
            }
        };
        // actual code
        for( auto const& [vertex_a, edge_list] : _topology._adjacency_map ){
            _vertices[vertex_a] = std::make_shared<Graph_Vertex>(vertex_a, _diagram);
        }
        for( auto const& [vertex_a, edge_list] : _topology._adjacency_map ){
            for( auto const& [vertex_b, edges] : edge_list )
            {
                if( vertex_a > vertex_b )
                {
                    _edges[vertex_a][vertex_b] = _edges[vertex_b][vertex_a];
                    continue;
                }
                for( auto const& edge_id : edges )
                {
                    Topology_Edge edge = _topology._edges[edge_id];
                    Edge_Ptr ptr = std::make_shared<Graph_Edge>(edge_id, _diagram, get_particle_ptr(edge_id));
                    switch( edge.type() ){
                        case Type::INCOMING:
                            _incoming.push_back(ptr);
                            break;
                        case Type::OUTGOING:
                            _outgoing.push_back(ptr);
                            break;
                        case Type::VIRTUAL:
                            _virtual.push_back(ptr);
                    }
                    if( edge.to().id() == vertex_a ){
                        if( edge.type() != Type::OUTGOING ){
                            ptr->front(_vertices[vertex_a]);
                            _vertices[vertex_a]->back(ptr);
                        }
                        if( edge.type() != Type::INCOMING ){
                            ptr->back(_vertices[vertex_b]);
                            _vertices[vertex_b]->front(ptr);
                        }
                    }
                    else
                    {
                        if( edge.type() != Type::OUTGOING )
                        {
                            ptr->front(_vertices[vertex_b]);
                            _vertices[vertex_b]->back(ptr);
                        }
                        if( edge.type() != Type::INCOMING ){
                            ptr->back(_vertices[vertex_a]);
                            _vertices[vertex_a]->front(ptr);
                        }
                    }
                    _edges[vertex_a][vertex_b].push_back(ptr);
                }
            }
        }
	}




}