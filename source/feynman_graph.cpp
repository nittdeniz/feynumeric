#include "feynman_diagram.hpp"
#include "feynman_graph.hpp"
#include "format.hpp"
#include "four_vector.hpp"
#include "messages.hpp"
#include "particle.hpp"
#include "utility.hpp"
#include "vertex.hpp"

#include <functional>
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
		// local helper function
		auto get_particle_ptr = [&](std::size_t edge_id)
		{
			switch( _topology._edge_list[edge_id].direction )
			{
				case Direction::INCOMING:
					return incoming_list[std::find(_topology._incoming_edges.begin(), _topology._incoming_edges.end(), edge_id) - _topology._incoming_edges.begin()];
				case Direction::OUTGOING:
					return outgoing_list[std::find(_topology._outgoing_edges.begin(), _topology._outgoing_edges.end(), edge_id) - _topology._outgoing_edges.begin()];
				case Direction::VIRTUAL:
					return virtual_list[std::find(_topology._virtual_edges.begin(), _topology._virtual_edges.end(), edge_id) - _topology._virtual_edges.begin()];
				default:
					critical_error("Invalid enum value for Direction in Feynman_Graph::create_graph.");
			}
		};

		// actual code
		for( auto const& [vertex_a, edge_list] : _topology._adjacency_map ){
			_vertices[vertex_a] = std::make_shared<Feynman_Graph::Vertex>(_diagram);
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
					Edge_Ptr ptr = std::make_shared<Feynman_Graph::Edge>(_diagram, get_particle_ptr(edge_id));
					if( _topology._edge_list[edge_id].from == vertex_a )
					{
						// add pointer shortcut
						if( contains(_topology._outgoing_edges, edge_id) )
						{
							_outgoing.push_back(ptr);
						}
						else if( contains(_topology._incoming_edges, edge_id) )
						{
							_incoming.push_back(ptr);
						}
						else{
							_virtual.push_back(ptr);
						}
						// link graph
						if( !contains(_topology._outgoing_edges, edge_id) )
						{
							ptr->front(_vertices[vertex_b]);
							_vertices[vertex_b]->back(ptr);
						}
						if( !contains(_topology._incoming_edges, edge_id) )
						{
							ptr->back(_vertices[vertex_a]);
							_vertices[vertex_a]->front(ptr);
						}
					}
					else
					{
						if( !contains(_topology._outgoing_edges, edge_id) )
						{
							ptr->front(_vertices[vertex_a]);
							_vertices[vertex_a]->back(ptr);
						}
						if( !contains(_topology._incoming_edges, edge_id) )
						{
							ptr->back(_vertices[vertex_b]);
							_vertices[vertex_b]->front(ptr);
						}
					}
					_edges[vertex_a][vertex_b].push_back(ptr);
				}
			}
		}


	}

	Feynman_Graph::Edge::Edge(Feynman_Diagram* diagram, Particle_Ptr const& P)
	: _diagram(diagram)
	, _particle(P)
	{

	}

	Feynman_Graph::Edge::Edge(Feynman_Graph::Edge const& edge)
	: std::enable_shared_from_this<Edge>(edge)
	, _diagram(edge._diagram)
	, _particle(edge._particle)
	, _front(edge._front)
	, _back(edge._front)
	, _relative_momentum(edge._relative_momentum)
	, _spin(edge._spin)
	, _lorentz_indices(edge._lorentz_indices)
	{

	}

	Feynman_Graph::Edge& Feynman_Graph::Edge::operator=(Feynman_Graph::Edge const& edge)
	{
		_diagram = edge._diagram;
		_particle = edge._particle;
		_front = edge._front;
		_back = edge._back;
		_relative_momentum = edge._relative_momentum;
		_spin = edge._spin;
		_lorentz_indices = edge._lorentz_indices;
		return *this;
	}

	bool Feynman_Graph::Edge::is_incoming() const
	{
		return _front != nullptr && _back == nullptr;
	}

	bool Feynman_Graph::Edge::is_outgoing() const
	{
		return _back != nullptr && _front == nullptr;
	}

	bool Feynman_Graph::Edge::is_virtual() const
	{
		return _back != nullptr && _front != nullptr;
	}

	std::string Feynman_Graph::Vertex::particles_to_string() const
	{
		std::string str;
		for( auto const& e : all() )
		{
			str += e->particle()->name();
		}
		return str;
	}

	Feynman_Graph::Vertex_Ptr Feynman_Graph::Edge::front() const
	{
		return _front;
	}

	Feynman_Graph::Vertex_Ptr Feynman_Graph::Edge::back() const
	{
		return _back;
	}

	Particle_Ptr Feynman_Graph::Edge::particle() const
	{
		return _particle;
	}

	Four_Vector Feynman_Graph::Edge::four_momentum(Kinematics const& kin) const
	{
		Four_Vector result;
		for( std::size_t i = 0; i < _relative_momentum.n_rows(); ++i )
		{
			result += _relative_momentum.at(i) * kin.momentum(i);
		}
		return result;
	}

	void Feynman_Graph::Edge::front(Feynman_Graph::Vertex_Ptr const& v)
	{
		_front = v;
	}

	void Feynman_Graph::Edge::back(Feynman_Graph::Vertex_Ptr const& v)
	{
		_back = v;
	}

	Matrix Feynman_Graph::Edge::relative_momentum() const
	{
		return _relative_momentum;
	}

	void Feynman_Graph::Edge::relative_momentum(Matrix const& momentum)
	{
		_relative_momentum = momentum;
	}

	void Feynman_Graph::Edge::spin(Angular_Momentum_Ptr const& spin)
	{
		_spin = spin;
	}

	Angular_Momentum_Ptr Feynman_Graph::Edge::spin() const
	{
		return _spin;
	}

	void Feynman_Graph::Edge::add_lorentz_index(Lorentz_Index_Ptr const& index)
	{
		_lorentz_indices.push_back(index);
	}

	std::vector<Lorentz_Index_Ptr> Feynman_Graph::Edge::lorentz_indices() const
	{
		return _lorentz_indices;
	}

	std::vector<Lorentz_Index_Ptr> Feynman_Graph::Edge::lorentz_indices(Feynman_Graph::Vertex_Ptr const& ptr) const
	{
		if( contains(ptr->front(), shared_from_this()) )
		{
			return std::vector<Lorentz_Index_Ptr>({_lorentz_indices.begin(), _lorentz_indices.begin() + _lorentz_indices.size()/2});
		}
		return std::vector<Lorentz_Index_Ptr>({_lorentz_indices.begin() + _lorentz_indices.size()/2, _lorentz_indices.end()});
	}

	std::function<Matrix(Kinematics const&)> Feynman_Graph::Edge::feynman_rule()
	{
		using namespace std::placeholders;
		auto e = shared_from_this();
		if( this->is_incoming() )
		{
			return std::bind(_particle->feynman_incoming, e, _1);
		}
		if( this->is_virtual() )
		{
			return std::bind(_particle->feynman_virtual, e, _1);
		}
		if( this->is_outgoing() )
		{
			return std::bind(_particle->feynman_outgoing, e, _1);
		}
		critical_error("Edge::feynman_rule() control structure reached invalid point. Edge is undefined.");
	}

	Feynman_Graph::Vertex::Vertex(Feynman_Diagram* diagram)
	: _diagram(diagram)
	{

	}

	Feynman_Graph::Vertex::Vertex(Feynman_Graph::Vertex const& vertex)
	: std::enable_shared_from_this<Feynman_Graph::Vertex>(vertex)
	, _diagram(vertex._diagram)
	{

	}

	Feynman_Graph::Vertex& Feynman_Graph::Vertex::operator=(Feynman_Graph::Vertex const& vertex)
	{
		_diagram = vertex._diagram;
		return *this;
	}

	std::vector<Feynman_Graph::Edge_Ptr> Feynman_Graph::Vertex::front() const
	{
		return _front;
	}

	std::vector<Feynman_Graph::Edge_Ptr> Feynman_Graph::Vertex::back() const
	{
		return _back;
	}

	std::vector<Feynman_Graph::Edge_Ptr> Feynman_Graph::Vertex::all() const
	{
		std::vector<Feynman_Graph::Edge_Ptr> result;
		result.insert(result.end(), _front.cbegin(), _front.cend());
		result.insert(result.end(), _back.cbegin(), _back.cend());
		return result;
	}

	void Feynman_Graph::Edge::lorentz_indices(std::vector<Lorentz_Index_Ptr> const& list)
	{
		_lorentz_indices = list;
	}

	void Feynman_Graph::Vertex::front(Edge_Ptr const& e)
	{
		_front.push_back(e);
	}

	void Feynman_Graph::Vertex::back(Edge_Ptr const& e)
	{
		_back.push_back(e);
	}

	std::size_t Feynman_Graph::Vertex::hash() const
	{
		/*
		using PD = Feynumeric::Vertex::Particle_Direction;
		std::vector<PD> lst;
		for( auto edge_ptr : _back )
		{
			lst.push_back({edge_ptr->particle(), Direction::INCOMING});
		}
		for( auto edge_ptr : _front )
		{
			lst.push_back({edge_ptr->particle(), Direction::OUTGOING});
		}
		return canonical_hash(lst);
		 */
		return 0;
	}

	std::function<Matrix(Kinematics const&)> Feynman_Graph::Vertex::feynman_rule()
	{
		using namespace std::placeholders;
		/*
		auto optional_vertex = _diagram->Vertex_Manager()->find_vertex(this);
		if( !optional_vertex.has_value() )
		{
			std::stringstream stream;
			for( auto const& item : _back )
			{
				stream << "[" << item->particle()->name() << " (IN)]";
			}
			for( auto const& item : _front )
			{
				stream << "[" << item->particle()->name() << " (OUT)]";
			}
			critical_error(FORMAT("No suitable vertex found for configuration: {}", stream.str()));
		}
		auto vertex = optional_vertex.value();

		auto edge_ptrs = all();

		for( std::size_t k = 0; k < edge_ptrs.size(); ++k )
		{
			if( edge_ptrs[k]->is_virtual() )
			{
				auto indices = edge_ptrs[k]->lorentz_indices();
				auto dummy_ptr = std::make_shared<Feynman_Graph::Edge>(*edge_ptrs[k]);
				if( contains(_front, edge_ptrs[k]) )
				{
					dummy_ptr->lorentz_indices({indices.begin(), indices.begin() + indices.size()/2});
				}
				else
				{
					dummy_ptr->lorentz_indices({indices.begin() + indices.size()/2, indices.end()});
				}
				_diagram->_graph._dummies.push_back(dummy_ptr);
				edge_ptrs[k] = dummy_ptr;
			}
		}
		return std::bind(vertex->vertex_function(), _1, vertex->sort(edge_ptrs, shared_from_this()));
		 */
		critical_error("not implemented");
	}
}