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

	}

	Graph_Edge::Graph_Edge(std::size_t id, Feynman_Diagram* diagram, Particle_Ptr const& P)
	: _eid(id)
	  , _diagram(diagram)
	  , _particle(P)
	{
	}

	Graph_Edge::Graph_Edge(Graph_Edge const& edge)
	: std::enable_shared_from_this<Graph_Edge>(edge)
	, _eid(edge._eid)
	, _diagram(edge._diagram)
	, _particle(edge._particle)
	, _front(edge._front)
	, _back(edge._back)
	, _relative_momentum(edge._relative_momentum)
	, _spin(edge._spin)
	, _lorentz_indices(edge._lorentz_indices)
	{

	}

	Graph_Edge& Graph_Edge::operator=(Graph_Edge const& edge)
	{
		_eid = edge._eid;
		_diagram = edge._diagram;
		_particle = edge._particle;
		_front = edge._front;
		_back = edge._back;
		_relative_momentum = edge._relative_momentum;
		_spin = edge._spin;
		_lorentz_indices = edge._lorentz_indices;
		return *this;
	}

	bool Graph_Edge::is_incoming() const
	{
		return _front != nullptr && _back == nullptr;
	}

	bool Graph_Edge::is_outgoing() const
	{
		return _back != nullptr && _front == nullptr;
	}

	bool Graph_Edge::is_virtual() const
	{
		return _back != nullptr && _front != nullptr;
	}

	std::string Graph_Vertex::particles_to_string() const
	{
		std::string str;
		for( auto const& e : all() )
		{
			str += e->particle()->name();
		}
		return str;
	}

	std::shared_ptr<Graph_Vertex> Graph_Edge::front() const
	{
		return _front;
	}

	std::shared_ptr<Graph_Vertex> Graph_Edge::back() const
	{
		return _back;
	}

	Particle_Ptr Graph_Edge::particle() const
	{
		return _particle;
	}

	Four_Vector Graph_Edge::four_momentum(Kinematics const& kin) const
	{
		Four_Vector result;
		for( std::size_t i = 0; i < _relative_momentum.n_rows(); ++i )
		{
			result += _relative_momentum.at(i) * kin.momentum(i);
		}
		return result;
	}

	void Graph_Edge::front(std::shared_ptr<Graph_Vertex> const& v)
	{
		_front = v;
	}

	void Graph_Edge::back(std::shared_ptr<Graph_Vertex> const& v)
	{
		_back = v;
	}

	Matrix Graph_Edge::relative_momentum() const
	{
		return _relative_momentum;
	}

	void Graph_Edge::relative_momentum(Matrix const& momentum)
	{
		_relative_momentum = momentum;
	}

	void Graph_Edge::spin(Angular_Momentum_Ptr const& spin)
	{
		_spin = spin;
	}

	Angular_Momentum_Ptr Graph_Edge::spin() const
	{
		return _spin;
	}

	void Graph_Edge::add_lorentz_index(Lorentz_Index_Ptr const& index)
	{
		_lorentz_indices.push_back(index);
	}

	std::vector<Lorentz_Index_Ptr> Graph_Edge::lorentz_indices() const
	{
		return _lorentz_indices;
	}

	std::vector<Lorentz_Index_Ptr> Graph_Edge::lorentz_indices(std::shared_ptr<Graph_Vertex> const& ptr) const
	{
		if( is_virtual() ){
			if( contains(ptr->front(), shared_from_this()) ){
				return std::vector<Lorentz_Index_Ptr>({_lorentz_indices.begin(), _lorentz_indices.begin() + _lorentz_indices.size()/2});
			}
			return std::vector<Lorentz_Index_Ptr>({_lorentz_indices.begin() + _lorentz_indices.size()/2, _lorentz_indices.end()});
		}
		else{
			return _lorentz_indices;
		}
	}

	std::function<Matrix(Kinematics const&)> Graph_Edge::feynman_rule()
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

	Graph_Vertex::Graph_Vertex(std::size_t id, Feynman_Diagram* diagram)
	: _vid(id)
	  , _diagram(diagram)
	{

	}

	Graph_Vertex::Graph_Vertex(Graph_Vertex const& vertex)
	: std::enable_shared_from_this<Graph_Vertex>(vertex)
    , _vid(vertex._vid)
	, _diagram(vertex._diagram)
	, _front(vertex._front)
	, _back(vertex._back)
	, _feynman_rule(vertex._feynman_rule)
	{

	}

	Graph_Vertex& Graph_Vertex::operator=(Graph_Vertex const& vertex)
	{
		_vid = vertex._vid;
		_diagram = vertex._diagram;
		_front = vertex._front;
		_back = vertex._back;
		_feynman_rule = vertex._feynman_rule;
		return *this;
	}

	std::vector<std::shared_ptr<Graph_Edge>> Graph_Vertex::front() const
	{
		return _front;
	}

	std::vector<std::shared_ptr<Graph_Edge>> Graph_Vertex::back() const
	{
		return _back;
	}

	std::vector<std::shared_ptr<Graph_Edge>> Graph_Vertex::all() const
	{
		std::vector<std::shared_ptr<Graph_Edge>> result;
		result.insert(result.end(), _front.cbegin(), _front.cend());
		result.insert(result.end(), _back.cbegin(), _back.cend());
		return result;
	}

	void Graph_Edge::lorentz_indices(std::vector<Lorentz_Index_Ptr> const& list)
	{
		_lorentz_indices = list;
	}

	void Graph_Vertex::front(Edge_Ptr const& e)
	{
		_front.push_back(e);
	}

	void Graph_Vertex::back(Edge_Ptr const& e)
	{
		_back.push_back(e);
	}

	std::vector<Vertex::Particle_Direction> Graph_Vertex::particle_directions()
	{
		std::vector<Feynumeric::Vertex::Particle_Direction> result;
		result.reserve(_front.size() + _back.size());
		for( auto const& item : _back ){
			result.push_back({item->particle(), Edge_Direction::IN});
		}
		for( auto const& item : _front ){
			result.push_back({item->particle(), Edge_Direction::OUT});
		}
		std::sort(result.begin(), result.end());
		return result;
	}

	std::size_t Graph_Edge::id() const
	{
		return _eid;
	}

	std::size_t Graph_Vertex::id() const
	{
		return _vid;
	}

	std::function<Matrix(Kinematics const&)> Graph_Vertex::feynman_rule()
	{
		using namespace std::placeholders;
		auto pd = particle_directions();
		auto optional_vertex = _diagram->Vertex_Manager()->find(pd);
		if( !optional_vertex.has_value() ){
			std::stringstream stream;
			for( auto const& item : _back ){
				stream << "[" << item->particle()->name() << " (IN)]";
			}
			for( auto const& item : _front ){
				stream << "[" << item->particle()->name() << " (OUT)]";
			}
			critical_error(FORMAT("No suitable vertex found for configuration: {}", stream.str()));
		}
		auto vertex = optional_vertex.value();

		auto edge_ptrs = all();

		for( std::size_t k = 0; k < edge_ptrs.size(); ++k ){
			if( edge_ptrs[k]->is_virtual() ){
				auto indices = edge_ptrs[k]->lorentz_indices();
				auto dummy_ptr = std::make_shared<Graph_Edge>(*edge_ptrs[k]);
				// If it's a fermion we need to reverse the indices, because we're tracing the line backwards
				if( edge_ptrs[k]->particle()->is_true_fermion() ){
					if( contains(_front, edge_ptrs[k])){
						dummy_ptr->lorentz_indices({indices.begin() + indices.size() / 2, indices.end()});
					} else{
						dummy_ptr->lorentz_indices({indices.begin(), indices.begin() + indices.size() / 2});
					}
				}else{
					if( contains(_front, edge_ptrs[k])){
						dummy_ptr->lorentz_indices({indices.begin(), indices.begin() + indices.size() / 2});
					} else{
						dummy_ptr->lorentz_indices({indices.begin() + indices.size() / 2, indices.end()});
					}
				}
				_diagram->_graph._dummies.push_back(dummy_ptr);
				edge_ptrs[k] = dummy_ptr;
			}
		}

		std::vector<std::shared_ptr<Graph_Edge>> permuted;
		for( auto const& item : vertex.vertex._particle_directions ){
			for( auto const& edge : edge_ptrs ){
//				if( item.direction == Edge_Direction::IN  && edge->_front && edge->_front->id() != _vid )
//					continue;
//				if( item.direction == Edge_Direction::OUT && edge->_back  && edge->_back->id() != _vid )
//					continue;
				if( item.direction == Edge_Direction::IN && edge->_front != shared_from_this() )
					continue;
				if( item.direction == Edge_Direction::OUT && edge->_back != shared_from_this() )
					continue;
				auto current = edge->particle();
				while( current ){
					if( item.particle == current ){
						permuted.push_back(edge);
						edge_ptrs.erase(std::remove(edge_ptrs.begin(), edge_ptrs.end(), edge), edge_ptrs.end());
						goto end_loop;
					}
					current = current->parent();
				}
			}
			end_loop:;
		}
		return std::bind(vertex.vertex._vertex_function, _1, permuted);
	}
}