#include "graph_edge.hpp"
#include "graph_vertex.hpp"
#include "particle.hpp"
#include "utility.hpp"

namespace Feynumeric{
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
        critical_error("Topology_Edge::feynman_rule() control structure reached invalid point. Topology_Edge is undefined.");
    }
}