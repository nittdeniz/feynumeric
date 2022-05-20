#include <memory>

#include "feynman_diagram.hpp"
#include "graph_edge.hpp"
#include "graph_vertex.hpp"
#include "particle.hpp"
#include "utility.hpp"

namespace Feynumeric{


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

    std::vector<Particle_Direction> Graph_Vertex::particle_directions()
    {
        std::vector<Particle_Direction> result;
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

    std::string Graph_Vertex::particles_to_string() const
    {
        std::string str;
        for( auto const& e : all() )
        {
            str += e->particle()->name();
        }
        return str;
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