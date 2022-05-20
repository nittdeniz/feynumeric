#include "format.hpp"
#include "matrix.hpp"
#include "messages.hpp"
#include "topology.hpp"

#include <algorithm>
#include <set>

namespace Feynumeric{
    Topology_Vertex::Topology_Vertex(std::string const& label)
    : _label(label)
    {
        if( label.size() < 2 ){
            critical_error("Topology::Topology_Vertex must be at least 2 characters long.");
        }
        switch( label[0] ){
            case 'i':{
                _type = Type::INCOMING;
                break;
            }
            case 'v':{
                _type = Type::VIRTUAL;
                break;
            }
            case 'o':{
                _type = Type::OUTGOING;
                break;
            }
            default:{
                critical_error(FORMAT("Topology::Topology_Vertex first character must be 'i', 'o', or 'v'. Given: {}.", label[0]));
            }
        }
        _num = std::stoull(label.substr(1));
    }

    Topology_Vertex::Topology_Vertex(const Topology_Vertex &vertex)
    : _type(vertex._type)
    , _label(vertex._label)
    , _id(vertex._id)
    , _num(vertex._num)
    {

    }

    Topology_Vertex& Topology_Vertex::operator=(const Topology_Vertex &vertex)
    {
        _type = vertex._type;
        _label = vertex._label;
        _id = vertex._id;
        _num = vertex._num;
        return *this;
    }

    std::size_t Topology_Vertex::id() const
    {
        return _id;
    }

    std::size_t Topology_Vertex::num() const
    {
        return _num;
    }

    std::string Topology_Vertex::label() const
    {
        return _label;
    }

    Type Topology_Vertex::type() const
    {
        return _type;
    }

    void Topology_Vertex::id(std::size_t i)
    {
        _id = i;
    }

    Topology_Edge::Topology_Edge(const Topology_Vertex &f, const Topology_Vertex &t)
    : _from(f)
    , _to(t)
    {
        if( _from._type == Type::INCOMING && _to._type == Type::VIRTUAL ){
            _type = Type::INCOMING;
        }
        else if( _from._type == Type::VIRTUAL && _to._type == Type::VIRTUAL ){
            _type = Type::VIRTUAL;
        }
        else if( _from._type == Type::VIRTUAL && _to._type == Type::OUTGOING ){
            _type = Type::OUTGOING;
        }
        else{
            critical_error("Topology::Topology_Edge must be pair (i, v), (v, v), or (v, o).");
        }
    }



    Topology_Edge::Topology_Edge(std::string const& from, std::string const& to)
    : Topology_Edge(Topology_Vertex(from), Topology_Vertex(to))
    {

    }

    Topology_Edge::Topology_Edge(const Topology_Edge &copy)
    : _from(copy._from)
    , _to(copy._to)
    , _type(copy._type)
    {
    }

    Topology_Edge &Topology_Edge::operator=(const Topology_Edge &copy)
    {
        _from = copy._from;
        _to = copy._to;
        _type = copy._type;
        return *this;
    }

    Topology_Vertex& Topology_Edge::from()
    {
        return _from;
    }

    Topology_Vertex& Topology_Edge::to()
    {
        return _to;
    }

    Type Topology_Edge::type() const
    {
        return _type;
    }

    void Topology::validate()
    {
        if( _edges.empty() ){
            critical_error("Topology contains no edges.");
        }
        if( _incoming_edges.empty() ){
            critical_error("Topology contains no incoming edges.");
        }
        if( _outgoing_edges.empty() ){
            critical_error("Topology contains no outgoing edges.");
        }
        {
            std::vector<std::size_t> nums;
            for( auto const &edge : _incoming_edges )
            {
                nums.push_back(_edges[edge].from().num());
            }
//            std::sort(nums.begin(), nums.end());
            for( std::size_t i = 0; i < nums.size(); ++i ){
                if( nums[i] != i ){
                    critical_error(FORMAT("Topology: Incoming index with id {} must have id {}.", nums[i], i + 1));
                }
            }
        }
        {
            std::vector<std::size_t> nums;
            for( auto const &edge : _outgoing_edges )
            {
                nums.push_back(_edges[edge].to().num());
            }
//            std::sort(nums.begin(), nums.end());
            for( std::size_t i = 0; i < nums.size(); ++i ){
                if( nums[i] != i ){
                    critical_error(FORMAT("Topology: Outgoing index with id {} must have id {}.", nums[i], i));
                }
            }
        }
        {
            std::set<std::size_t> nums;
            for( auto const &edge : _virtual_edges )
            {
                nums.insert(_edges[edge].from().num());
                nums.insert(_edges[edge].to().num());
            }
            std::size_t i{0};
            for( auto it = nums.begin(); it != nums.end(); it++, i++){
                if( (*it) != i ){
                    critical_error(FORMAT("Topology: Topology_Vertex index with id {} must have id {}.", *it, i+1));
                }
            }
        }
        std::size_t const num_vertices = _adjacency_map.size();
        std::size_t const num_edges    = _edges.size();
        Matrix M(num_vertices, num_vertices, 1);
        for( auto const& [id_a, list_a] : _adjacency_map ){
            for( auto const& [id_b, list_b] : list_a ){
                if( !list_b.empty() ){
                    M(id_a, id_b) = 1;
                }
            }
        }
        Matrix N = M;
        for( std::size_t i = 0; i < num_edges; ++i ){
            N *= M;
        }
        for( auto const& elem : N ){
            if( elem == 0. ){
                critical_error("Topology is not connected.");
            }
        }
    }

    Topology::Topology(std::vector<Topology_Edge> edge_list)
    {
        std::size_t i{0};
        std::map<std::string, std::size_t> vertex_labels;
        for( auto& edge : edge_list ){
            auto pos = _edges.size();
            switch( edge.type() ){
                case Type::INCOMING:  _incoming_edges.push_back(pos); break;
                case Type::OUTGOING:  _outgoing_edges.push_back(pos); break;
                case Type::VIRTUAL:   _virtual_edges.push_back(pos); break;
            }
            auto const label1 = edge.from().label();
            if( !vertex_labels.contains(label1) ){
                vertex_labels[label1] = i++;
            }
            auto const label2 = edge.to().label();
            if( !vertex_labels.contains(label2) ){
                vertex_labels[label2] = i++;
            }
            auto const id1 = vertex_labels[label1];
            auto const id2 = vertex_labels[label2];
            _adjacency_map[id1][id2].push_back(pos);
            _adjacency_map[id2][id1].push_back(pos);
            edge.from().id(id1);
            edge.to().id(id2);
            _edges.push_back(edge);
        }
        std::sort(_incoming_edges.begin(), _incoming_edges.end(), [&](std::size_t const lhs, std::size_t const rhs){
            return _edges[lhs].from().num() < _edges[rhs].from().num();
        });
        std::sort(_outgoing_edges.begin(), _outgoing_edges.end(), [&](std::size_t const lhs, std::size_t const rhs){
            return _edges[lhs].to().num() < _edges[rhs].to().num();
        });
        std::sort(_virtual_edges.begin(), _virtual_edges.end(), [&](std::size_t const lhs, std::size_t const rhs){
            auto min1 = std::min(_edges[lhs].from().num(), _edges[lhs].to().num());
            auto min2 = std::min(_edges[rhs].from().num(), _edges[rhs].to().num());
            if( min1 == min2 ){
                return std::max(_edges[lhs].from().num(), _edges[lhs].to().num()) < std::max(_edges[rhs].from().num(), _edges[rhs].to().num());
            }
            return min1 < min2;
        });
        validate();
    }

    Topology::Topology(const Topology &copy)
    : _edges(copy._edges)
    , _incoming_edges(copy._incoming_edges)
    , _outgoing_edges(copy._outgoing_edges)
    , _virtual_edges(copy._virtual_edges)
    , _adjacency_map(copy._adjacency_map)
    {

    }

    Topology &Topology::operator=(const Topology &copy)
    {
        _edges          = copy._edges;
        _incoming_edges = copy._incoming_edges;
        _outgoing_edges = copy._outgoing_edges;
        _virtual_edges  = copy._virtual_edges;
        _adjacency_map  = copy._adjacency_map;
        return *this;
    }
}