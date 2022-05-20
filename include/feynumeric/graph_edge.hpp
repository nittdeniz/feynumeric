#ifndef FEYNUMERIC_GRAPH_EDGE_HPP
#define FEYNUMERIC_GRAPH_EDGE_HPP

#include <memory>
#include <vector>

#include "kinematics.hpp"
#include "matrix.hpp"
#include "topology.hpp"

namespace Feynumeric
{
    class Angular_Momentum;
    class Feynman_Diagram;
    class Graph_Vertex;
    class Lorentz_Index;
    class Particle;

    class Graph_Edge : public std::enable_shared_from_this<Graph_Edge>
    {
    private:
        std::size_t _eid;
        using Angular_Momentum_Ptr = std::shared_ptr<Angular_Momentum>;
        using Lorentz_Index_Ptr = std::shared_ptr<Lorentz_Index>;
        using Particle_Ptr = std::shared_ptr<Particle>;
        using Vertex_Ptr = std::shared_ptr<Graph_Vertex>;
        Feynman_Diagram *_diagram;
        Particle_Ptr _particle;
        Vertex_Ptr _front;
        Vertex_Ptr _back;
        Matrix _relative_momentum;
        Angular_Momentum_Ptr _spin;
        std::vector<Lorentz_Index_Ptr> _lorentz_indices;
        Topology_Edge _topology_edge;

        std::function< Matrix(Kinematics const &)

        >
        _feynman_rule;
    public:
        Graph_Edge(std::size_t id, Feynman_Diagram *diagram, Particle_Ptr const &P, Topology_Edge const& topo_edge);

        Graph_Edge(Graph_Edge const &edge);

        Graph_Edge &operator=(Graph_Edge const &edge);

        bool is_incoming() const;

        bool is_outgoing() const;

        bool is_virtual() const;

        Vertex_Ptr front() const;

        Vertex_Ptr back() const;

        std::size_t id() const;

        Particle_Ptr particle() const;

        Topology_Edge topology_edge() const;

        void spin(Angular_Momentum_Ptr const &spin);

        Angular_Momentum_Ptr spin() const;

        void add_lorentz_index(Lorentz_Index_Ptr const &index);

        void lorentz_indices(std::vector<Lorentz_Index_Ptr> const &list);

        std::vector<Lorentz_Index_Ptr> lorentz_indices() const;

        std::vector<Lorentz_Index_Ptr> lorentz_indices(Vertex_Ptr const &ptr) const;

        Four_Vector four_momentum(Kinematics const &) const;

        void front(Vertex_Ptr const &v);

        void back(Vertex_Ptr const &v);

        Matrix relative_momentum() const;

        void relative_momentum(Matrix const &momentum);

        std::function< Matrix(Kinematics const &)

        >

        feynman_rule();

        friend class Graph_Vertex;
    };
}
#endif //FEYNUMERIC_GRAPH_EDGE_HPP
