#ifndef FEYNUMERIC_GRAPH_VERTEX_HPP
#define FEYNUMERIC_GRAPH_VERTEX_HPP

#include <functional>
#include <memory>
#include <vector>

#include "kinematics.hpp"
#include "matrix.hpp"
#include "particle_direction.hpp"


namespace Feynumeric{
    class Feynman_Diagram;
    class Graph_Edge;

    class Graph_Vertex : public std::enable_shared_from_this<Graph_Vertex>
    {
    private:
        using Edge_Ptr = std::shared_ptr<Graph_Edge>;
        std::size_t _vid;
        Feynman_Diagram* _diagram;
        std::vector<Edge_Ptr> _front;
        std::vector<Edge_Ptr> _back;
        std::function<Matrix()> _feynman_rule;
        std::vector<Particle_Direction> particle_directions();

    public:
        explicit Graph_Vertex(std::size_t id, Feynman_Diagram* diagram);
        Graph_Vertex(Graph_Vertex const& vertex);
        Graph_Vertex& operator=(Graph_Vertex const& vertex);

        std::vector<Edge_Ptr> front() const;
        std::vector<Edge_Ptr> back() const;
        std::vector<Edge_Ptr> all() const;

        std::size_t id() const;

        void front(Edge_Ptr const& e);
        void back(Edge_Ptr const& e);

        std::string particles_to_string() const;
        std::function<Matrix(Kinematics const& kin)> feynman_rule();
    };
}
#endif //FEYNUMERIC_GRAPH_VERTEX_HPP
