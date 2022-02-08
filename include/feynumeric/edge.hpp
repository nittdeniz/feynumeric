/*
#ifndef FEYNUMERIC_EDGE_HPP
#define FEYNUMERIC_EDGE_HPP

#include <functional>
#include <memory>
#include <optional>
#include <vector>

#include "angular_momentum.hpp"
#include "lorentz_index.hpp"
#include "matrix.hpp"
#include "momentum.hpp"

namespace Feynumeric
{
	class Diagram;

    class Edge;
    using Edge_Ptr = std::shared_ptr<Edge>;

    class Particle;
    using Particle_Ptr = std::shared_ptr<Particle>;

    class Edge
    {
    public:
        enum class Type
        {
            UNDEFINED,
            INCOMING,
            OUTGOING,
            VIRTUAL
        };
    private:
        std::size_t _a, _b;
        Type _type;
        Particle_Ptr _particle;
        Matrix _momentum;
        std::vector<Edge*> _neighbours;

        std::vector<Lorentz_Index_Ptr> _assigned_indices;
        Angular_Momentum_Ptr _angular_momentum;

        Diagram* _diagram;



    public:
        Edge(std::size_t a, std::size_t b, Type type=Type::UNDEFINED, Particle_Ptr particle = nullptr);
        Edge(Edge const& edge);
        Edge& operator=(Edge const& edge);

        void set_diagram(Diagram* diagram);

        std::size_t a() const;
        std::size_t b() const;

//        void a(Vertex_Id new_a);
//        void b(Vertex_Id new_b);

        void particle(Particle_Ptr new_particle);
        Particle_Ptr particle() const;

        void add_neighbour(Edge* neighbour);
        void clear_neighbours();

        std::vector<Lorentz_Index_Ptr> get_lorentz_indices() const;

        void assign_lorentz_index(Lorentz_Index_Ptr const&);
        void clear_lorentz_indices();

        void assign_angular_momentum(Angular_Momentum_Ptr const&);


        std::vector<Edge*> neighbours();

        bool is_incoming() const;
        bool is_outgoing() const;
        bool is_virtual() const;
        bool is_undefined() const;

        void momentum(Matrix const& momentum);
        Matrix momentum() const;
        Momentum_Func four_momentum() const;
        Angular_Momentum_Ptr spin() const;
        std::string to_string() const;

        std::function<Matrix(Kinematics const&)> feynman_rule();

        friend bool shares_vertex(Edge* lhs, Edge* rhs);
        friend std::optional<std::size_t> shared_vertex(Edge* lhs, Edge* rhs);

        friend std::ostream& operator<<(std::ostream& out, Edge const& edge);
        friend bool operator==(Edge const& lhs, Edge const& rhs);
        friend bool operator!=(Edge const& lhs, Edge const& rhs);

        friend class Graph;
    };

    bool shares_vertex(Edge* lhs, Edge* rhs);
    std::optional<std::size_t> shared_vertex(Edge* lhs, Edge* rhs);

    std::ostream& operator<<(std::ostream& out, Edge const& edge);
    bool operator==(Edge const& lhs, Edge const& rhs);
    bool operator!=(Edge const& lhs, Edge const& rhs);

}

#endif // FEYNUMERIC_EDGE_HPP*/