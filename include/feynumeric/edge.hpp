#ifndef FEYNUMERIC_EDGE_HPP
#define FEYNUMERIC_EDGE_HPP

#include <functional>
#include <memory>
#include <optional>
#include <vector>

#include "matrix.hpp"

namespace Feynumeric
{
    class Edge;
    using Edge_Ptr = std::shared_ptr<Edge>;

    class Particle;
    using Particle_Ptr = std::shared_ptr<Particle>;

    struct Edge_Id
    {
        std::size_t id;
        Edge_Id() : id(666666){}
        explicit Edge_Id(std::size_t i) : id(i){}
        Edge_Id(Edge_Id const& other) : id(other.id){}
        Edge_Id& operator=(Edge_Id const& other){ id = other.id; return *this;}
        operator std::size_t() const
        {
            return id;
        }
        inline friend bool operator==(Edge_Id a, Edge_Id b)
        {
            return a.id == b.id;
        }
    };
    struct Vertex_Id
    {
        std::size_t id;
        Vertex_Id() : id(666666){}
        explicit Vertex_Id(std::size_t i) : id(i){}
        Vertex_Id(Vertex_Id const& other) : id(other.id){}
        Vertex_Id& operator=(Vertex_Id const& other){ id = other.id; return *this;}
        operator std::size_t() const
        {
            return id;
        }
        inline friend bool operator==(Vertex_Id a, Vertex_Id b)
        {
            return a.id == b.id;
        }
    };

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
        Vertex_Id _a, _b;
        Type _type;
        Particle_Ptr _particle;
        Matrix _momentum;
        vector<Edge_Id> _neighbour_ids;

        vector<std::size_t> _assigned_indices;
        std::size_t _angular_momentum_index;



    public:
        Edge(Vertex_Id a, Vertex_Id b = Vertex_Id{0}, Type type=Type::UNDEFINED, Particle_Ptr particle = nullptr);
        Edge(std::size_t a, std::size_t b, Type type=Type::UNDEFINED);
        Edge(Edge const& edge);
        Edge& operator=(Edge const& edge);

        Vertex_Id a() const;
        Vertex_Id b() const;

//        void a(Vertex_Id new_a);
//        void b(Vertex_Id new_b);

        void particle(Particle_Ptr new_particle);
        Particle_Ptr particle() const;

        void add_neighbour(Edge_Id neighbour);

        std::vector<std::size_t> get_lorentz_indices() const;

        void assign_lorentz_index(std::size_t id);
        void clear_lorentz_indices();

        void assign_angular_momentum(std::size_t id);


        vector<Edge_Id> neighbour_ids() const;

        bool is_incoming() const;
        bool is_outgoing() const;
        bool is_virtual() const;
        bool is_undefined() const;

        void momentum(Matrix const& momentum);
        Matrix momentum() const;
        std::string to_string() const;

        std::function<Matrix()> feynman_rule() const;

        friend bool shares_vertex(Edge const& lhs, Edge const& rhs);
        friend std::optional<Vertex_Id> shared_vertex(Edge const& lhs, Edge const& rhs);

        friend std::ostream& operator<<(std::ostream& out, Edge const& edge);
        friend bool operator==(Edge const& lhs, Edge const& rhs);
        friend bool operator!=(Edge const& lhs, Edge const& rhs);

        friend class Graph;
    };

    bool shares_vertex(Edge const& lhs, Edge const& rhs);
    std::optional<Vertex_Id> shared_vertex(Edge const& lhs, Edge const& rhs);

    std::ostream& operator<<(std::ostream& out, Edge const& edge);
    bool operator==(Edge const& lhs, Edge const& rhs);
    bool operator!=(Edge const& lhs, Edge const& rhs);

}

#endif // FEYNUMERIC_EDGE_HPP