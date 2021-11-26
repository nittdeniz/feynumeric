#ifndef FEYNCALC_EDGE_HPP
#define FEYNCALC_EDGE_HPP

#include <vector>

#include "particle.hpp"

namespace Feyncalc
{
    using std::size_t;
    using std::vector;
    class Edge;
    using Edge_Ptr = std::shared_ptr<Edge>;
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
        size_t _a, _b;
        Type _type;
        Particle_Ptr _particle;
        vector<int> _momentum;
        vector<size_t> _neighbours;
    public:
        Edge(size_t a = 0, size_t b = 0, Type type=Type::UNDEFINED, Particle_Ptr particle = nullptr, vector<int> momentum = {});
        Edge(Edge const& edge);
        Edge& operator=(Edge const& edge);

        size_t a() const;
        size_t b() const;

        void a(size_t new_a);
        void b(size_t new_b);

        Edge sorted() const;

        void type(Type new_type);

        void particle(Particle_Ptr new_particle);
        Particle_Ptr particle() const;

        void add_neighbour(size_t neighbour);

        vector<size_t> neighbours() const;

        bool is_incoming() const;
        bool is_outgoing() const;
        bool is_virtual() const;
        bool is_undefined() const;

        Matrix momentum(vector<Matrix> const& momenta) const;

        friend bool shares_vertex(Edge const& lhs, Edge const& rhs);

        friend std::ostream& operator<<(std::ostream& out, Edge const& edge);
        friend bool operator==(Edge const& lhs, Edge const& rhs);
        friend bool operator!=(Edge const& lhs, Edge const& rhs);
    };

    bool shares_vertex(Edge const& lhs, Edge const& rhs);

    std::ostream& operator<<(std::ostream& out, Edge const& edge);
    bool operator==(Edge const& lhs, Edge const& rhs);
    bool operator!=(Edge const& lhs, Edge const& rhs);

}

#endif // FEYNCALC_EDGE_HPP