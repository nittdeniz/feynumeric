#ifndef FEYNCALC_VERTEX_MANAGER_HPP
#define FEYNCALC_VERTEX_MANAGER_HPP

#include <functional>
#include <map>
#include <memory>
#include <utility>
#include <string>

#include "diagram.hpp"
#include "edge.hpp"
#include "function_types.hpp"
#include "matrix.hpp"
#include "particle.hpp"

namespace Feyncalc
{
    class Vertex_Manager
    {
    public:
        enum class Direction : unsigned int
        {
            IN = 0x01,
            OUT = 0x02,
            BOTH = 0x03
        };
        using String_Direction_Pair = std::pair<std::string, Direction>;
    private:
        std::map<
                std::string,
                std::map<std::vector<Direction>, Vertex_Function>
                > _vertex_rules;

        std::pair<std::string, std::vector<Direction>> generate_keys(std::vector<String_Direction_Pair> pairs) const;
    public:
        Vertex_Manager();
        Vertex_Manager(Vertex_Manager const& other);
        Vertex_Manager& operator=(Vertex_Manager const& other);

        void add_vertex( Vertex_Function const& function, std::vector<std::pair<Particle_Ptr, Direction>> const& particles);

        Vertex_Function get_vertex(std::vector<String_Direction_Pair> const& pairs) const;
        Vertex_Function get_vertex_function(Vertex_Id vertex_id, std::vector<Edge> const& edges) const;
    };

    inline std::string to_string(Vertex_Manager::Direction const& direction)
    {
        switch( direction )
        {
            case Vertex_Manager::Direction::IN:
                return "in";
            case Vertex_Manager::Direction::OUT:
                return "out";
            case Vertex_Manager::Direction::BOTH:
                return "both";
        }
    }

    using Vertex_Manager_Ptr = std::shared_ptr<Vertex_Manager>;
}
#endif // FEYNCALC_VERTEX_MANAGER_HPP