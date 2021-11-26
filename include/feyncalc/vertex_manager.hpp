#ifndef FEYNCALC_VERTEX_MANAGER_HPP
#define FEYNCALC_VERTEX_MANAGER_HPP

#include <functional>
#include <map>
#include <utility>
#include <string>

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
        using Vertex_Function = std::function<Matrix()>;
        using String_Direction_Pair = std::pair<std::string, Direction>;
    private:
        std::map<
                std::string,
                std::map<std::vector<Direction>, Vertex_Function>
                > _vertex_rules;

        std::pair<std::string, std::vector<Direction>> generate_keys(std::vector<String_Direction_Pair> pairs) const;
    public:


        void add_vertex( Vertex_Function const& function, std::vector<std::pair<Particle_Ptr, Direction>> const& particles);

        Vertex_Function get_vertex(std::vector<String_Direction_Pair> const& pairs) const;
    };
}
#endif // FEYNCALC_VERTEX_MANAGER_HPP