#ifndef FEYNCALC_VERTEX_MANAGER_HPP
#define FEYNCALC_VERTEX_MANAGER_HPP

#include <functional>
#include <map>
#include <string>

#include "matrix.hpp"
#include "particle.hpp"

namespace Feyncalc
{
    class Vertex_Manager
    {
        using Vertex_Function = std::function<Matrix()>;
    private:
        std::map<std::string, Vertex_Function> _vertex_rules;
    public:
        void add_vertex(std::vector<Particle_Ptr>&& particles, Vertex_Function&& function);

        Vertex_Function get_vertex(std::vector<Particle_Ptr> const& particles);
    };
}
#endif // FEYNCALC_VERTEX_MANAGER_HPP