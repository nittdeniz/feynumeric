#ifndef FEYNUMERIC_PARTICLE_DIRECTION_HPP
#define FEYNUMERIC_PARTICLE_DIRECTION_HPP

#include "edge_direction.hpp"
#include <memory>

namespace Feynumeric{
    class Particle;
    struct Particle_Direction
    {
        using Particle_Ptr = std::shared_ptr<Particle>;
        Particle_Ptr particle;
        Edge_Direction direction;
        Particle_Direction(Particle_Ptr ptr, Edge_Direction dir = Edge_Direction::ANY);
        Particle_Direction(Particle_Direction const& pd);
        Particle_Direction& operator=(Particle_Direction const& pd);
    };
    bool operator<(Particle_Direction const& a, Particle_Direction const& b);
}

#endif //FEYNUMERIC_PARTICLE_DIRECTION_HPP
