#include "particle.hpp"
#include "particle_direction.hpp"

namespace Feynumeric{
    Particle_Direction::Particle_Direction(Particle_Ptr ptr, Edge_Direction dir)
            : particle(ptr)
            , direction(dir)
    {

    }

    Particle_Direction::Particle_Direction(Particle_Direction const& pd)
            : particle(pd.particle)
            , direction(pd.direction)
    {

    }

    Particle_Direction& Particle_Direction::operator=(Particle_Direction const& pd)
    {
        particle = pd.particle;
        direction = pd.direction;
        return *this;
    }

    bool operator<(Particle_Direction const& a, Particle_Direction const& b)
    {
        if( a.particle->name() == b.particle->name() )
        {
            return a.direction < b.direction;
        }
        return a.particle->name() < b.particle->name();
    }
}