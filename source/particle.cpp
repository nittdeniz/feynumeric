#include "feyncalc/particle.hpp"

namespace Feyncalc
{
    Particle::Particle(string&& name, Type type, double mass, int charge, Angular_Momentum spin)
    : _name(std::move(name))
    , _type(type)
    , _mass(mass)
    , _charge(charge)
    , _spin(std::move(spin))
    {

    }
}
