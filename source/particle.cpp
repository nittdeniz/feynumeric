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

    bool Particle::is_fermion() const
    {
        return _type == Type::Particle && _spin.is_half_odd_integer();
    }

    bool Particle::is_anti_fermion() const
    {
        return _type == Type::AntiParticle && _spin.is_half_odd_integer();
    }
}
