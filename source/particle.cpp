#include "feynumeric/particle.hpp"

namespace Feynumeric
{
    Particle::Particle(std::string&& name, Type type, double mass, double charge, Angular_Momentum spin)
    : _name(std::move(name))
    , _type(type)
    , _mass(mass)
    , _charge(charge)
    , _spin(std::move(spin))
    {
		_spin.m(spin.j());
    }

    std::any Particle::user_data(std::string key) const
    {
        return _user_data.at(key);
    }

    std::string Particle::name() const
    {
        return _name;
    }

    double Particle::mass() const
    {
        return _mass;
    }

    int Particle::charge() const
    {
        return _charge;
    }

    Angular_Momentum Particle::spin() const
    {
        return _spin;
    }

    Angular_Momentum Particle::isospin() const
    {
        return _isospin;
    }

    void Particle::user_data(std::string key, std::any data)
    {
        _user_data[key] = data;
    }

    bool Particle::is_fermion() const
    {
        return _type == Type::Particle
            && _spin.is_half_odd_integer();
    }

    bool Particle::is_anti_fermion() const
    {
        return _type == Particle::Type::AntiParticle
            && _spin.is_half_odd_integer();
    }

    unsigned int Particle::n_lorentz_indices() const
    {
        return static_cast<unsigned int>(_spin.j());
    }

    std::ostream &operator<<(std::ostream& out, const Particle &p)
    {
        out << p._name;
        return out;
    }

    bool is_fermion(const Particle &particle)
    {
        return particle.is_fermion();
    }

    bool is_anti_fermion(const Particle &particle)
    {
        return particle.is_anti_fermion();
    }
}
