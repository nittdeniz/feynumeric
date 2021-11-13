#ifndef FEYNCALC_PARTICLE_HPP
#define FEYNCALC_PARTICLE_HPP

#include <functional>
#include <memory>
#include <string>

#include "angular_momentum.hpp"
#include "momentum.hpp"

namespace Feyncalc
{
    using std::string;
    class Particle
    {
    public:
        enum class Type
        {
            Majorana,
            Particle,
            AntiParticle
        };
    private:
        string _name;
        Type _type;
        double _mass;
        int _charge;
        Angular_Momentum _spin;
        Angular_Momentum _isospin;
        std::function<double(double)> _width;
    public:
        Particle(string&& name, Type type, double mass = 0, int charge = 0, Angular_Momentum spin = 0);
        bool is_fermion() const;
        bool is_anti_fermion() const;
    };

    using Particle_Ptr = std::shared_ptr<Particle>;
}

#endif // FEYNCALC_PARTICLE_HPP