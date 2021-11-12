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
    private:
        string _name;
        double mass;
        std::function<double(double)> _width;
        int _charge;
        Angular_Momentum _spin;
        Angular_Momentum _isospin;
    public:
    };

    using Particle_Ptr = std::shared_ptr<Particle>;
}

#endif // FEYNCALC_PARTICLE_HPP