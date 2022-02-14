#ifndef Feynumeric_KINEMATICS_HPP
#define Feynumeric_KINEMATICS_HPP

#include "momentum.hpp"
#include <vector>

namespace Feynumeric
{
    class Kinematics
    {
    private:
        double _sqrt_s;
        std::size_t _n_in, _n_out;
        std::vector<Four_Momentum> _momenta;
    public:
        Kinematics(double sqrt_s, std::size_t n_in, std::size_t n_out);
        double sqrt_s() const;
        Four_Momentum const& incoming(std::size_t i) const;
        Four_Momentum const& outgoing(std::size_t i) const;
        Four_Momentum const& momentum(std::size_t i) const;
        void incoming(std::size_t i, Four_Momentum const& p);
        void outgoing(std::size_t i, Four_Momentum const& p);
    };
}

#endif // Feynumeric_KINEMATICS_HPP