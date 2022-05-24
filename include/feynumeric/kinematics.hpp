#ifndef Feynumeric_KINEMATICS_HPP
#define Feynumeric_KINEMATICS_HPP

#include "four_vector.hpp"
#include "topology.hpp"

#include <vector>

namespace Feynumeric
{
    class Kinematics
    {
    private:
        double _sqrt_s;
        std::size_t _n_in, _n_out;
        std::vector<Four_Vector> _momenta;
    public:
        Kinematics(double sqrt_s, std::size_t n_in, std::size_t n_out);
        double sqrt_s() const;
        Four_Vector const& incoming(std::size_t i) const;
        Four_Vector const& outgoing(std::size_t i) const;
        Four_Vector const& momentum(std::size_t i) const;
        void incoming(std::size_t i, Four_Vector const& p);
        void outgoing(std::size_t i, Four_Vector const& p);
    };
}

#endif // Feynumeric_KINEMATICS_HPP