#include <cmath>

#include "phase_space.hpp"
#include "units.hpp"


namespace Feynumeric{
    double phase_space2(std::size_t N_polarisations, double sqrt_s, double qout, double qin)
    {
        using namespace Feynumeric::Units;
        return
                1. / N_polarisations * 1. / ( 32 * M_PI * sqrt_s * sqrt_s) * qout / qin * 1._hbarc * 1._hbarc;
    }

    double phase_space3(std::size_t N_polarisations, double sqrt_s, double qout, double qin)
    {
        using namespace Feynumeric::Units;
        return
                1. / N_polarisations * 1. / ( 32 * M_PI * sqrt_s * sqrt_s) * qout / qin * 1._hbarc * 1._hbarc;
    }
}
