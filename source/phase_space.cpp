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

    double phase_space3(double s, double qout, double qin, double qout_star)
    {
        using namespace Feynumeric::Units;
        return std::pow(2*M_PI, -5.) * 1./(32 * s) * qout/qin * qout_star * 1._hbarc * 1._hbarc;
    }
}
