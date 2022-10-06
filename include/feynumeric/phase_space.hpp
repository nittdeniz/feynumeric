#ifndef FEYNUMERIC_PHASE_SPACE_HPP
#define FEYNUMERIC_PHASE_SPACE_HPP

namespace Feynumeric{
    double phase_space2(std::size_t N_polarisations, double sqrt_s, double qout, double qin);
    double phase_space3(double s, double qout, double qin, double qout_star);
}

#endif //FEYNUMERIC_PHASE_SPACE_HPP
