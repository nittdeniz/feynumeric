#ifndef FEYNCALC_MOMENTUM_HPP
#define FEYNCALC_MOMENTUM_HPP

#include "matrix.hpp"

namespace Feyncalc
{
    double kallen_lambda(double a, double b, double c);
    double momentum(double M, double m1, double m2);
    Matrix four_momentum(double mass, double q, double cos_theta);

}
#endif // FEYNCALC_MOMENTUM_HPP

