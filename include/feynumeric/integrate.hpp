#ifndef Feynumeric_INTEGRATE_HPP
#define Feynumeric_INTEGRATE_HPP

#include <functional>

namespace Feynumeric
{
    double integrate(std::function<double(double)> f, double a, double b);

    double integrate(std::function<double(double)> f, double f_a, double f_b, double f_c, double a, double b, double c);
}

#endif // Feynumeric_INTEGRATE_HPP