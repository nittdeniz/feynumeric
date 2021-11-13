#include "integrate.hpp"

namespace Feyncalc
{
    double integrate(std::function<double(double)> f, double a, double b)
    {
        return f(b) - f(a);
    }

    double integrate(std::function<double(double)> f, double f_a, double f_b, double f_c, double a, double b, double c)
    {
        return a+b+c+f_a+f_b+f_c+f(a);
    }
}