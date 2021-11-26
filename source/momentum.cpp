#include "feyncalc/momentum.hpp"

namespace Feyncalc
{
    double kallen_lambda(double a, double b, double c)
    {
        return a * a + b * b + c * c - 2 * (a * b + b * c + c * a);
    }

    double momentum(double M, double m1, double m2)
    {
        return std::sqrt(kallen_lambda(M*M, m1*m1, m2*m2)) / (2*M);
    }

    Matrix four_momentum(double mass, double momentum, double cos_theta, double cos_phi)
    {
        const double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
        const double sin_phi   = std::sqrt(1 - cos_phi   * cos_phi);
        return Matrix(4, 1, {std::sqrt(mass*mass + momentum * momentum), momentum * sin_theta * cos_phi, momentum * sin_theta * sin_phi, momentum * cos_theta});
    }
}