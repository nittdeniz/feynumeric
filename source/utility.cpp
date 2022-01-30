#include "utility.hpp"

namespace Feynumeric
{
    bool almost_identical(double a, double b, double epsilon)
    {
        return std::abs(a-b) < std::abs(epsilon * std::min(a, b));
    }

    bool almost_identical(Complex aa, Complex bb, double epsilon)
    {
        double a = (aa * std::conj(aa)).real();
        double b = (bb * std::conj(bb)).real();
        return std::abs(a-b) < std::abs(epsilon * std::min(a, b));
    }

    Complex dot3(Matrix const& a, Matrix const& b)
    {
        return a.at(0) * b.at(0) + a.at(1) * b.at(1) + a.at(2) * b.at(2);
    }

    Complex dot4(Matrix const& a, Matrix const& b)
    {
    return a.at(0) * b.at(0) - a.at(1) * b.at(1) - a.at(2) * b.at(2) - a.at(3) * b.at(3);
    }
}

