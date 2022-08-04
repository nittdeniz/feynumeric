#ifndef Feynumeric_INTEGRATE_HPP
#define Feynumeric_INTEGRATE_HPP

#include <functional>

namespace Feynumeric
{
    double integrate(std::function<double(double)> const& f, double const left, double const right, double const epsilon=1.e-6);
	double integrate(std::function<double(double)> const& f, double const left, double const f_left, double const right, double const f_right, double const epsilon =1.e-6);
    double integrate(std::function<double(double)> const& f, double const left, double const f_left, double const mid, double const f_mid, double const right, double const f_right, double const epsilon);
}

#endif // Feynumeric_INTEGRATE_HPP