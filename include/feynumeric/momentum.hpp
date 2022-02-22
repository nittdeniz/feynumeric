
#ifndef Feynumeric_MOMENTUM_HPP
#define Feynumeric_MOMENTUM_HPP

#include "matrix.hpp"
#include "four_vector.hpp"

namespace Feynumeric
{
	double kallen_lambda(double a, double b, double c);

	double momentum(double M, double m1, double m2);

	Four_Vector four_momentum(double q, double m, double cos_theta = 1., double cos_phi = 1.);
}
#endif // Feynumeric_MOMENTUM_HPP