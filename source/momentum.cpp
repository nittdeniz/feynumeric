#include "feynumeric/momentum.hpp"
#include "three_vector.hpp"

namespace Feynumeric
{
    double kallen_lambda(double a, double b, double c)
    {
        return a * a + b * b + c * c - 2 * (a * b + b * c + c * a);
    }

    double momentum(double M, double m1, double m2)
    {
        return std::sqrt(kallen_lambda(M*M, m1*m1, m2*m2)) / (2*M);
    }

	Four_Vector four_momentum(double q, double m, double cos_theta, double cos_phi)
	{
		double E = std::sqrt(q*q + m*m);
		auto p = Three_Vector::from_spherical(q, cos_theta, cos_phi);
		return Four_Vector(E, p);
	}
}
