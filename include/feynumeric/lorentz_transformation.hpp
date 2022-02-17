#ifndef FEYNUMERIC_LORENTZ_TRANSFORMATION_HPP
#define FEYNUMERIC_LORENTZ_TRANSFORMATION_HPP

#include "feynumeric/utility.hpp"
#include "feynumeric/matrix.hpp"
#include "feynumeric/momentum.hpp"

namespace Feynumeric
{
    Matrix boost(Four_Momentum const& p, Matrix const& four_vector);
    Matrix boost(Four_Momentum const& p, Four_Momentum const& q);

    Matrix rotateX(Four_Momentum const& p, double cos_theta);
	Matrix rotateX(Matrix const& p, double cos_theta);

	Matrix rotateY(Four_Momentum const& p, double cos_theta);
	Matrix rotateY(Matrix const& p, double cos_theta);

	Matrix rotateZ(Four_Momentum const& p, double cos_theta);
	Matrix rotateZ(Matrix const& p, double cos_theta);
}

#endif // FEYNUMERIC_LORENTZ_TRANSFORMATION_HPP