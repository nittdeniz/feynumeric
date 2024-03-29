#ifndef FEYNUMERIC_LORENTZ_TRANSFORMATION_HPP
#define FEYNUMERIC_LORENTZ_TRANSFORMATION_HPP

#include "feynumeric/utility.hpp"
#include "feynumeric/matrix.hpp"
#include "feynumeric/four_vector.hpp"

namespace Feynumeric
{
    Matrix boost(Four_Vector const& p, Matrix const& four_vector);
    Matrix boost(Four_Vector const& p, Four_Vector const& q);

    Matrix rotateX(double cos_theta);
	Matrix rotateY(double cos_theta);
	Matrix rotateZ(double cos_theta);
}

#endif // FEYNUMERIC_LORENTZ_TRANSFORMATION_HPP