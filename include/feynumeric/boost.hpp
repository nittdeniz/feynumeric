#ifndef FEYNUMERIC_BOOST_HPP
#define FEYNUMERIC_BOOST_HPP

#include "feynumeric/utility.hpp"
#include "feynumeric/matrix.hpp"

namespace Feynumeric
{
    Matrix boost(Matrix const& p, Matrix const& a);
}

#endif // FEYNUMERIC_BOOST_HPP