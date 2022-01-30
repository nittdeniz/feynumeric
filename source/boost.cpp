#include "feynumeric/boost.hpp"

namespace Feynumeric
{
    Matrix boost(Matrix const& p, Matrix const& a)
    {
        if(
                ( (p.n_cols() == 1 && p.n_rows() == 4) || (p.n_cols() == 4 && p.n_rows() == 1) )
                &&  ( (a.n_cols() == 1 && a.n_rows() == 4) || (a.n_cols() == 4 && a.n_rows() == 1) )
                )
        {
            Matrix beta(1, 3, {p.at(1), p.at(2), p.at(3)});
            beta /= p.at(0);
            double beta_squared = dot3(beta, beta).real();
            double gamma = 1./std::sqrt(1-beta_squared);

            Matrix boost(4, 4);
            boost(0, 0) = gamma;
            boost(0, 1) = -gamma * beta.at(0);
            boost(0, 2) = -gamma * beta.at(1);
            boost(0, 3) = -gamma * beta.at(2);
            boost(1, 1) = 1. + (gamma - 1) * beta.at(0) * beta.at(0);
            boost(1, 0) = -gamma * beta.at(0);
            boost(2, 0) = -gamma * beta.at(1);
            boost(3, 0) = -gamma * beta.at(2);

            if( a.n_cols() == 1 )
            {
                return boost * a;
            }
            return boost * a.T();
        }
        critical_error("Invalid boost vector.\n");
    }
}

