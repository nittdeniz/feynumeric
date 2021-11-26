#ifndef FEYNCALC_DIRAC_HPP
#define FEYNCALC_DIRAC_HPP

#include <array>

#include "angular_momentum.hpp"
#include "matrix.hpp"

namespace Feyncalc
{
    extern std::array<Matrix, 4> GA;
    extern std::array<Matrix, 4> GAC;

    Matrix GS(Matrix const& a);

    Matrix dirac_sigma(Matrix const& a, Matrix const& b);

    Matrix epsilon1(Matrix const& p, Angular_Momentum const& lambda);
    Matrix epsilon2(Matrix const& p, Angular_Momentum const& lambda);

    Matrix u12(Matrix const& p, Angular_Momentum const& s);
    Matrix ubar12(Matrix const& p, Angular_Momentum const& s);

    Matrix u32(Matrix const& p, Angular_Momentum const& s);
    Matrix ubar32(Matrix const& p, Angular_Momentum const& s);

    Matrix u52(Matrix const& p, Angular_Momentum const& s);
    Matrix ubar52(Matrix const& p, Angular_Momentum const& s);

    Matrix v12(Matrix const& p, Angular_Momentum const& s);
    Matrix vbar12(Matrix const& p, Angular_Momentum const& s);

    Matrix v32(Matrix const& p, Angular_Momentum const& s);
    Matrix vbar32(Matrix const& p, Angular_Momentum const& s);

    Matrix v52(Matrix const& p, Angular_Momentum const& s);
    Matrix vbar52(Matrix const& p, Angular_Momentum const& s);
}

#endif // FEYNCALC_DIRAC_HPP


