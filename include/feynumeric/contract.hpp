#ifndef FEYNUMERIC_CONTRACT_HPP
#define FEYNUMERIC_CONTRACT_HPP

#include <memory>

#define CONTRACT_MATRIX(EXPR, VAR) [&](){ \
    Matrix result(4, 4);                      \
    Lorentz_Index_Ptr VAR = std::make_shared<Lorentz_Index>();\
    result += EXPR; ++(*VAR);            \
    result += EXPR; ++(*VAR);            \
    result += EXPR; ++(*VAR);            \
    result += EXPR;            \
    return result;}()
#define CONTRACT_SCALAR(EXPR, VAR) [&](){ \
    Complex result;                      \
    Lorentz_Index_Ptr VAR = std::make_shared<Lorentz_Index>();\
    result += EXPR; ++(*VAR);            \
    result += EXPR; ++(*VAR);            \
    result += EXPR; ++(*VAR);            \
    result += EXPR;            \
    return result;}()
#endif