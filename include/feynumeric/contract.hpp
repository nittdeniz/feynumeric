#ifndef FEYNUMERIC_CONTRACT_HPP
#define FEYNUMERIC_CONTRACT_HPP

#include <memory>

#define CONTRACT_MATRIX(EXPR, VAR) [&](){ \
    Matrix result_$$123$$(4, 4);                      \
    Lorentz_Index_Ptr VAR = std::make_shared<Lorentz_Index>();\
    result_$$123$$ += EXPR; ++(*VAR);            \
    result_$$123$$ += EXPR; ++(*VAR);            \
    result_$$123$$ += EXPR; ++(*VAR);            \
    result_$$123$$ += EXPR;            \
    return result_$$123$$;}()
#define CONTRACT_SCALAR(EXPR, VAR) [&](){ \
    Complex result_$$123$$;                      \
    Lorentz_Index_Ptr VAR = std::make_shared<Lorentz_Index>();\
    result_$$123$$ += EXPR; ++(*VAR);            \
    result_$$123$$ += EXPR; ++(*VAR);            \
    result_$$123$$ += EXPR; ++(*VAR);            \
    result_$$123$$ += EXPR;            \
    return result_$$123$$;}()
#endif