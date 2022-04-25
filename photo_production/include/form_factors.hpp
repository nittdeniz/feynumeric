#ifndef FORM_FACTORS_HPP
#define FORM_FACTORS_HPP

#include <functional>

using FORM_FACTOR_FUNCTION = std::function<double(double, double, double, int)>;

extern FORM_FACTOR_FUNCTION identity;
extern FORM_FACTOR_FUNCTION moniz;
extern FORM_FACTOR_FUNCTION  manley;
extern FORM_FACTOR_FUNCTION cassing;
extern FORM_FACTOR_FUNCTION cutkosky;


#endif