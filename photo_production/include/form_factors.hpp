#ifndef FORM_FACTORS_HPP
#define FORM_FACTORS_HPP

#include <functional>
#include <feynumeric/particle.hpp>

using FORM_FACTOR_FUNCTION = std::function<double(Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double virtual_mass)>;

extern FORM_FACTOR_FUNCTION identity;
extern FORM_FACTOR_FUNCTION moniz;
extern FORM_FACTOR_FUNCTION manley;
extern FORM_FACTOR_FUNCTION cassing;
extern FORM_FACTOR_FUNCTION cutkosky;
extern FORM_FACTOR_FUNCTION breit_wigner;
extern FORM_FACTOR_FUNCTION dyson_factor_32;
extern FORM_FACTOR_FUNCTION dyson_factor_12;
#endif