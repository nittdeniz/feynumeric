#ifndef FORM_FACTORS_HPP
#define FORM_FACTORS_HPP

#include <functional>
#include <map>
#include <string>
#include <feynumeric/particle.hpp>

using FORM_FACTOR_FUNCTION = std::function<double(Feynumeric::Particle_Ptr const& R, Feynumeric::Particle_Ptr const& N, Feynumeric::Particle_Ptr const& pi, double virtual_mass)>;

extern FORM_FACTOR_FUNCTION identity;
extern FORM_FACTOR_FUNCTION moniz;
extern FORM_FACTOR_FUNCTION manley;
extern FORM_FACTOR_FUNCTION cassing;
extern FORM_FACTOR_FUNCTION cutkosky;
extern FORM_FACTOR_FUNCTION breit_wigner;
extern FORM_FACTOR_FUNCTION gaussian;
extern FORM_FACTOR_FUNCTION multipol_gauss;

extern std::string const CMD_FORM_FACTOR_NONE;
extern std::string const CMD_FORM_FACTOR_CASSING;
extern std::string const CMD_FORM_FACTOR_CUTKOSKY;
extern std::string const CMD_FORM_FACTOR_MANLEY;
extern std::string const CMD_FORM_FACTOR_MONIZ;
extern std::string const CMD_FORM_FACTOR_BREIT_WIGNER;
extern std::string const CMD_FORM_FACTOR_GAUSSIAN;
extern std::string const CMD_FORM_FACTOR_MULTIPOL_GAUSS;

extern std::map<std::string, FORM_FACTOR_FUNCTION> ff_dict;
#endif