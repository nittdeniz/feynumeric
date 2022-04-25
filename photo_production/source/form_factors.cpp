#include "form_factors.hpp"

FORM_FACTOR_FUNCTION identity = [](double, double, double, double){
	return 1;
};

FORM_FACTOR_FUNCTION moniz = [](double beta, double q0, double q, int l){
	return 1;
};

FORM_FACTOR_FUNCTION manley = [](double beta, double q0, double q, int l){
	return 1;
};

FORM_FACTOR_FUNCTION cassing = [](double beta, double q0, double q, int l){
	return 1;
};

FORM_FACTOR_FUNCTION cutkosky = [](double, double q0, double q, int l){
	return 1;
};

