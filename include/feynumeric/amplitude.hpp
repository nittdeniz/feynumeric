#ifndef FEYNUMERIC_AMPLITUDE_HPP
#define FEYNUMERIC_AMPLITUDE_HPP

#include "functions.hpp"
#include "polynomial.hpp"

#include <vector>

namespace Feynumeric
{
	template<std::size_t N>
	class Amplitude
	{
		std::vector<FPolynomial<N>> _fpolynomial;
		Amplitude(std::vector<std::vector<Polynomial>> const& polynomials)
		{

		}
	};
}

#endif