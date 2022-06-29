#ifndef Feynumeric_COMPLEX_HPP
#define Feynumeric_COMPLEX_HPP

#include "format.hpp"

#include <complex>
#include <iostream>

namespace Feynumeric
{
    using Complex = std::complex<double>;
    using namespace std::literals::complex_literals;

    inline std::ostream& operator<<(std::ostream& out, Complex const& c){
    	if( c.imag() == 0 ){
    		if( c.real() < 0 ){
    			return out << FORMAT("({:.16f})", c.real());
    		}
    		return out << FORMAT("{:.16f}", c.real());
    	}
    	if( c.real() == 0 ){
    		if( c.imag() < 0 ){
    			return out << FORMAT("({:.16f}I)", c.imag());
    		}
			return out << FORMAT("{:.16f}I", c.imag());
    	}
    	if( c.imag() < 0 ){
    		return out << FORMAT("({:.16f}{:.16f}I)", c.real(), c.imag());
    	}
    	return out << FORMAT("({:.16f}+{:.16f}I)", c.real(), c.imag());
    }
}

#endif // Feynumeric_COMPLEX_HPP