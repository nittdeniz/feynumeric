#ifndef Feynumeric_COMPLEX_HPP
#define Feynumeric_COMPLEX_HPP

#include <complex>
#include <iostream>

namespace Feynumeric
{
    using Complex = std::complex<double>;
    using namespace std::literals::complex_literals;

    inline std::ostream& operator<<(std::ostream& out, Complex const& c){
    	if( c.imag() == 0 ){
    		return out << c.real();
    	}
    	if( c.real() == 0 ){
			return out << c.imag() << "I";
    	}
    	if( c.imag() < 0 ){
    		return out << c.real() << c.imag() << "I";
    	}
    	return out << "(" << c.real() << "+" << c.imag() << "I)";
    }
}

#endif // Feynumeric_COMPLEX_HPP