#include "lorentz_index.hpp"

#include <stdexcept>

namespace Feynumeric
{
	Lorentz_Index::operator unsigned int() const{
		return _value;
	}

	Lorentz_Index& Lorentz_Index::operator++(){
		_value = (_value+1)%4;
		return *this;
	}

	Lorentz_Index& Lorentz_Index::operator--(){
		_value = (_value-1)%4;
		return *this;
	}

	void Lorentz_Index::reset(){
		_value = 0;
	}
}