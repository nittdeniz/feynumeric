#ifndef Feynumeric_LORENTZ_INDEX_HPP
#define Feynumeric_LORENTZ_INDEX_HPP

namespace Feynumeric
{
    class Lorentz_Index
    {
    private:
		unsigned int _value;
    public:
    	operator unsigned int() const;
    	Lorentz_Index& operator++();
    	Lorentz_Index& operator--();
    	void reset();
    };
}

#endif // Feynumeric_LORENTZ_INDEX_HPP