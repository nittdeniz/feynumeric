#ifndef FEYNUMERIC_SUM_HPP
#define FEYNUMERIC_SUM_HPP

#include <functional>

namespace Feynumeric
{
    template<typename T>
    inline T sum(std::function<T(int)> f, int const start, int const end, int const step_size = 1)
    {
        T result = 0.;
        for( int i = start; i <= end; i += step_size ){
            result += f(i);
        }
        return result;
    }
}

#endif // FEYNUMERIC_SUM_HPP