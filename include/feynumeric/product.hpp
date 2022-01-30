#ifndef FEYNUMERIC_PRODUCT_HPP
#define FEYNUMERIC_PRODUCT_HPP

namespace Feynumeric
{
    template<typename T>
    inline T product(std::function<T(int)> f, int const start, int const end, int const step_size = 1)
    {
        T result = 1.;
        for( int i = start; i <= end; i += step_size ){
            result *= f(i);
        }
        return result;
    }
}
#endif // FEYNUMERIC_PRODUCT_HPP