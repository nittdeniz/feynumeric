#ifndef Feynumeric_RANGE_HPP
#define Feynumeric_RANGE_HPP

namespace Feynumeric
{
    class Range
    {
    private:
        double _start = 0.
             , _end   = 0.
             , _delta = 0.;
    public:
        Range(double start, double end, double delta = 1.);

    };
}

#endif // Feynumeric_RANGE_HPP