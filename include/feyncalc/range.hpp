#ifndef FEYNCALC_RANGE_HPP
#define FEYNCALC_RANGE_HPP

namespace Feyncalc
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

#endif // FEYNCALC_RANGE_HPP