#ifndef FEYNCALC_ANGULAR_MOMENTUM_HPP
#define FEYNCALC_ANGULAR_MOMENTUM_HPP

namespace Feyncalc
{
    class Angular_Momentum
    {
    private:
        double _value = 0.;
    public:
        Angular_Momentum(double value = 0.);
        operator double() const;
        static bool is_valid_spin(double value);
    };

}

#endif // FEYNCALC_ANGULAR_MOMENTUM_HPP