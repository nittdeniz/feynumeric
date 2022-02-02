#ifndef Feynumeric_ANGULAR_MOMENTUM_HPP
#define Feynumeric_ANGULAR_MOMENTUM_HPP

namespace Feynumeric
{
    class Angular_Momentum
    {
    private:
        double _j = 0.;
        double _m = 0.;
    public:
        Angular_Momentum(double j = 0., double m = 0.);
        Angular_Momentum(Angular_Momentum const& J);
        Angular_Momentum& operator=(Angular_Momentum const& J);
        static bool is_valid_spin(double j);
        bool is_half_odd_integer() const;

        double j() const;
        double m() const;
        void m(double new_m);
    };

}

#endif // Feynumeric_ANGULAR_MOMENTUM_HPP