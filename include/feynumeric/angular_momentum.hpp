#ifndef Feynumeric_ANGULAR_MOMENTUM_HPP
#define Feynumeric_ANGULAR_MOMENTUM_HPP

#include <memory>

namespace Feynumeric
{
    class Angular_Momentum
    {
    private:
        double _j = 0.;
        double _m = 0.;
        bool _massless;
    public:
        Angular_Momentum(double j = 0., double m = 0., bool massless = false);
        Angular_Momentum(Angular_Momentum const& J);
        Angular_Momentum& operator=(Angular_Momentum const& J);

        static bool is_valid_spin(double j);

        bool is_half_odd_integer() const;
        bool massless() const;

        double j() const;
        double m() const;
        void m(double new_m);

        void reset();

        std::size_t n_states() const;
        friend Angular_Momentum& operator++(Angular_Momentum& J);
    };

	Angular_Momentum& operator++(Angular_Momentum& J);

    using Angular_Momentum_Ptr = std::shared_ptr<Angular_Momentum>;

}

#endif // Feynumeric_ANGULAR_MOMENTUM_HPP