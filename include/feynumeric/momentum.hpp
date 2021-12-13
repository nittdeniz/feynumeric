#ifndef Feynumeric_MOMENTUM_HPP
#define Feynumeric_MOMENTUM_HPP

#include "matrix.hpp"

namespace Feynumeric
{
    double kallen_lambda(double a, double b, double c);
    double momentum(double M, double m1, double m2);
    Matrix four_momentum(double mass, double q, double cos_theta);

    class Three_Momentum : public Matrix
    {
    public:
        Three_Momentum();
        Three_Momentum(Three_Momentum const& other);
        Three_Momentum& operator=(Three_Momentum const& other);

        double x() const;
        double y() const;
        double z() const;
        void x(double x);
        void y(double y);
        void z(double z);

        double dot() const;
    };

    class Four_Momentum : public Matrix
    {
    public:
        Four_Momentum();
        Four_Momentum(Four_Momentum const& other);
        Four_Momentum& operator=(Four_Momentum const& other);

        double E() const;
        double x() const;
        double y() const;
        double z() const;

        void E(double E);
        void x(double x);
        void y(double y);
        void z(double z);

        double dot() const;
    };

}
#endif // Feynumeric_MOMENTUM_HPP

