
#ifndef Feynumeric_MOMENTUM_HPP
#define Feynumeric_MOMENTUM_HPP

#include "matrix.hpp"
#include "four_vector.hpp"

namespace Feynumeric
{
	double kallen_lambda(double a, double b, double c);

	double momentum(double M, double m1, double m2);

	Four_Vector four_momentum(double q, double m, double cos_theta, double cos_phi);
}
    //Matrix four_momentum(double mass, double q, double cos_theta);
/*
    class Three_Momentum : public Matrix
    {
    public:
        Three_Momentum();
        Three_Momentum(double x, double y, double z);
        Three_Momentum(Three_Momentum const& other);
        Three_Momentum& operator=(Three_Momentum const& other);

        double x() const;
        double y() const;
        double z() const;
        void x(double x);
        void y(double y);
        void z(double z);

        double dot() const;

        friend Three_Momentum operator+(Three_Momentum const& lhs, Three_Momentum const& rhs);
	    friend Three_Momentum operator-(Three_Momentum const& lhs, Three_Momentum const& rhs);
        template<typename T>
        friend Three_Momentum operator/(Three_Momentum const& lhs, T const& rhs);
	    template<typename T>
	    friend Three_Momentum operator*(Three_Momentum const& lhs, T const& rhs);
	    template<typename T>
	    friend Three_Momentum operator*(T const& lhs, Three_Momentum const& rhs);
    };

	Three_Momentum operator+(Three_Momentum const& lhs, Three_Momentum const& rhs);
	Three_Momentum operator-(Three_Momentum const& lhs, Three_Momentum const& rhs);

	template<typename T>
	Three_Momentum operator*(T const& lhs, Three_Momentum const& rhs)
	{
		return rhs * lhs;
	}

    template<typename T>
    Three_Momentum operator*(Three_Momentum const& lhs, T const& rhs)
    {
    	Three_Momentum copy(lhs);
    	copy._data[0] *= rhs;
    	copy._data[1] *= rhs;
    	copy._data[2] *= rhs;
    	return copy;
    }

    template<typename T>
    Three_Momentum operator/(Three_Momentum const& lhs, T const& rhs)
    {
		auto inv = 1./rhs;
		return lhs * inv;
    }

    class Four_Momentum : public Matrix
    {
    public:
        Four_Momentum();
        Four_Momentum(double q, double m, double cos_theta=1, double cos_phi=1);
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

        double squared() const;

        double dot() const;

        friend Four_Momentum operator+(Four_Momentum const& lhs, Four_Momentum const& rhs);
        friend Four_Momentum operator-(Four_Momentum const& lhs);
	    friend Four_Momentum operator-(Four_Momentum const& lhs, Four_Momentum const& rhs);
	    template<typename T>
	    friend Four_Momentum operator*(T const& lhs, Four_Momentum const& rhs);
	    template<typename T>
	    friend Four_Momentum operator*(Four_Momentum const& lhs, T const& rhs);

	    template<typename T>
	    friend Four_Momentum operator/(Four_Momentum const& lhs, T const& rhs);
    };

	Four_Momentum operator+(Four_Momentum const& lhs, Four_Momentum const& rhs);
	Four_Momentum operator-(Four_Momentum const& lhs);
	Four_Momentum operator-(Four_Momentum const& lhs, Four_Momentum const& rhs);
	template<typename T>
	Four_Momentum operator*(T const& lhs, Four_Momentum const& rhs){
		Four_Momentum copy(rhs);
		copy._data[0] *= lhs;
		copy._data[1] *= lhs;
		copy._data[2] *= lhs;
		copy._data[3] *= lhs;
		return copy;
	}
	template<typename T>
	Four_Momentum operator*(Four_Momentum const& lhs, T const& rhs)
	{
		return rhs * lhs;
	}

	template<typename T>
	Four_Momentum operator/(Four_Momentum const& lhs, T const& rhs)
	{
		Four_Momentum copy(lhs);
		copy._data[0] /= rhs;
		copy._data[1] /= rhs;
		copy._data[2] /= rhs;
		copy._data[3] /= rhs;
		return copy;
	}
}
*/
#endif // Feynumeric_MOMENTUM_HPP